## Harry Clifford and Richard D. Emes
## School of Veterinary Medicine and Science, University of Nottingham. UK. 2011
## Please cite Clifford et al Frontiers in Genetics 2011

# http://journal.frontiersin.org/Journal/10.3389/fgene.2011.00088/abstract


##This script takes a table of beta-values obtained from infinium assays in .csv format, and clusters the assay results based on similarities
##Hierarchical and fuzzyclustering are the methods used, and within these, various distance and linkage methods are used
##The total ouput will include 32 clustering results, and a second output with the best 5 (scored by silhouette width) will also be returned
##The user will be asked to enter the number of clusters required in the output, and locations for input and output data
##The script also asks whether to include P-values in dendrograms. This runs with pvclust, which may take a much longer time
##This script will then write the ouput plots in a pdf file

install.packages("clv")
install.packages("foreach")
install.packages("iterators")
install.packages("pvclust")
##selects data for use, number of clusters to calculate, and output pathway
winDialog("ok", "Please select your data")
data <- read.csv(file.choose(), header=TRUE,row.names=1) #the csv file is assumed to have row AND column names
clusternumber <- as.numeric(winDialogString("How many clusters would you like to highlight?","2"))
winDialog("ok", "Please select the location for your output")
output <- choose.dir()
FullOutput <- paste(output,"AllResults.pdf", sep="\\")
Top5Output <- paste(output,"Top5Results.pdf", sep="\\")
##asks the user whether to include p-values in results (by using pvclust)
Top5Results.P <- winDialog("yesno", "Would you like to include P values in your Top 5 Results?")
AllResults.P <- winDialog("yesno", "Would you like to include P values in your Total Results?")
##transposes data to make it appropriate for use, then loads several of the packages
transposed.data <- t(data)
library(foreach)
library(cluster)
library(pvclust)
##uses looping to generate 8 different distance matrices, using different distance measures
distances <- foreach(a=rep(c("euclidean","maximum","manhattan","canberra"),8)) %do% dist(transposed.data,method=a,diag=TRUE,upper=TRUE)
##generates character vectors containing linkage and distance methods to be used with pvclust (should pvclust actually need to be run)
pvclust.linkages <- c(rep(c("ward","single","complete","average","mcquitty","median","centroid"),each=4))
pvclust.distances <- c(rep(c("euclidean","maximum","manhattan","canberra"),7))
##uses looping to run hierarchical clustering and fuzzy clustering with varying methods
hclust.clusters <- foreach(a=rep(c("ward","single","complete","average","mcquitty","median","centroid"),each=4),
Cluster=distances) %do% hclust(Cluster,method=a)
fuzzyclust.clusters <- foreach(a=distances[1:4]) %do% fanny(a,k=clusternumber,memb.exp=1.1)
##if the user said yes for obtaining p-values in AllResults, uses looping to run pvclusts, otherwise keeps hierarchical results instead
if(AllResults.P == "YES"){
pvclust.clusters <- foreach(a=pvclust.linkages,b=pvclust.distances) %do% pvclust(data,method.hclust=a,method.dist=b)
clusters <- c(pvclust.clusters,fuzzyclust.clusters)}else{
clusters <- c(hclust.clusters,fuzzyclust.clusters)}
rect.clusters <- c(hclust.clusters,fuzzyclust.clusters)
##uses looping to find the cluster vectors
clustvectors <- c(foreach(a=hclust.clusters) %do% cutree(a,k=clusternumber),foreach(b=fuzzyclust.clusters) %do% b$clustering)
##applies cluster vectors to find average silhouette widths, then enters these into a matrix
silwidths <- foreach(a=clustvectors,b=distances) %do% summary(silhouette(a,b))$avg.width
sil.values<-matrix(silwidths,ncol=32)
colnames(sil.values) <- c(1:32)
##orders the clusters by silhouette width then uses looping to find which are the best 5 clusters
top5.clusters <- invisible(foreach(a=1:5) %do% eval(parse(text=(names(sil.values[,])[order(unlist(sil.values),decreasing=TRUE)[a]]))))
##creates a character vector of plot titles in the order they will be created
plot.titles <- c("Hierarchical Ward Euclidean","Hierarchical Ward Maximum","Hierarchical Ward Manhattan","Hierarchical Ward Canberra",
"Hierarchical Single Euclidean","Hierarchical Single Maximum","Hierarchical Single Manhattan","Hierarchical Single Canberra",
"Hierarchical Complete Euclidean","Hierarchical Complete Maximum","Hierarchical Complete Manhattan","Hierarchical Complete Canberra",
"Hierarchical Average Euclidean","Hierarchical Average Maximum","Hierarchical Average Manhattan","Hierarchical Average Canberra",
"Hierarchical Mcquitty Euclidean","Hierarchical Mcquitty Maximum","Hierarchical Mcquitty Manhattan","Hierarchical Mcquitty Canberra",
"Hierarchical Median Euclidean","Hierarchical Median Maximum","Hierarchical Median Manhattan","Hierarchical Median Canberra",
"Hierarchical Centroid Euclidean","Hierarchical Centroid Maximum","Hierarchical Centroid Manhattan","Hierarchical Centroid Canberra",
"FuzzyCluster Euclidean","FuzzyCluster Maximum","FuzzyCluster Manhattan","FuzzyCluster Canberra")
##runs the best 5 clusters again with fuzzyclust, hierarchical clust or pvclust (depending on which is appropriate) by looping
top5.plots <- foreach (a=top5.clusters) %do% {
if(a>28){
fuzzyclust.clusters[[a-28]]}else{
if(Top5Results.P == "YES"){
pvclust(data,method.hclust=pvclust.linkages[a],method.dist=pvclust.distances[a])}else{
hclust.clusters[[a]]}}}
##builds a pdf of all results with looping, and adds titles and labels (where appropriate)
pdf(FullOutput)
invisible(foreach(a=clusters,b=plot.titles,c=silwidths,d=rect.clusters,e=c(1:32)) %do% {
plot(a, main="")
title(b, line=2, font.main=2)
title("Silhouette Width:", line=1, font.main=10)
title(c, line=0, font.main=10)
if(e<29){
rect.hclust(d,k=clusternumber)}}) #this adds rectangles around clusters in dendrograms only
invisible(dev.off())
##builds another pdf of the top 5 results with looping, once again adding titles and labels where appropriate
pdf(Top5Output)
invisible(foreach(a=1:5,b=top5.clusters,c=top5.plots) %do% {
plot(c,main="")
title(a, line=3, font.main=2)
title("Average Silhouette Width:", line=1, font.main=10)
title((sort(unlist(sil.values), decreasing=TRUE)[a]), line=0, font.main=10)
colnames(sil.values) <- c(1:32)
title(plot.titles[[b]], line=2, font.main=2)
if(b<29){
rect.hclust(rect.clusters[[b]],k=clusternumber)}})
invisible(dev.off())