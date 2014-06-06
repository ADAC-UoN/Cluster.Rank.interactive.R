Cluster.Rank.interactive.R
==========================

Cluster Rank - compare and rank hierarchical clustering methods

## Harry Clifford and Richard D. Emes
## School of Veterinary Medicine and Science, University of Nottingham. UK. 2011
## Please cite Clifford et al Frontiers in Genetics 2011

# http://journal.frontiersin.org/Journal/10.3389/fgene.2011.00088/abstract

Developed for illumina Methylation arrays (27k & 450k arrays) but will work with any tabulated numerical data 


Comparison of clustering methods for investigation of genome-wide methylation array data
Harry Clifford1, Frank Wessely1, Satish Pendurthi1,2 and Richard D. Emes1*

    1 School of Veterinary Medicine and Science, University of Nottingham, Nottingham, UK
    2 School of Contemporary Studies, University of Abertay, Dundee, UK

The use of genome-wide methylation arrays has proved very informative to investigate both clinical and biological questions in human epigenomics. The use of clustering methods either for exploration of these data or to compare to an a priori grouping, e.g., normal versus disease allows assessment of groupings of data without user bias. However no consensus on the methods to use for clustering of methylation array approaches has been reached. To determine the most appropriate clustering method for analysis of illumina array methylation data, a collection of data sets was simulated and used to compare clustering methods. Both hierarchical clustering and non-hierarchical clustering methods (k-means, k-medoids, and fuzzy clustering algorithms) were compared using a range of distance and linkage methods. As no single method consistently outperformed others across different simulations, we propose a method to capture the best clustering outcome based on an additional measure, the silhouette width. This approach produced a consistently higher cluster accuracy compared to using any one method in isolation.
