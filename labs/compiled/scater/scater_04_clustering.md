---
author: "Åsa Björklund  &  Paulo Czarnewski"
date: "Sept 13, 2019"
output:
  html_document:
    self_contained: true
    highlight: tango
    df_print: paged
    toc: yes
    toc_float:
      collapsed: false
      smooth_scroll: true
    toc_depth: 3
    keep_md: yes
    fig_caption: true
  html_notebook:
    self_contained: true
    highlight: tango
    df_print: paged
    toc: yes
    toc_float:
      collapsed: false
      smooth_scroll: true
    toc_depth: 3
---


<style>
h1, .h1, h2, .h2, h3, .h3, h4, .h4 { margin-top: 50px }
p.caption {font-size: 0.9em;font-style: italic;color: grey;margin-right: 10%;margin-left: 10%;text-align: justify}
</style>

# Clustering


In this tutorial we will continue the analysis of the integrated dataset. We will use the integrated PCA to perform the clustering. First we will construct a $k$-nearest neighbour graph in order to perform a clustering on the graph. We will also show how to perform hierarchical clustering and k-means clustering on PCA space.

Let's first load all necessary libraries and also the integrated dataset from the previous step.


```r
suppressPackageStartupMessages({
    library(scater)
    library(scran)
    library(cowplot)
    library(ggplot2)
    library(rafalib)
    library(pheatmap)
    library(igraph)
})

sce <- readRDS("data/3pbmc_qc_dr_int.rds")
```

## Graph clustering
***

The procedure of clustering on a Graph can be generalized as 3 main steps:

1) Build a kNN graph from the data

2) Prune spurious connections from kNN graph (optional step). This is a SNN graph.

3) Find groups of cells that maximizes the connections within the group compared other groups.

### Building kNN / SNN graph


The first step into graph clustering is to construct a k-nn graph, in case you don't have one. For this, we will use the PCA space. Thus, as done for dimensionality reduction, we will use ony the top *N* PCA dimensions for this purpose (the same used for computing UMAP / tSNE).


```r
# These 2 lines are for demonstration purposes only
g <- buildKNNGraph(sce, k = 30, use.dimred = "MNN", assay.type = "RNA")
sce@reducedDims$KNN <- igraph::as_adjacency_matrix(g)

# These 2 lines are the most recommended
g <- buildSNNGraph(sce, k = 30, use.dimred = "MNN", assay.type = "RNA")
sce@reducedDims$SNN <- as_adjacency_matrix(g, attr = "weight")
```

We can take a look at the kNN graph. It is a matrix where every connection between cells is represented as $1$s. This is called a **unweighted** graph (default in Seurat). Some cell connections can however have more importance than others, in that case the scale of the graph from $0$ to a maximum distance. Usually, the smaller the distance, the closer two points are, and stronger is their connection. This is called a **weighted** graph. Both weighted and unweighted graphs are suitable for clustering, but clustering on unweighted graphs is faster for large datasets (> 100k cells).


```r
# plot the KNN graph
pheatmap(sce@reducedDims$KNN[1:200, 1:200], col = c("white", "black"), border_color = "grey90", 
    legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2)
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
# or the SNN graph
pheatmap(sce@reducedDims$SNN[1:200, 1:200], col = colorRampPalette(c("white", "yellow", 
    "red", "black"))(20), border_color = "grey90", legend = T, cluster_rows = F, 
    cluster_cols = F, fontsize = 2)
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

As you can see, the way Scran computes the SNN graph is different to Seurat. It gives edges to all cells that shares a neighbor, but weights the edges by how similar the neighbors are. Hence, the SNN graph has more edges than the KNN graph.


### Clustering on a graph


Once the graph is built, we can now perform graph clustering. The clustering is done respective to a resolution which can be interpreted as how coarse you want your cluster to be. Higher resolution means higher number of clusters.



```r
g <- buildSNNGraph(sce, k = 5, use.dimred = "MNN", assay.type = "RNA")
sce$louvain_SNNk5 <- factor(cluster_louvain(g)$membership)

g <- buildSNNGraph(sce, k = 10, use.dimred = "MNN", assay.type = "RNA", )
sce$louvain_SNNk10 <- factor(cluster_louvain(g)$membership)

g <- buildSNNGraph(sce, k = 15, use.dimred = "MNN", assay.type = "RNA")
sce$louvain_SNNk15 <- factor(cluster_louvain(g)$membership)

plot_grid(ncol = 3, plotReducedDim(sce, use_dimred = "UMAP_on_MNN", colour_by = "louvain_SNNk5", 
    add_ticks = F) + ggplot2::ggtitle(label = "louvain_SNNk5"), plotReducedDim(sce, 
    use_dimred = "UMAP_on_MNN", colour_by = "louvain_SNNk10", add_ticks = F) + ggplot2::ggtitle(label = "louvain_SNNk10"), 
    plotReducedDim(sce, use_dimred = "UMAP_on_MNN", colour_by = "louvain_SNNk15", 
        add_ticks = F) + ggplot2::ggtitle(label = "louvain_SNNk15"))
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


## K-means clustering
***

K-means is a generic clustering algorithm that has been used in many application areas. In R, it can be applied via the kmeans function. Typically, it is applied to a reduced dimension representation of the expression data (most often PCA, because of the interpretability of the low-dimensional distances). We need to define the number of clusters in advance. Since the results depend on the initialization of the cluster centers, it is typically recommended to run K-means with multiple starting configurations (via the nstart argument).


```r
sce$kmeans_5 <- factor(kmeans(x = sce@reducedDims$MNN, centers = 5)$cluster)
sce$kmeans_10 <- factor(kmeans(x = sce@reducedDims$MNN, centers = 10)$cluster)
sce$kmeans_15 <- factor(kmeans(x = sce@reducedDims$MNN, centers = 15)$cluster)

plot_grid(ncol = 3, plotReducedDim(sce, use_dimred = "UMAP_on_MNN", colour_by = "kmeans_5", 
    add_ticks = F) + ggplot2::ggtitle(label = "KMeans5"), plotReducedDim(sce, use_dimred = "UMAP_on_MNN", 
    colour_by = "kmeans_10", add_ticks = F) + ggplot2::ggtitle(label = "KMeans10"), 
    plotReducedDim(sce, use_dimred = "UMAP_on_MNN", colour_by = "kmeans_15", add_ticks = F) + 
        ggplot2::ggtitle(label = "KMeans15"))
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


## Hierarchical clustering
***

### Defining distance between cells

The base R `stats` package already contains a function `dist` that calculates distances between all pairs of samples. Since we want to compute distances between samples, rather than among genes, we need to transpose the data before applying it to the `dist` function. This can be done by simply adding the transpose function `t()` to the data. The distance methods available  in `dist` are: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".


```r
d <- dist(sce@reducedDims$MNN, method = "euclidean")
```

As you might have realized, correlation is not a method implemented in the `dist` function. However, we can create our own distances and transform them to a distance object. We can first compute sample correlations using the `cor` function.
As you already know, correlation range from -1 to 1, where 1 indicates that two samples are closest, -1 indicates that two samples are the furthest and 0 is somewhat in between. This, however, creates a problem in defining distances because a distance of 0 indicates that two samples are closest, 1 indicates that two samples are the furthest and distance of -1 is not meaningful. We thus need to transform the correlations to a positive scale (a.k.a. **adjacency**):

\[adj = \frac{1- cor}{2}\]

Once we transformed the correlations to a 0-1 scale, we can simply convert it to a distance object using `as.dist` function. The transformation does not need to have a maximum of 1, but it is more intuitive to have it at 1, rather than at any other number.


```r
# Compute sample correlations
sample_cor <- cor(Matrix::t(sce@reducedDims$MNN))

# Transform the scale from correlations
sample_cor <- (1 - sample_cor)/2

# Convert it to a distance object
d2 <- as.dist(sample_cor)
```

### Clustering cells

After having calculated the distances between samples calculated, we can now proceed with the hierarchical clustering per-se. We will use the function `hclust` for this purpose, in which we can simply run it with the distance objects created above. The methods available are: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid". It is possible to plot the dendrogram for all cells, but this is very time consuming and we will omit for this tutorial.


```r
# euclidean
h_euclidean <- hclust(d, method = "ward.D2")

# correlation
h_correlation <- hclust(d2, method = "ward.D2")
```

 Once your dendrogram is created, the next step is to define which samples belong to a particular cluster. After identifying the dendrogram, we can now literally cut the tree at a fixed threshold (with `cutree`) at different levels to define the clusters. We can either define the number of clusters or decide on a height. We can simply try different clustering levels.


```r
#euclidean distance
sce$hc_euclidean_5 <- factor( cutree(h_euclidean,k = 5) )
sce$hc_euclidean_10 <- factor( cutree(h_euclidean,k = 10) )
sce$hc_euclidean_15 <- factor( cutree(h_euclidean,k = 15) )

#correlation distance
sce$hc_corelation_5 <- factor( cutree(h_correlation,k = 5) )
sce$hc_corelation_10 <- factor( cutree(h_correlation,k = 10) )
sce$hc_corelation_15 <- factor( cutree(h_correlation,k = 15) )


plot_grid(ncol = 3,
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "hc_euclidean_5",add_ticks = F)+
    ggplot2::ggtitle(label ="HC_euclidean_5"),
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "hc_euclidean_10",add_ticks = F)+
    ggplot2::ggtitle(label ="HC_euclidean_10"),
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "hc_euclidean_15",add_ticks = F)+
    ggplot2::ggtitle(label ="HC_euclidean_15"),
  
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "hc_corelation_5",add_ticks = F)+
    ggplot2::ggtitle(label ="HC_correlation_5"),
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "hc_corelation_10",add_ticks = F)+
    ggplot2::ggtitle(label ="HC_correlation_10"),
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "hc_corelation_15",add_ticks = F)+
    ggplot2::ggtitle(label ="HC_correlation_15")
)
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


Finally, lets save the integrated data for further analysis.


```r
saveRDS(sce, "data/3pbmc_qc_dr_int_cl.rds")
```

## Check QC-stats
By now you should know how to plot different features onto your data. Take the QC metrics that were calculated in the first exercise, that should be stored in your data object, and plot it onto your UMAP and as violin plots per cluster using the clustering method of your choice. For example, plot number of UMIS, detected genes, percent mitochondrial reads. 
Then, check carefully if there is any bias in how your data is separated due to quality metrics. Could it be explained biologically, or could you have technical bias there?


### Session Info
***


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS  10.15.2
## 
## Matrix products: default
## BLAS/LAPACK: /Users/asbj/miniconda3/envs/sc_course/lib/R/lib/libRblas.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] igraph_1.2.4.1              pheatmap_1.0.12            
##  [3] rafalib_1.0.0               cowplot_1.0.0              
##  [5] scran_1.10.1                scater_1.10.1              
##  [7] ggplot2_3.2.1               SingleCellExperiment_1.4.0 
##  [9] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
## [11] BiocParallel_1.16.6         matrixStats_0.55.0         
## [13] Biobase_2.42.0              GenomicRanges_1.34.0       
## [15] GenomeInfoDb_1.18.1         IRanges_2.16.0             
## [17] S4Vectors_0.20.1            BiocGenerics_0.28.0        
## [19] RJSONIO_1.3-1.2             optparse_1.6.4             
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.1            dynamicTreeCut_1.63-1    edgeR_3.24.3            
##  [4] viridisLite_0.3.0        DelayedMatrixStats_1.4.0 assertthat_0.2.1        
##  [7] statmod_1.4.32           GenomeInfoDbData_1.2.0   vipor_0.4.5             
## [10] yaml_2.2.0               pillar_1.4.2             lattice_0.20-38         
## [13] glue_1.3.1               limma_3.38.3             digest_0.6.23           
## [16] RColorBrewer_1.1-2       XVector_0.22.0           colorspace_1.4-1        
## [19] htmltools_0.4.0          Matrix_1.2-17            plyr_1.8.4              
## [22] pkgconfig_2.0.3          zlibbioc_1.28.0          purrr_0.3.3             
## [25] scales_1.0.0             HDF5Array_1.10.1         getopt_1.20.3           
## [28] tibble_2.1.3             withr_2.1.2              lazyeval_0.2.2          
## [31] magrittr_1.5             crayon_1.3.4             evaluate_0.14           
## [34] beeswarm_0.2.3           tools_3.5.1              formatR_1.7             
## [37] stringr_1.4.0            Rhdf5lib_1.4.3           munsell_0.5.0           
## [40] locfit_1.5-9.1           compiler_3.5.1           rlang_0.4.2             
## [43] rhdf5_2.26.2             grid_3.5.1               RCurl_1.95-4.12         
## [46] BiocNeighbors_1.0.0      labeling_0.3             bitops_1.0-6            
## [49] rmarkdown_1.17           gtable_0.3.0             reshape2_1.4.3          
## [52] R6_2.4.1                 gridExtra_2.3            knitr_1.26              
## [55] dplyr_0.8.3              stringi_1.4.3            ggbeeswarm_0.6.0        
## [58] Rcpp_1.0.3               tidyselect_0.2.5         xfun_0.11
```
