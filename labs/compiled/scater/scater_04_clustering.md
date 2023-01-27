---
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 27, 2023'
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
editor_options:
  chunk_output_type: console
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

sce <- readRDS("data/results/covid_qc_dr_int.rds")
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
g <- buildKNNGraph(sce, k = 30, use.dimred = "MNN")
reducedDim(sce, "KNN") <- igraph::as_adjacency_matrix(g)

# These 2 lines are the most recommended
g <- buildSNNGraph(sce, k = 30, use.dimred = "MNN")
reducedDim(sce, "SNN") <- as_adjacency_matrix(g, attr = "weight")
```

We can take a look at the kNN graph. It is a matrix where every connection between cells is represented as $1$s. This is called a **unweighted** graph (default in Seurat). Some cell connections can however have more importance than others, in that case the scale of the graph from $0$ to a maximum distance. Usually, the smaller the distance, the closer two points are, and stronger is their connection. This is called a **weighted** graph. Both weighted and unweighted graphs are suitable for clustering, but clustering on unweighted graphs is faster for large datasets (> 100k cells).


```r
# plot the KNN graph
pheatmap(reducedDim(sce, "KNN")[1:200, 1:200], col = c("white", "black"), border_color = "grey90",
    legend = F, cluster_rows = F, cluster_cols = F, fontsize = 2)
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
# or the SNN graph
pheatmap(reducedDim(sce, "SNN")[1:200, 1:200], col = colorRampPalette(c("white",
    "yellow", "red", "black"))(20), border_color = "grey90", legend = T, cluster_rows = F,
    cluster_cols = F, fontsize = 2)
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-3-2.png)<!-- -->

As you can see, the way Scran computes the SNN graph is different to Seurat. It gives edges to all cells that shares a neighbor, but weights the edges by how similar the neighbors are. Hence, the SNN graph has more edges than the KNN graph.


### Clustering on a graph


Once the graph is built, we can now perform graph clustering. The clustering is done respective to a resolution which can be interpreted as how coarse you want your cluster to be. Higher resolution means higher number of clusters.



```r
g <- buildSNNGraph(sce, k = 5, use.dimred = "MNN")
sce$louvain_SNNk5 <- factor(cluster_louvain(g)$membership)

g <- buildSNNGraph(sce, k = 10, use.dimred = "MNN")
sce$louvain_SNNk10 <- factor(cluster_louvain(g)$membership)

g <- buildSNNGraph(sce, k = 15, use.dimred = "MNN")
sce$louvain_SNNk15 <- factor(cluster_louvain(g)$membership)

plot_grid(ncol = 3, plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = "louvain_SNNk5") +
    ggplot2::ggtitle(label = "louvain_SNNk5"), plotReducedDim(sce, dimred = "UMAP_on_MNN",
    colour_by = "louvain_SNNk10") + ggplot2::ggtitle(label = "louvain_SNNk10"), plotReducedDim(sce,
    dimred = "UMAP_on_MNN", colour_by = "louvain_SNNk15") + ggplot2::ggtitle(label = "louvain_SNNk15"))
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


We can now use the `clustree` package to visualize how cells are distributed between clusters depending on resolution.



```r
suppressPackageStartupMessages(library(clustree))

clustree(sce, prefix = "louvain_SNNk")
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## K-means clustering
***

K-means is a generic clustering algorithm that has been used in many application areas. In R, it can be applied via the kmeans function. Typically, it is applied to a reduced dimension representation of the expression data (most often PCA, because of the interpretability of the low-dimensional distances). We need to define the number of clusters in advance. Since the results depend on the initialization of the cluster centers, it is typically recommended to run K-means with multiple starting configurations (via the nstart argument).


```r
sce$kmeans_5 <- factor(kmeans(x = reducedDim(sce, "MNN"), centers = 5)$cluster)
sce$kmeans_10 <- factor(kmeans(x = reducedDim(sce, "MNN"), centers = 10)$cluster)
sce$kmeans_15 <- factor(kmeans(x = reducedDim(sce, "MNN"), centers = 15)$cluster)

plot_grid(ncol = 3, plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = "kmeans_5") +
    ggplot2::ggtitle(label = "KMeans5"), plotReducedDim(sce, dimred = "UMAP_on_MNN",
    colour_by = "kmeans_10") + ggplot2::ggtitle(label = "KMeans10"), plotReducedDim(sce,
    dimred = "UMAP_on_MNN", colour_by = "kmeans_15") + ggplot2::ggtitle(label = "KMeans15"))
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


```r
clustree(sce, prefix = "kmeans_")
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

## Hierarchical clustering
***

### Defining distance between cells

The base R `stats` package already contains a function `dist` that calculates distances between all pairs of samples. Since we want to compute distances between samples, rather than among genes, we need to transpose the data before applying it to the `dist` function. This can be done by simply adding the transpose function `t()` to the data. The distance methods available  in `dist` are: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".


```r
d <- dist(reducedDim(sce, "MNN"), method = "euclidean")
```

As you might have realized, correlation is not a method implemented in the `dist` function. However, we can create our own distances and transform them to a distance object. We can first compute sample correlations using the `cor` function.
As you already know, correlation range from -1 to 1, where 1 indicates that two samples are closest, -1 indicates that two samples are the furthest and 0 is somewhat in between. This, however, creates a problem in defining distances because a distance of 0 indicates that two samples are closest, 1 indicates that two samples are the furthest and distance of -1 is not meaningful. We thus need to transform the correlations to a positive scale (a.k.a. **adjacency**):

\[adj = \frac{1- cor}{2}\]

Once we transformed the correlations to a 0-1 scale, we can simply convert it to a distance object using `as.dist` function. The transformation does not need to have a maximum of 1, but it is more intuitive to have it at 1, rather than at any other number.


```r
# Compute sample correlations
sample_cor <- cor(Matrix::t(reducedDim(sce, "MNN")))

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
  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "hc_euclidean_5")+
    ggplot2::ggtitle(label ="HC_euclidean_5"),
  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "hc_euclidean_10")+
    ggplot2::ggtitle(label ="HC_euclidean_10"),
  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "hc_euclidean_15")+
    ggplot2::ggtitle(label ="HC_euclidean_15"),

  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "hc_corelation_5")+
    ggplot2::ggtitle(label ="HC_correlation_5"),
  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "hc_corelation_10")+
    ggplot2::ggtitle(label ="HC_correlation_10"),
  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "hc_corelation_15")+
    ggplot2::ggtitle(label ="HC_correlation_15")
)
```

![](scater_04_clustering_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


Finally, lets save the integrated data for further analysis.


```r
saveRDS(sce, "data/results/covid_qc_dr_int_cl.rds")
```

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

By now you should know how to plot different features onto your data. Take the QC metrics that were calculated in the first exercise, that should be stored in your data object, and plot it as violin plots per cluster using the clustering method of your choice. For example, plot number of UMIS, detected genes, percent mitochondrial reads.

Then, check carefully if there is any bias in how your data is separated due to quality metrics. Could it be explained biologically, or could you have technical bias there?
</div>


### Session Info
***


```r
sessionInfo()
```

```
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Big Sur/Monterey 10.16
## 
## Matrix products: default
## BLAS/LAPACK: /Users/asabjor/miniconda3/envs/scRNAseq2023/lib/libopenblasp-r0.3.21.dylib
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] clustree_0.5.0              ggraph_2.1.0               
##  [3] igraph_1.3.5                pheatmap_1.0.12            
##  [5] rafalib_1.0.0               cowplot_1.1.1              
##  [7] scran_1.22.1                scater_1.22.0              
##  [9] ggplot2_3.4.0               scuttle_1.4.0              
## [11] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
## [13] Biobase_2.54.0              GenomicRanges_1.46.1       
## [15] GenomeInfoDb_1.30.1         IRanges_2.28.0             
## [17] S4Vectors_0.32.4            BiocGenerics_0.40.0        
## [19] MatrixGenerics_1.6.0        matrixStats_0.63.0         
## [21] RJSONIO_1.3-1.7             optparse_1.7.3             
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7              RColorBrewer_1.1-3       
##  [3] backports_1.4.1           tools_4.1.3              
##  [5] bslib_0.4.2               utf8_1.2.2               
##  [7] R6_2.5.1                  irlba_2.3.5.1            
##  [9] vipor_0.4.5               DBI_1.1.3                
## [11] colorspace_2.1-0          withr_2.5.0              
## [13] tidyselect_1.2.0          gridExtra_2.3            
## [15] compiler_4.1.3            cli_3.6.0                
## [17] BiocNeighbors_1.12.0      formatR_1.14             
## [19] DelayedArray_0.20.0       labeling_0.4.2           
## [21] sass_0.4.5                checkmate_2.1.0          
## [23] scales_1.2.1              stringr_1.5.0            
## [25] digest_0.6.31             rmarkdown_2.20           
## [27] XVector_0.34.0            pkgconfig_2.0.3          
## [29] htmltools_0.5.4           sparseMatrixStats_1.6.0  
## [31] highr_0.10                fastmap_1.1.0            
## [33] limma_3.50.3              rlang_1.0.6              
## [35] DelayedMatrixStats_1.16.0 farver_2.1.1             
## [37] jquerylib_0.1.4           generics_0.1.3           
## [39] jsonlite_1.8.4            BiocParallel_1.28.3      
## [41] dplyr_1.0.10              RCurl_1.98-1.9           
## [43] magrittr_2.0.3            BiocSingular_1.10.0      
## [45] GenomeInfoDbData_1.2.7    Matrix_1.5-3             
## [47] Rcpp_1.0.10               ggbeeswarm_0.7.1         
## [49] munsell_0.5.0             fansi_1.0.4              
## [51] viridis_0.6.2             lifecycle_1.0.3          
## [53] stringi_1.7.12            yaml_2.3.7               
## [55] edgeR_3.36.0              MASS_7.3-58.2            
## [57] zlibbioc_1.40.0           grid_4.1.3               
## [59] parallel_4.1.3            ggrepel_0.9.2            
## [61] dqrng_0.3.0               lattice_0.20-45          
## [63] graphlayouts_0.8.4        beachmat_2.10.0          
## [65] locfit_1.5-9.7            metapod_1.2.0            
## [67] knitr_1.41                pillar_1.8.1             
## [69] ScaledMatrix_1.2.0        glue_1.6.2               
## [71] evaluate_0.20             tweenr_2.0.2             
## [73] vctrs_0.5.2               polyclip_1.10-4          
## [75] purrr_1.0.1               tidyr_1.2.1              
## [77] gtable_0.3.1              getopt_1.20.3            
## [79] assertthat_0.2.1          cachem_1.0.6             
## [81] ggforce_0.4.1             xfun_0.36                
## [83] rsvd_1.0.5                tidygraph_1.2.2          
## [85] viridisLite_0.4.1         tibble_3.1.8             
## [87] beeswarm_0.4.0            cluster_2.1.4            
## [89] bluster_1.4.0             statmod_1.5.0
```
