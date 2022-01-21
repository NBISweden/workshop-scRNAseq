---
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 14, 2022'
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
if (!require(clustree)) {
    install.packages("clustree", dependencies = FALSE)
}
```

```
## Loading required package: clustree
```

```
## Loading required package: ggraph
```

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
# install.packages('clustree')
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
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS/LAPACK: /Users/asbj/miniconda3/envs/scRNAseq2022_tmp/lib/libopenblasp-r0.3.18.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] igraph_1.2.11               pheatmap_1.0.12            
##  [3] clustree_0.4.4              ggraph_2.0.5               
##  [5] reticulate_1.22             harmony_1.0                
##  [7] Rcpp_1.0.8                  umap_0.2.7.0               
##  [9] rafalib_1.0.0               scDblFinder_1.8.0          
## [11] DoubletFinder_2.0.3         org.Hs.eg.db_3.14.0        
## [13] AnnotationDbi_1.56.1        cowplot_1.1.1              
## [15] scran_1.22.0                scater_1.22.0              
## [17] ggplot2_3.3.5               scuttle_1.4.0              
## [19] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
## [21] Biobase_2.54.0              GenomicRanges_1.46.0       
## [23] GenomeInfoDb_1.30.0         IRanges_2.28.0             
## [25] S4Vectors_0.32.0            BiocGenerics_0.40.0        
## [27] MatrixGenerics_1.6.0        matrixStats_0.61.0         
## [29] RJSONIO_1.3-1.6             optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                tidyselect_1.1.1         
##   [3] RSQLite_2.2.8             htmlwidgets_1.5.4        
##   [5] grid_4.1.2                BiocParallel_1.28.0      
##   [7] Rtsne_0.15                munsell_0.5.0            
##   [9] ScaledMatrix_1.2.0        codetools_0.2-18         
##  [11] ica_1.0-2                 xgboost_1.5.0.1          
##  [13] statmod_1.4.36            future_1.23.0            
##  [15] miniUI_0.1.1.1            batchelor_1.10.0         
##  [17] withr_2.4.3               colorspace_2.0-2         
##  [19] highr_0.9                 knitr_1.37               
##  [21] Seurat_4.0.6              ROCR_1.0-11              
##  [23] tensor_1.5                listenv_0.8.0            
##  [25] labeling_0.4.2            GenomeInfoDbData_1.2.7   
##  [27] polyclip_1.10-0           bit64_4.0.5              
##  [29] farver_2.1.0              rprojroot_2.0.2          
##  [31] parallelly_1.30.0         vctrs_0.3.8              
##  [33] generics_0.1.1            xfun_0.29                
##  [35] R6_2.5.1                  graphlayouts_0.8.0       
##  [37] ggbeeswarm_0.6.0          rsvd_1.0.5               
##  [39] locfit_1.5-9.4            hdf5r_1.3.5              
##  [41] bitops_1.0-7              spatstat.utils_2.3-0     
##  [43] cachem_1.0.6              DelayedArray_0.20.0      
##  [45] assertthat_0.2.1          promises_1.2.0.1         
##  [47] scales_1.1.1              beeswarm_0.4.0           
##  [49] gtable_0.3.0              beachmat_2.10.0          
##  [51] globals_0.14.0            goftest_1.2-3            
##  [53] tidygraph_1.2.0           rlang_0.4.12             
##  [55] splines_4.1.2             lazyeval_0.2.2           
##  [57] checkmate_2.0.0           spatstat.geom_2.3-1      
##  [59] BiocManager_1.30.16       yaml_2.2.1               
##  [61] reshape2_1.4.4            abind_1.4-5              
##  [63] backports_1.4.1           httpuv_1.6.5             
##  [65] tools_4.1.2               ellipsis_0.3.2           
##  [67] spatstat.core_2.3-2       jquerylib_0.1.4          
##  [69] RColorBrewer_1.1-2        ggridges_0.5.3           
##  [71] plyr_1.8.6                sparseMatrixStats_1.6.0  
##  [73] zlibbioc_1.40.0           purrr_0.3.4              
##  [75] RCurl_1.98-1.5            openssl_1.4.6            
##  [77] rpart_4.1-15              deldir_1.0-6             
##  [79] pbapply_1.5-0             viridis_0.6.2            
##  [81] zoo_1.8-9                 SeuratObject_4.0.4       
##  [83] ggrepel_0.9.1             cluster_2.1.2            
##  [85] here_1.0.1                magrittr_2.0.1           
##  [87] RSpectra_0.16-0           data.table_1.14.2        
##  [89] scattermore_0.7           ResidualMatrix_1.4.0     
##  [91] lmtest_0.9-39             RANN_2.6.1               
##  [93] fitdistrplus_1.1-6        patchwork_1.1.1          
##  [95] mime_0.12                 evaluate_0.14            
##  [97] xtable_1.8-4              gridExtra_2.3            
##  [99] compiler_4.1.2            tibble_3.1.6             
## [101] KernSmooth_2.23-20        crayon_1.4.2             
## [103] htmltools_0.5.2           mgcv_1.8-38              
## [105] later_1.2.0               tidyr_1.1.4              
## [107] DBI_1.1.2                 tweenr_1.0.2             
## [109] formatR_1.11              MASS_7.3-55              
## [111] Matrix_1.4-0              getopt_1.20.3            
## [113] parallel_4.1.2            metapod_1.2.0            
## [115] pkgconfig_2.0.3           plotly_4.10.0            
## [117] spatstat.sparse_2.1-0     vipor_0.4.5              
## [119] bslib_0.3.1               dqrng_0.3.0              
## [121] XVector_0.34.0            stringr_1.4.0            
## [123] digest_0.6.29             sctransform_0.3.3        
## [125] RcppAnnoy_0.0.19          spatstat.data_2.1-2      
## [127] Biostrings_2.62.0         rmarkdown_2.11           
## [129] leiden_0.3.9              uwot_0.1.11              
## [131] edgeR_3.36.0              DelayedMatrixStats_1.16.0
## [133] curl_4.3.2                shiny_1.7.1              
## [135] lifecycle_1.0.1           nlme_3.1-155             
## [137] jsonlite_1.7.2            BiocNeighbors_1.12.0     
## [139] askpass_1.1               viridisLite_0.4.0        
## [141] limma_3.50.0              fansi_1.0.0              
## [143] pillar_1.6.4              lattice_0.20-45          
## [145] KEGGREST_1.34.0           fastmap_1.1.0            
## [147] httr_1.4.2                survival_3.2-13          
## [149] glue_1.6.0                remotes_2.4.2            
## [151] png_0.1-7                 bluster_1.4.0            
## [153] bit_4.0.4                 ggforce_0.3.3            
## [155] stringi_1.7.6             sass_0.4.0               
## [157] blob_1.2.2                BiocSingular_1.10.0      
## [159] memoise_2.0.1             dplyr_1.0.7              
## [161] irlba_2.3.5               future.apply_1.8.1
```
