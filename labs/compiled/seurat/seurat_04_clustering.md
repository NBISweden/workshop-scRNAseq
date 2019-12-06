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
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(rafalib)
})

alldata <- readRDS("data/3pbmc_qc_dr_int.rds")
```

## Graph clustering
***

The procedure of clustering on a Graph can be generalized as 3 main steps:

1) Build a kNN graph from the data

2) Prune spurious connections from kNN graph (optional step). This is a SNN graph.

3) Find groups of cells that maximizes the connections within the group compared other groups.

### Building kNN / SNN graph


The first step into graph clustering is to construct a k-nn graph, in case you don't have one. For this, we will use the PCA space. Thus, as done for dimensionality reduction, we will use ony the top *N* PCA dimensions for this purpose (the same used for computing UMAP / tSNE).

As we can see above, the **Seurat** function `FindNeighbors` already computes both the KNN and SNN graphs, in which we can control the minimal percentage of shared neighbours to be kept. See `?FindNeighbors` for additional options.


```r
alldata <- FindNeighbors(alldata,
                         reduction = "PCA_on_CCA",
                         dims = 1:30,
                         k.param = 60,
                         prune.SNN = 1/15)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
names(alldata@graphs)
```

```
## [1] "CCA_nn"  "CCA_snn"
```

We can take a look at the kNN graph. It is a matrix where every connection between cells is represented as $1$s. This is called a **unweighted** graph (default in Seurat). Some cell connections can however have more importance than others, in that case the scale of the graph from $0$ to a maximum distance. Usually, the smaller the distance, the closer two points are, and stronger is their connection. This is called a **weighted** graph. Both weighted and unweighted graphs are suitable for clustering, but clustering on unweighted graphs is faster for large datasets (> 100k cells).


```r
pheatmap(alldata@graphs$CCA_nn[1:200,1:200],
         col=c("white","black"),border_color = "grey90",
         legend = F,cluster_rows = F,cluster_cols = F,fontsize = 2)
```

![](seurat_04_clustering_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

### Clustering on a graph


Once the graph is built, we can now perform graph clustering. The clustering is done respective to a resolution which can be interpreted as how coarse you want your cluster to be. Higher resolution means higher number of clusters.

In **Seurat**, the function `FindClusters` will do a graph-based clustering using "Louvain" algorithim by default (`algorithm = 1`). TO use the leiden algorithm, you need to set it to `algorithm = 4`. See `?FindClusters` for additional options.


```r
alldata <- FindClusters(alldata, graph.name = "CCA_snn", resolution = .5 , algorithm = 1)
alldata <- FindClusters(alldata, graph.name = "CCA_snn", resolution = 1  , algorithm = 1)
alldata <- FindClusters(alldata, graph.name = "CCA_snn", resolution = 2  , algorithm = 1)

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "CCA_snn_res.0.5")+ggtitle("louvain_0.5"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "CCA_snn_res.1")+ggtitle("louvain_1"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "CCA_snn_res.2")+ggtitle("louvain_2")
)
```

![](seurat_04_clustering_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


## K-means clustering
***

K-means is a generic clustering algorithm that has been used in many application areas. In R, it can be applied via the kmeans function. Typically, it is applied to a reduced dimension representation of the expression data (most often PCA, because of the interpretability of the low-dimensional distances). We need to define the number of clusters in advance. Since the results depend on the initialization of the cluster centers, it is typically recommended to run K-means with multiple starting configurations (via the nstart argument).


```r
alldata$kmeans_5 <- kmeans(x = alldata@reductions[["PCA_on_CCA"]]@cell.embeddings,centers = 5)$cluster
alldata$kmeans_10 <- kmeans(x = alldata@reductions[["PCA_on_CCA"]]@cell.embeddings,centers = 10)$cluster
alldata$kmeans_15 <- kmeans(x = alldata@reductions[["PCA_on_CCA"]]@cell.embeddings,centers = 15)$cluster

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "kmeans_5")+ggtitle("kmeans_5"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "kmeans_10")+ggtitle("kmeans_10"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "kmeans_15")+ggtitle("kmeans_15")
)
```

![](seurat_04_clustering_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


## Hierarchical clustering
***

### Defining distance between cells

The base R `stats` package already contains a function `dist` that calculates distances between all pairs of samples. Since we want to compute distances between samples, rather than among genes, we need to transpose the data before applying it to the `dist` function. This can be done by simply adding the transpose function `t()` to the data. The distance methods available  in `dist` are: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".


```r
d <- dist( alldata@reductions[["PCA_on_CCA"]]@cell.embeddings,
           method="euclidean")
```

As you might have realized, correlation is not a method implemented in the `dist` function. However, we can create our own distances and transform them to a distance object. We can first compute sample correlations using the `cor` function.
As you already know, correlation range from -1 to 1, where 1 indicates that two samples are closest, -1 indicates that two samples are the furthest and 0 is somewhat in between. This, however, creates a problem in defining distances because a distance of 0 indicates that two samples are closest, 1 indicates that two samples are the furthest and distance of -1 is not meaningful. We thus need to transform the correlations to a positive scale (a.k.a. **adjacency**):

\[adj = \frac{1- cor}{2}\]

Once we transformed the correlations to a 0-1 scale, we can simply convert it to a distance object using `as.dist` function. The transformation does not need to have a maximum of 1, but it is more intuitive to have it at 1, rather than at any other number.


```r
#Compute sample correlations
sample_cor <- cor( Matrix::t(alldata@reductions[["PCA_on_CCA"]]@cell.embeddings) )

#Transform the scale from correlations
sample_cor <- (1 - sample_cor) / 2

#Convert it to a distance object
d2 <- as.dist(sample_cor)
```

### Defining distance between cells

After having calculated the distances between samples calculated, we can now proceed with the hierarchical clustering per-se. We will use the function `hclust` for this purpose, in which we can simply run it with the distance objects created above. The methods available are: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid". It is possible to plot the dendrogram for all cells, but this is very time consuming and we will omit for this tutorial.


```r
#euclidean
h_euclidean <- hclust(d, method="ward.D2")

#correlation
h_correlation <- hclust(d2, method="ward.D2")
```

 Once your dendrogram is created, the next step is to define which samples belong to a particular cluster. After identifying the dendrogram, we can now literally cut the tree at a fixed threshold (with `cutree`) at different levels to define the clusters. We can either define the number of clusters or decide on a height. We can simply try different clustering levels.


```r
#euclidean distance
alldata$hc_euclidean_5 <- cutree(h_euclidean,k = 5)
alldata$hc_euclidean_10 <- cutree(h_euclidean,k = 10)
alldata$hc_euclidean_15 <- cutree(h_euclidean,k = 15)

#correlation distance
alldata$hc_corelation_5 <- cutree(h_correlation,k = 5)
alldata$hc_corelation_10 <- cutree(h_correlation,k = 10)
alldata$hc_corelation_15 <- cutree(h_correlation,k = 15)


plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "hc_euclidean_5")+ggtitle("hc_euc_5"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "hc_euclidean_10")+ggtitle("hc_euc_10"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "hc_euclidean_15")+ggtitle("hc_euc_15"),
  
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "hc_corelation_5")+ggtitle("hc_cor_5"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "hc_corelation_10")+ggtitle("hc_cor_10"),
  DimPlot(alldata, reduction = "UMAP_on_CCA", group.by = "hc_corelation_15")+ggtitle("hc_cor_15")
)
```

![](seurat_04_clustering_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


Finally, lets save the integrated data for further analysis.


```r
saveRDS(alldata,"data/3pbmc_qc_dr_int_cl.rds")
```


### Additional Materials
***
<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 20px;}
</style>
<div class = "blue">
If you would like to learn more about building graphs and clustering, please follow this other tutorial
</div>


### Session Info
***


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS  10.15
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
##  [1] rafalib_1.0.0               pheatmap_1.0.12            
##  [3] scran_1.10.1                SingleCellExperiment_1.4.0 
##  [5] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
##  [7] matrixStats_0.55.0          Biobase_2.42.0             
##  [9] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1        
## [11] IRanges_2.16.0              S4Vectors_0.20.1           
## [13] BiocGenerics_0.28.0         BiocParallel_1.16.6        
## [15] ggplot2_3.2.1               cowplot_1.0.0              
## [17] Matrix_1.2-17               Seurat_3.0.1               
## [19] RJSONIO_1.3-1.2             optparse_1.6.4             
## 
## loaded via a namespace (and not attached):
##   [1] ggbeeswarm_0.6.0         Rtsne_0.15               colorspace_1.4-1        
##   [4] ggridges_0.5.1           dynamicTreeCut_1.63-1    XVector_0.22.0          
##   [7] BiocNeighbors_1.0.0      listenv_0.7.0            npsurv_0.4-0            
##  [10] getopt_1.20.3            ggrepel_0.8.1            bit64_0.9-7             
##  [13] codetools_0.2-16         splines_3.5.1            R.methodsS3_1.7.1       
##  [16] lsei_1.2-0               scater_1.10.1            knitr_1.26              
##  [19] zeallot_0.1.0            jsonlite_1.6             ica_1.0-2               
##  [22] cluster_2.1.0            png_0.1-7                R.oo_1.23.0             
##  [25] HDF5Array_1.10.1         sctransform_0.2.0        compiler_3.5.1          
##  [28] httr_1.4.1               backports_1.1.5          assertthat_0.2.1        
##  [31] lazyeval_0.2.2           limma_3.38.3             htmltools_0.4.0         
##  [34] tools_3.5.1              rsvd_1.0.2               igraph_1.2.4.1          
##  [37] gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.0  
##  [40] RANN_2.6.1               reshape2_1.4.3           dplyr_0.8.3             
##  [43] Rcpp_1.0.3               vctrs_0.2.0              gdata_2.18.0            
##  [46] ape_5.3                  nlme_3.1-141             DelayedMatrixStats_1.4.0
##  [49] gbRd_0.4-11              lmtest_0.9-37            xfun_0.11               
##  [52] stringr_1.4.0            globals_0.12.4           lifecycle_0.1.0         
##  [55] irlba_2.3.3              gtools_3.8.1             statmod_1.4.32          
##  [58] future_1.15.1            edgeR_3.24.3             MASS_7.3-51.4           
##  [61] zlibbioc_1.28.0          zoo_1.8-6                scales_1.0.0            
##  [64] rhdf5_2.26.2             RColorBrewer_1.1-2       yaml_2.2.0              
##  [67] reticulate_1.13          pbapply_1.4-2            gridExtra_2.3           
##  [70] stringi_1.4.3            venn_1.7                 caTools_1.17.1.2        
##  [73] bibtex_0.4.2             Rdpack_0.11-0            SDMTools_1.1-221.1      
##  [76] rlang_0.4.2              pkgconfig_2.0.3          bitops_1.0-6            
##  [79] evaluate_0.14            lattice_0.20-38          Rhdf5lib_1.4.3          
##  [82] ROCR_1.0-7               purrr_0.3.3              htmlwidgets_1.5.1       
##  [85] labeling_0.3             bit_1.1-14               tidyselect_0.2.5        
##  [88] plyr_1.8.4               magrittr_1.5             R6_2.4.1                
##  [91] gplots_3.0.1.1           pillar_1.4.2             withr_2.1.2             
##  [94] fitdistrplus_1.0-14      survival_2.44-1.1        RCurl_1.95-4.12         
##  [97] tibble_2.1.3             future.apply_1.3.0       tsne_0.1-3              
## [100] crayon_1.3.4             hdf5r_1.2.0              KernSmooth_2.23-15      
## [103] plotly_4.9.1             rmarkdown_1.17           viridis_0.5.1           
## [106] locfit_1.5-9.1           grid_3.5.1               data.table_1.11.6       
## [109] metap_1.1                digest_0.6.23            tidyr_1.0.0             
## [112] R.utils_2.9.0            munsell_0.5.0            beeswarm_0.2.3          
## [115] viridisLite_0.3.0        vipor_0.4.5
```
