---
title: "Scater/Scran: Dimensionality reduction"
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

# Dimensionality reduction

Paulo Czarnewski


<style>
h1, .h1, h2, .h2, h3, .h3, h4, .h4 { margin-top: 50px }
p.caption {font-size: 0.9em;font-style: italic;color: grey;margin-right: 10%;margin-left: 10%;text-align: justify}
</style>

## Data preparation
***

First, let's load all necessary libraries and the QC-filtered dataset from the previous step.


```r
suppressPackageStartupMessages({
    library(scater)
    library(scran)
    library(cowplot)
    library(ggplot2)
    library(rafalib)
    library(umap)
})

sce <- readRDS("data/results/covid_qc.rds")
```

### Feature selection

Next, we first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.


```r
sce <- computeSumFactors(sce, sizes = c(20, 40, 60, 80))
sce <- logNormCounts(sce)
var.out <- modelGeneVar(sce, method = "loess")
hvgs = getTopHVGs(var.out, n = 2000)

mypar(1, 2)
# plot mean over TOTAL variance Visualizing the fit:
fit.var <- metadata(var.out)
plot(fit.var$mean, fit.var$var, xlab = "Mean of log-expression", ylab = "Variance of log-expression")
curve(fit.var$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

# Select 1000 top variable genes
hvg.out <- getTopHVGs(var.out, n = 1000)

# highligt those cells in the plot
cutoff <- rownames(var.out) %in% hvg.out
points(fit.var$mean[cutoff], fit.var$var[cutoff], col = "red", pch = 16, cex = 0.6)


# plot mean over BIOLOGICAL variance
plot(var.out$mean, var.out$bio, pch = 16, cex = 0.4, xlab = "Mean log-expression",
    ylab = "Variance of log-expression")
lines(c(min(var.out$mean), max(var.out$mean)), c(0, 0), col = "dodgerblue", lwd = 2)
points(var.out$mean[cutoff], var.out$bio[cutoff], col = "red", pch = 16, cex = 0.6)
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

### Z-score transformation

Now that the data is prepared, we now proceed with PCA. Since each gene has a different expression level, it means that genes with higher expression values will naturally have higher variation that will be captured by PCA. This means that we need to somehow give each gene a similar weight when performing PCA (see below). The common practice is to center and scale each gene before performing PCA. This exact scaling is called Z-score normalization it is very useful for PCA, clustering and plotting heatmaps. <br>Additionally, we can use regression to remove any unwanted sources of variation from the dataset, such as `cell cycle`, `sequencing depth`, `percent mitocondria`. This is achieved by doing a generalized linear regression using these parameters as covariates in the model. Then the residuals of the model are taken as the "regressed data". Although perhaps not in the best way, batch effect regression can also be done here.

By default variables are scaled in the PCA step and is not done separately. But it could be acheieved by running the commads below:


```r
# sce@assays$data@listData$scaled.data <-
# apply(exprs(sce)[rownames(hvg.out),,drop=FALSE],2,function(x) scale(x,T,T))
# rownames(sce@assays$data@listData$scaled.data) <- rownames(hvg.out)
```


## PCA
***

Performing PCA has many useful applications and interpretations, which much depends on the data used. In the case of life sciences, we want to segregate samples based on gene expression patterns in the data.

As said above, we use the `logcounts` and then set `scale_features` to TRUE in order to scale each gene.


```r
# runPCA and specify the variable genes to use for dim reduction with
# subset_row
sce <- runPCA(sce, exprs_values = "logcounts", ncomponents = 50, subset_row = hvg.out,
    scale = TRUE)
```

We then plot the first principal components.


```r
plot_grid(ncol = 3, plotReducedDim(sce, dimred = "PCA", colour_by = "sample", ncomponents = 1:2,
    point_size = 0.6), plotReducedDim(sce, dimred = "PCA", colour_by = "sample",
    ncomponents = 3:4, point_size = 0.6), plotReducedDim(sce, dimred = "PCA", colour_by = "sample",
    ncomponents = 5:6, point_size = 0.6))
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

To identify which genes (Seurat) or metadata paramters (Scater/Scran) contribute the most to each PC, one can retreive the loading matrix information. Unfortunatelly this is not implemented in Scater/Scran, so you will need to compute PCA using `logcounts`.


```r
plot_grid(ncol = 2, plotExplanatoryPCs(sce))
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

We can also plot the amount of variance explained by each PC.


```r
mypar()
plot(attr(reducedDim(sce, "PCA"), "percentVar")[1:50] * 100, type = "l", ylab = "% variance",
    xlab = "Principal component #")
points(attr(reducedDim(sce, "PCA"), "percentVar")[1:50] * 100, pch = 21, bg = "grey",
    cex = 0.5)
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Based on this plot, we can see that the top 8 PCs retain a lot of information, while other PCs contain pregressivelly less. However, it is still advisable to use more PCs since they might contain informaktion about rare cell types (such as platelets and DCs in this dataset)

## tSNE
***

We can now run [BH-tSNE](https://arxiv.org/abs/1301.3342).


```r
set.seed(42)
sce <- runTSNE(sce, dimred = "PCA", n_dimred = 30, perplexity = 30, name = "tSNE_on_PCA")
```

We can now plot the tSNE colored per dataset. We can clearly see the effect of batches present in the dataset.


```r
plot_grid(ncol = 3, plotReducedDim(sce, dimred = "tSNE_on_PCA", colour_by = "sample"))
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


## UMAP
***

We can now run [UMAP](https://arxiv.org/abs/1802.03426) for cell embeddings.


```r
sce <- runUMAP(sce, dimred = "PCA", n_dimred = 30, ncomponents = 2, name = "UMAP_on_PCA")
# see ?umap and ?runUMAP for more info
```

Another usefullness of UMAP is that it is not limitted by the number of dimensions the data cen be reduced into (unlike tSNE). We can simply reduce the dimentions altering the `n.components` parameter.


```r
sce <- runUMAP(sce, dimred = "PCA", n_dimred = 30, ncomponents = 10, name = "UMAP10_on_PCA")
# see ?umap and ?runUMAP for more info
```

We can now plot the UMAP colored per dataset. Although less distinct as in the tSNE, we still see quite an effect of the different batches in the data.


```r
plot_grid(ncol = 3, plotReducedDim(sce, dimred = "UMAP_on_PCA", colour_by = "sample") +
    ggplot2::ggtitle(label = "UMAP_on_PCA"), plotReducedDim(sce, dimred = "UMAP10_on_PCA",
    colour_by = "sample", ncomponents = 1:2) + ggplot2::ggtitle(label = "UMAP10_on_PCA"),
    plotReducedDim(sce, dimred = "UMAP10_on_PCA", colour_by = "sample", ncomponents = 3:4) +
        ggplot2::ggtitle(label = "UMAP10_on_PCA"))
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-12-1.png)<!-- -->


## Using ScaledData and graphs for DR
***

Althought running a sencond dimmensionality reduction (i.e tSNE or UMAP) on PCA would be a standard approach (because it allows higher computation efficiency), the options are actually limiteless. Below we will show a couple of other common options such as running directly on the scaled data (which was used for PCA) or on a graph built from scaled data. We will show from now on only UMAP, but the same applies for tSNE.

### Using ScaledData for UMAP

To run tSNE or UMAP on the scaled data, one firts needs to select the number of variables to use. This is because including dimentions that do contribute to the separation of your cell types will in the end mask those differences. Another reason for it is because running with all genes/features also will take longer or might be computationally unfeasible. Therefore we will use the scaled data of the highly variable genes.


```r
sce <- runUMAP(sce, exprs_values = "logcounts", name = "UMAP_on_ScaleData")
```

To run tSNE or UMAP on the a graph, we first need to build a graph from the data. In fact, both tSNE and UMAP first build a graph from the data using a specified distance metrix and then optimize the embedding. Since a graph is just a matrix containing distances from cell to cell and as such, you can run either UMAP or tSNE using any other distance metric desired. Euclidean and Correlation are ususally the most commonly used.

### Using a Graph for UMAP


```r
# Build Graph
nn <- RANN::nn2(reducedDim(sce, "PCA"), k = 30)
names(nn) <- c("idx", "dist")
g <- buildKNNGraph(sce, k = 30, use.dimred = "PCA")
reducedDim(sce, "KNN") <- igraph::as_adjacency_matrix(g)


# Run UMAP and rename it for comparisson temp <- umap::umap.defaults
try(reducedDim(sce, "UMAP_on_Graph") <- NULL)
reducedDim(sce, "UMAP_on_Graph") <- uwot::umap(X = NULL, n_components = 2, nn_method = nn)
```


We can now plot the UMAP comparing both on PCA vs ScaledSata vs Graph.


```r
plot_grid(ncol = 3, plotReducedDim(sce, dimred = "UMAP_on_PCA", colour_by = "sample") +
    ggplot2::ggtitle(label = "UMAP_on_PCA"), plotReducedDim(sce, dimred = "UMAP_on_ScaleData",
    colour_by = "sample") + ggplot2::ggtitle(label = "UMAP_on_ScaleData"), plotReducedDim(sce,
    dimred = "UMAP_on_Graph", colour_by = "sample") + ggplot2::ggtitle(label = "UMAP_on_Graph"))
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

## Ploting genes of interest
***


Let's plot some marker genes for different celltypes onto the embedding. Some genes are:

Markers	| Cell Type
--- | ---
CD3E	| T cells
CD3E CD4	| CD4+ T cells
CD3E CD8A	| CD8+ T cells
GNLY, NKG7	| NK cells
MS4A1	| B cells
CD14, LYZ, CST3, MS4A7	| CD14+ Monocytes
FCGR3A, LYZ, CST3, MS4A7	| FCGR3A+  Monocytes
FCER1A, CST3 | DCs


```r
plotlist <- list()
for (i in c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7",
    "FCGR3A", "CST3", "FCER1A")) {
    plotlist[[i]] <- plotReducedDim(sce, dimred = "UMAP_on_PCA", colour_by = i, by_exprs_values = "logcounts") +
        scale_fill_gradientn(colours = colorRampPalette(c("grey90", "orange3", "firebrick",
            "firebrick", "red", "red"))(10)) + ggtitle(label = i) + theme(plot.title = element_text(size = 20))
}
plot_grid(ncol = 3, plotlist = plotlist)
```

![](scater_02_dim_reduction_files/figure-html/unnamed-chunk-16-1.png)<!-- -->


We can finally save the object for use in future steps.

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

Select some of your dimensionality reductions and plot some of the QC stats that were calculated in the previous lab. Can you see if some of the separation in your data is driven by quality of the cells?
</div>


```r
saveRDS(sce, "data/results/covid_qc_dm.rds")
```

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
##  [1] umap_0.2.7.0                rafalib_1.0.0              
##  [3] scDblFinder_1.8.0           DoubletFinder_2.0.3        
##  [5] org.Hs.eg.db_3.14.0         AnnotationDbi_1.56.1       
##  [7] cowplot_1.1.1               scran_1.22.0               
##  [9] scater_1.22.0               ggplot2_3.3.5              
## [11] scuttle_1.4.0               SingleCellExperiment_1.16.0
## [13] SummarizedExperiment_1.24.0 Biobase_2.54.0             
## [15] GenomicRanges_1.46.0        GenomeInfoDb_1.30.0        
## [17] IRanges_2.28.0              S4Vectors_0.32.0           
## [19] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
## [21] matrixStats_0.61.0          RJSONIO_1.3-1.6            
## [23] optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                reticulate_1.22          
##   [3] tidyselect_1.1.1          RSQLite_2.2.8            
##   [5] htmlwidgets_1.5.4         grid_4.1.2               
##   [7] BiocParallel_1.28.0       Rtsne_0.15               
##   [9] munsell_0.5.0             ScaledMatrix_1.2.0       
##  [11] codetools_0.2-18          ica_1.0-2                
##  [13] xgboost_1.5.0.1           statmod_1.4.36           
##  [15] future_1.23.0             miniUI_0.1.1.1           
##  [17] withr_2.4.3               colorspace_2.0-2         
##  [19] highr_0.9                 knitr_1.37               
##  [21] Seurat_4.0.6              ROCR_1.0-11              
##  [23] tensor_1.5                listenv_0.8.0            
##  [25] labeling_0.4.2            GenomeInfoDbData_1.2.7   
##  [27] polyclip_1.10-0           bit64_4.0.5              
##  [29] farver_2.1.0              parallelly_1.30.0        
##  [31] vctrs_0.3.8               generics_0.1.1           
##  [33] xfun_0.29                 R6_2.5.1                 
##  [35] ggbeeswarm_0.6.0          rsvd_1.0.5               
##  [37] locfit_1.5-9.4            hdf5r_1.3.5              
##  [39] bitops_1.0-7              spatstat.utils_2.3-0     
##  [41] cachem_1.0.6              DelayedArray_0.20.0      
##  [43] assertthat_0.2.1          promises_1.2.0.1         
##  [45] scales_1.1.1              beeswarm_0.4.0           
##  [47] gtable_0.3.0              beachmat_2.10.0          
##  [49] globals_0.14.0            goftest_1.2-3            
##  [51] rlang_0.4.12              splines_4.1.2            
##  [53] lazyeval_0.2.2            spatstat.geom_2.3-1      
##  [55] BiocManager_1.30.16       yaml_2.2.1               
##  [57] reshape2_1.4.4            abind_1.4-5              
##  [59] httpuv_1.6.5              tools_4.1.2              
##  [61] ellipsis_0.3.2            spatstat.core_2.3-2      
##  [63] jquerylib_0.1.4           RColorBrewer_1.1-2       
##  [65] ggridges_0.5.3            Rcpp_1.0.8               
##  [67] plyr_1.8.6                sparseMatrixStats_1.6.0  
##  [69] zlibbioc_1.40.0           purrr_0.3.4              
##  [71] RCurl_1.98-1.5            openssl_1.4.6            
##  [73] rpart_4.1-15              deldir_1.0-6             
##  [75] pbapply_1.5-0             viridis_0.6.2            
##  [77] zoo_1.8-9                 SeuratObject_4.0.4       
##  [79] ggrepel_0.9.1             cluster_2.1.2            
##  [81] magrittr_2.0.1            RSpectra_0.16-0          
##  [83] data.table_1.14.2         scattermore_0.7          
##  [85] lmtest_0.9-39             RANN_2.6.1               
##  [87] fitdistrplus_1.1-6        patchwork_1.1.1          
##  [89] mime_0.12                 evaluate_0.14            
##  [91] xtable_1.8-4              gridExtra_2.3            
##  [93] compiler_4.1.2            tibble_3.1.6             
##  [95] KernSmooth_2.23-20        crayon_1.4.2             
##  [97] htmltools_0.5.2           mgcv_1.8-38              
##  [99] later_1.2.0               tidyr_1.1.4              
## [101] DBI_1.1.2                 formatR_1.11             
## [103] MASS_7.3-55               Matrix_1.4-0             
## [105] getopt_1.20.3             parallel_4.1.2           
## [107] metapod_1.2.0             igraph_1.2.11            
## [109] pkgconfig_2.0.3           plotly_4.10.0            
## [111] spatstat.sparse_2.1-0     vipor_0.4.5              
## [113] bslib_0.3.1               dqrng_0.3.0              
## [115] XVector_0.34.0            stringr_1.4.0            
## [117] digest_0.6.29             sctransform_0.3.3        
## [119] RcppAnnoy_0.0.19          spatstat.data_2.1-2      
## [121] Biostrings_2.62.0         rmarkdown_2.11           
## [123] leiden_0.3.9              uwot_0.1.11              
## [125] edgeR_3.36.0              DelayedMatrixStats_1.16.0
## [127] curl_4.3.2                shiny_1.7.1              
## [129] lifecycle_1.0.1           nlme_3.1-155             
## [131] jsonlite_1.7.2            BiocNeighbors_1.12.0     
## [133] askpass_1.1               viridisLite_0.4.0        
## [135] limma_3.50.0              fansi_1.0.0              
## [137] pillar_1.6.4              lattice_0.20-45          
## [139] KEGGREST_1.34.0           fastmap_1.1.0            
## [141] httr_1.4.2                survival_3.2-13          
## [143] glue_1.6.0                remotes_2.4.2            
## [145] png_0.1-7                 bluster_1.4.0            
## [147] bit_4.0.4                 stringi_1.7.6            
## [149] sass_0.4.0                blob_1.2.2               
## [151] BiocSingular_1.10.0       memoise_2.0.1            
## [153] dplyr_1.0.7               irlba_2.3.5              
## [155] future.apply_1.8.1
```
