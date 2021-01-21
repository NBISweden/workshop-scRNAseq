---
title: "Scater/Scran: Dimensionality reduction"
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 21, 2021'
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
hvgs = getTopHVGs(dec, n = 2000)

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
# runPCA and specify the variable genes to use for dim reduction with subset_row
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


```r
saveRDS(sce, "data/results/covid_qc_dm.rds")
```

### Session Info
***


```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Catalina 10.15.5
## 
## Matrix products: default
## BLAS/LAPACK: /Users/paulo.czarnewski/.conda/envs/scRNAseq2021/lib/libopenblasp-r0.3.12.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] umap_0.2.7.0                rafalib_1.0.0              
##  [3] scDblFinder_1.4.0           biomaRt_2.46.0             
##  [5] DoubletFinder_2.0.3         org.Hs.eg.db_3.12.0        
##  [7] AnnotationDbi_1.52.0        cowplot_1.1.1              
##  [9] scran_1.18.0                scater_1.18.0              
## [11] ggplot2_3.3.3               SingleCellExperiment_1.12.0
## [13] SummarizedExperiment_1.20.0 Biobase_2.50.0             
## [15] GenomicRanges_1.42.0        GenomeInfoDb_1.26.0        
## [17] IRanges_2.24.0              S4Vectors_0.28.0           
## [19] BiocGenerics_0.36.0         MatrixGenerics_1.2.0       
## [21] matrixStats_0.57.0          RJSONIO_1.3-1.4            
## [23] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] reticulate_1.18           tidyselect_1.1.0         
##   [3] RSQLite_2.2.2             htmlwidgets_1.5.3        
##   [5] grid_4.0.3                BiocParallel_1.24.0      
##   [7] Rtsne_0.15                munsell_0.5.0            
##   [9] codetools_0.2-18          ica_1.0-2                
##  [11] statmod_1.4.35            xgboost_1.3.0.1          
##  [13] future_1.21.0             miniUI_0.1.1.1           
##  [15] withr_2.4.0               colorspace_2.0-0         
##  [17] knitr_1.30                Seurat_3.2.3             
##  [19] ROCR_1.0-11               tensor_1.5               
##  [21] listenv_0.8.0             labeling_0.4.2           
##  [23] GenomeInfoDbData_1.2.4    polyclip_1.10-0          
##  [25] bit64_4.0.5               farver_2.0.3             
##  [27] parallelly_1.23.0         vctrs_0.3.6              
##  [29] generics_0.1.0            xfun_0.20                
##  [31] BiocFileCache_1.14.0      R6_2.5.0                 
##  [33] ggbeeswarm_0.6.0          rsvd_1.0.3               
##  [35] locfit_1.5-9.4            hdf5r_1.3.3              
##  [37] bitops_1.0-6              spatstat.utils_1.20-2    
##  [39] DelayedArray_0.16.0       assertthat_0.2.1         
##  [41] promises_1.1.1            scales_1.1.1             
##  [43] beeswarm_0.2.3            gtable_0.3.0             
##  [45] beachmat_2.6.0            globals_0.14.0           
##  [47] goftest_1.2-2             rlang_0.4.10             
##  [49] splines_4.0.3             lazyeval_0.2.2           
##  [51] yaml_2.2.1                reshape2_1.4.4           
##  [53] abind_1.4-5               httpuv_1.5.5             
##  [55] tools_4.0.3               ellipsis_0.3.1           
##  [57] RColorBrewer_1.1-2        ggridges_0.5.3           
##  [59] Rcpp_1.0.6                plyr_1.8.6               
##  [61] sparseMatrixStats_1.2.0   progress_1.2.2           
##  [63] zlibbioc_1.36.0           purrr_0.3.4              
##  [65] RCurl_1.98-1.2            prettyunits_1.1.1        
##  [67] rpart_4.1-15              openssl_1.4.3            
##  [69] deldir_0.2-9              pbapply_1.4-3            
##  [71] viridis_0.5.1             zoo_1.8-8                
##  [73] ggrepel_0.9.1             cluster_2.1.0            
##  [75] magrittr_2.0.1            data.table_1.13.6        
##  [77] RSpectra_0.16-0           scattermore_0.7          
##  [79] lmtest_0.9-38             RANN_2.6.1               
##  [81] fitdistrplus_1.1-3        hms_1.0.0                
##  [83] patchwork_1.1.1           mime_0.9                 
##  [85] evaluate_0.14             xtable_1.8-4             
##  [87] XML_3.99-0.5              gridExtra_2.3            
##  [89] compiler_4.0.3            tibble_3.0.5             
##  [91] KernSmooth_2.23-18        crayon_1.3.4             
##  [93] htmltools_0.5.1           mgcv_1.8-33              
##  [95] later_1.1.0.1             tidyr_1.1.2              
##  [97] DBI_1.1.1                 formatR_1.7              
##  [99] dbplyr_2.0.0              MASS_7.3-53              
## [101] rappdirs_0.3.1            Matrix_1.3-2             
## [103] getopt_1.20.3             igraph_1.2.6             
## [105] pkgconfig_2.0.3           plotly_4.9.3             
## [107] scuttle_1.0.0             xml2_1.3.2               
## [109] vipor_0.4.5               dqrng_0.2.1              
## [111] XVector_0.30.0            stringr_1.4.0            
## [113] digest_0.6.27             sctransform_0.3.2        
## [115] RcppAnnoy_0.0.18          spatstat.data_1.7-0      
## [117] rmarkdown_2.6             leiden_0.3.6             
## [119] uwot_0.1.10               edgeR_3.32.0             
## [121] DelayedMatrixStats_1.12.0 curl_4.3                 
## [123] shiny_1.5.0               lifecycle_0.2.0          
## [125] nlme_3.1-151              jsonlite_1.7.2           
## [127] BiocNeighbors_1.8.0       viridisLite_0.3.0        
## [129] askpass_1.1               limma_3.46.0             
## [131] pillar_1.4.7              lattice_0.20-41          
## [133] fastmap_1.0.1             httr_1.4.2               
## [135] survival_3.2-7            glue_1.4.2               
## [137] remotes_2.2.0             spatstat_1.64-1          
## [139] png_0.1-7                 bluster_1.0.0            
## [141] bit_4.0.4                 stringi_1.5.3            
## [143] blob_1.2.1                BiocSingular_1.6.0       
## [145] memoise_1.1.0             dplyr_1.0.3              
## [147] irlba_2.3.3               future.apply_1.7.0
```
