---
title: #INTEG_TITLE:
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


<style>
h1, .h1, h2, .h2, h3, .h3, h4, .h4 { margin-top: 50px }
p.caption {font-size: 0.9em;font-style: italic;color: grey;margin-right: 10%;margin-left: 10%;text-align: justify}
</style>

In this tutorial we will look at different ways of integrating multiple single cell RNA-seq datasets. We will explore two different methods to correct for batch effects across datasets. We will also look at a quantitative measure to assess the quality of the integrated data. Seurat uses the data integration method presented in Comprehensive Integration of Single Cell Data, while Scran and Scanpy use a mutual Nearest neighbour method (MNN). Below you can find a list of the most recent methods for single data integration:

Markdown | Language | Library | Ref
--- | --- | --- | ---
CCA | R | Seurat | [Cell](https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub)
MNN | R/Python | Scater/Scanpy | [Nat. Biotech.](https://www.nature.com/articles/nbt.4091)
Conos | R | conos | [Nat. Methods](https://www.nature.com/articles/s41592-019-0466-z?error=cookies_not_supported&code=5680289b-6edb-40ad-9934-415dac4fdb2f)
Scanorama | Python | scanorama | [Nat. Biotech.](https://www.nature.com/articles/s41587-019-0113-3)

Let's first load necessary libraries and the data saved in the previous lab.


```r
suppressPackageStartupMessages({
    library(scater)
    library(scran)
    library(cowplot)
    library(ggplot2)
    library(rafalib)
    library(venn)
})

sce <- readRDS("data/results/covid_qc_dm.rds")
print(reducedDims(sce))
```

```
## List of length 8
## names(8): PCA UMAP tSNE_on_PCA ... UMAP_on_ScaleData KNN UMAP_on_Graph
```

We split the combined object into a list, with each dataset as an element. We perform standard preprocessing (log-normalization), and identify variable features individually for each dataset based on a variance stabilizing transformation ("vst").


```r
sce.list <- lapply(unique(sce$sample), function(x) {
    x <- sce[, sce$sample == x]
})


mypar(1, 3)
hvgs_per_dataset <- lapply(sce.list, function(x) {
    x <- computeSumFactors(x, sizes = c(20, 40, 60, 80))
    x <- logNormCounts(x)
    var.out <- modelGeneVar(x, method = "loess")
    hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.2), ]
    hvg.out <- hvg.out[order(hvg.out$bio, decreasing = TRUE), ]
    return(rownames(hvg.out))
})
```

```
## Warning in (function (x, sizes, min.mean = NULL, positive = FALSE, scaling =
## NULL) : encountered negative size factor estimates

## Warning in (function (x, sizes, min.mean = NULL, positive = FALSE, scaling =
## NULL) : encountered negative size factor estimates

## Warning in (function (x, sizes, min.mean = NULL, positive = FALSE, scaling =
## NULL) : encountered negative size factor estimates

## Warning in (function (x, sizes, min.mean = NULL, positive = FALSE, scaling =
## NULL) : encountered negative size factor estimates

## Warning in (function (x, sizes, min.mean = NULL, positive = FALSE, scaling =
## NULL) : encountered negative size factor estimates
```

```r
names(hvgs_per_dataset) <- unique(sce$sample)

venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
    cexil = 1, lwd = 1, col = "white", borders = NA)
```

![](scater_03_integration_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

The mutual nearest neighbors (MNN) approach within the scran package utilizes a novel approach to adjust for batch effects. The `fastMNN()` function returns a representation of the data with reduced dimensionality, which can be used in a similar fashion to other lower-dimensional representations such as PCA. In particular, this representation can be used for downstream methods such as clustering. The BNPARAM can be used to specify the specific nearest neighbors method to use from the BiocNeighbors package. Here we make use of the [Annoy library](https://github.com/spotify/annoy) via the `BiocNeighbors::AnnoyParam()` argument. We save the reduced-dimension MNN representation into the reducedDims slot of our sce object.


```r
mnn_out <- batchelor::fastMNN(sce, subset.row = unique(unlist(hvgs_per_dataset)), 
    batch = factor(sce$sample), k = 20, d = 50)
```

**NOTE**: `fastMNN()` does not produce a batch-corrected expression matrix.


```r
mnn_out <- t(reducedDim(mnn_out, "corrected"))
colnames(mnn_out) <- unlist(lapply(sce.list, function(x) {
    colnames(x)
}))
mnn_out <- mnn_out[, colnames(sce)]
rownames(mnn_out) <- paste0("dim", 1:50)
reducedDim(sce, "MNN") <- t(mnn_out)
```

We can observe that a new assay slot is now created under the name `MNN`.


```r
reducedDims(sce)
```

```
## List of length 9
## names(9): PCA UMAP tSNE_on_PCA UMAP_on_PCA ... KNN UMAP_on_Graph MNN
```

Thus, the result from `fastMNN()` should solely be treated as a reduced dimensionality representation, suitable for direct plotting, TSNE/UMAP, clustering, and trajectory analysis that relies on such results.


```r
set.seed(42)
sce <- runTSNE(sce, dimred = "MNN", n_dimred = 50, perplexity = 30, name = "tSNE_on_MNN")
sce <- runUMAP(sce, dimred = "MNN", n_dimred = 50, ncomponents = 2, name = "UMAP_on_MNN")
```

We can now plot the un-integrated and the integrated space reduced dimensions.


```r
plot_grid(ncol = 3,
  plotReducedDim(sce,dimred = "PCA",colour_by = "sample", point_size = 0.6)+ ggplot2::ggtitle(label ="PCA"),
  plotReducedDim(sce,dimred = "tSNE_on_PCA",colour_by = "sample", point_size = 0.6)+ ggplot2::ggtitle(label ="tSNE_on_PCA"),
  plotReducedDim(sce,dimred = "UMAP_on_PCA",colour_by = "sample",point_size = 0.6)+ ggplot2::ggtitle(label ="UMAP_on_PCA"),
  
  plotReducedDim(sce,dimred = "MNN",colour_by = "sample", point_size = 0.6)+ ggplot2::ggtitle(label ="MNN"),
  plotReducedDim(sce,dimred = "tSNE_on_MNN",colour_by = "sample", point_size = 0.6)+ ggplot2::ggtitle(label ="tSNE_on_MNN"),
  plotReducedDim(sce,dimred = "UMAP_on_MNN",colour_by = "sample", point_size = 0.6)+ ggplot2::ggtitle(label ="UMAP_on_MNN")
)
```

![](scater_03_integration_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

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
    plotlist[[i]] <- plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = i, by_exprs_values = "logcounts", 
        point_size = 0.6) + scale_fill_gradientn(colours = colorRampPalette(c("grey90", 
        "orange3", "firebrick", "firebrick", "red", "red"))(10)) + ggtitle(label = i) + 
        theme(plot.title = element_text(size = 20))
}
plot_grid(ncol = 3, plotlist = plotlist)
```

![](scater_03_integration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

#INTEG_R1:

#INTEG_R2:


```r
library(harmony)
```

```
## Loading required package: Rcpp
```

```r
reducedDimNames(sce)

sce <- RunHarmony(sce, group.by.vars = "sample", reduction.save = "harmony", reduction = "PCA", 
    dims.use = 1:50)
```

```
## Harmony 1/10
```

```
## Harmony 2/10
```

```
## Harmony 3/10
```

```
## Harmony 4/10
```

```
## Harmony 5/10
```

```
## Harmony converged after 5 iterations
```

```r
# Here we use all PCs computed from Harmony for UMAP calculation
sce <- runUMAP(sce, dimred = "harmony", n_dimred = 50, ncomponents = 2, name = "UMAP_on_Harmony")
```

```
##  [1] "PCA"               "UMAP"              "tSNE_on_PCA"      
##  [4] "UMAP_on_PCA"       "UMAP10_on_PCA"     "UMAP_on_ScaleData"
##  [7] "KNN"               "UMAP_on_Graph"     "MNN"              
## [10] "tSNE_on_MNN"       "UMAP_on_MNN"
```


#INTEG_R3:

#INTEG_R4:


```r
hvgs <- unique(unlist(hvgs_per_dataset))

scelist <- list()
genelist <- list()
for (i in 1:length(sce.list)) {
    scelist[[i]] <- t(as.matrix(logcounts(sce.list[[i]])[hvgs, ]))
    genelist[[i]] <- hvgs
}

lapply(scelist, dim)
```

```
## [[1]]
## [1] 870 733
## 
## [[2]]
## [1] 581 733
## 
## [[3]]
## [1] 1035  733
## 
## [[4]]
## [1] 1029  733
## 
## [[5]]
## [1] 1136  733
## 
## [[6]]
## [1] 1070  733
```

#INTEG_R5:


```r
library(reticulate)
scanorama <- import("scanorama")

integrated.data <- scanorama$integrate(datasets_full = scelist, genes_list = genelist)

intdimred <- do.call(rbind, integrated.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)
rownames(intdimred) <- colnames(logcounts(sce))

# Add standard deviations in order to draw Elbow Plots in Seurat
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)
attr(intdimred, "varExplained") <- stdevs

reducedDim(sce, "Scanorama_PCA") <- intdimred

# Here we use all PCs computed from Scanorama for UMAP calculation
sce <- runUMAP(sce, dimred = "Scanorama_PCA", n_dimred = 50, ncomponents = 2, name = "UMAP_on_Scanorama")
```

#INTEG_R6:


```r
p1 <- plotReducedDim(sce, dimred = "UMAP_on_PCA", colour_by = "sample", point_size = 0.6) + 
    ggplot2::ggtitle(label = "UMAP_on_PCA")
p2 <- plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = "sample", point_size = 0.6) + 
    ggplot2::ggtitle(label = "UMAP_on_MNN")
p3 <- plotReducedDim(sce, dimred = "UMAP_on_Harmony", colour_by = "sample", point_size = 0.6) + 
    ggplot2::ggtitle(label = "UMAP_on_Harmony")
p4 <- plotReducedDim(sce, dimred = "UMAP_on_Scanorama", colour_by = "sample", point_size = 0.6) + 
    ggplot2::ggtitle(label = "UMAP_on_Scanorama")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + Seurat::NoLegend() + Seurat::NoAxes(), 
    p2 + Seurat::NoLegend() + Seurat::NoAxes(), p3 + Seurat::NoLegend() + Seurat::NoAxes(), 
    p4 + Seurat::NoLegend() + Seurat::NoAxes(), nrow = 2), leg, ncol = 2, widths = c(8, 
    2))
```

![](scater_03_integration_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

#INTEG_R7:


Finally, lets save the integrated data for further analysis.


```r
saveRDS(sce, "data/results/covid_qc_dr_int.rds")
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
##  [1] reticulate_1.18             harmony_1.0                
##  [3] Rcpp_1.0.6                  venn_1.9                   
##  [5] umap_0.2.7.0                rafalib_1.0.0              
##  [7] scDblFinder_1.4.0           biomaRt_2.46.0             
##  [9] DoubletFinder_2.0.3         org.Hs.eg.db_3.12.0        
## [11] AnnotationDbi_1.52.0        cowplot_1.1.1              
## [13] scran_1.18.0                scater_1.18.0              
## [15] ggplot2_3.3.3               SingleCellExperiment_1.12.0
## [17] SummarizedExperiment_1.20.0 Biobase_2.50.0             
## [19] GenomicRanges_1.42.0        GenomeInfoDb_1.26.0        
## [21] IRanges_2.24.0              S4Vectors_0.28.0           
## [23] BiocGenerics_0.36.0         MatrixGenerics_1.2.0       
## [25] matrixStats_0.57.0          RJSONIO_1.3-1.4            
## [27] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] tidyselect_1.1.0          RSQLite_2.2.2            
##   [3] htmlwidgets_1.5.3         grid_4.0.3               
##   [5] BiocParallel_1.24.0       Rtsne_0.15               
##   [7] munsell_0.5.0             codetools_0.2-18         
##   [9] ica_1.0-2                 statmod_1.4.35           
##  [11] xgboost_1.3.0.1           future_1.21.0            
##  [13] miniUI_0.1.1.1            batchelor_1.6.0          
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
##  [59] plyr_1.8.6                sparseMatrixStats_1.2.0  
##  [61] progress_1.2.2            zlibbioc_1.36.0          
##  [63] purrr_0.3.4               RCurl_1.98-1.2           
##  [65] prettyunits_1.1.1         rpart_4.1-15             
##  [67] openssl_1.4.3             deldir_0.2-9             
##  [69] pbapply_1.4-3             viridis_0.5.1            
##  [71] zoo_1.8-8                 ggrepel_0.9.1            
##  [73] cluster_2.1.0             magrittr_2.0.1           
##  [75] data.table_1.13.6         RSpectra_0.16-0          
##  [77] scattermore_0.7           ResidualMatrix_1.0.0     
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
## [109] vipor_0.4.5               admisc_0.11              
## [111] dqrng_0.2.1               XVector_0.30.0           
## [113] stringr_1.4.0             digest_0.6.27            
## [115] sctransform_0.3.2         RcppAnnoy_0.0.18         
## [117] spatstat.data_1.7-0       rmarkdown_2.6            
## [119] leiden_0.3.6              uwot_0.1.10              
## [121] edgeR_3.32.0              DelayedMatrixStats_1.12.0
## [123] curl_4.3                  shiny_1.5.0              
## [125] lifecycle_0.2.0           nlme_3.1-151             
## [127] jsonlite_1.7.2            BiocNeighbors_1.8.0      
## [129] viridisLite_0.3.0         askpass_1.1              
## [131] limma_3.46.0              pillar_1.4.7             
## [133] lattice_0.20-41           fastmap_1.0.1            
## [135] httr_1.4.2                survival_3.2-7           
## [137] glue_1.4.2                remotes_2.2.0            
## [139] spatstat_1.64-1           png_0.1-7                
## [141] bluster_1.0.0             bit_4.0.4                
## [143] stringi_1.5.3             blob_1.2.1               
## [145] BiocSingular_1.6.0        memoise_1.1.0            
## [147] dplyr_1.0.3               irlba_2.3.3              
## [149] future.apply_1.7.0
```

