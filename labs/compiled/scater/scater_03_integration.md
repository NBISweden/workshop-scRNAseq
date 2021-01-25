---
title: #INTEG_TITLE:
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 22, 2021'
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
```

```
## Warning: package 'venn' was built under R version 3.6.3
```

```r
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
## Warning in FUN(...): encountered negative size factor estimates

## Warning in FUN(...): encountered negative size factor estimates

## Warning in FUN(...): encountered negative size factor estimates

## Warning in FUN(...): encountered negative size factor estimates
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
```

```
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
## Scale for 'fill' is already present. Adding another scale for 'fill', which
## will replace the existing scale.
```

```r
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

```
## Warning: package 'Rcpp' was built under R version 3.6.3
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
## [1] 901 444
## 
## [[2]]
## [1] 598 444
## 
## [[3]]
## [1] 1052  444
## 
## [[4]]
## [1] 1062  444
## 
## [[5]]
## [1] 1175  444
## 
## [[6]]
## [1] 1108  444
```

#INTEG_R5:


```r
library(reticulate)
```

```
## Warning: package 'reticulate' was built under R version 3.6.3
```

```r
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
## R version 3.6.1 (2019-07-05)
## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
## Running under: Ubuntu 20.04 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /home/czarnewski/miniconda3/envs/scRNAseq2021/lib/libopenblasp-r0.3.10.so
## 
## locale:
##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] reticulate_1.18             harmony_1.0                
##  [3] Rcpp_1.0.6                  venn_1.9                   
##  [5] umap_0.2.7.0                rafalib_1.0.0              
##  [7] scDblFinder_1.1.8           DoubletFinder_2.0.3        
##  [9] org.Hs.eg.db_3.10.0         AnnotationDbi_1.48.0       
## [11] cowplot_1.1.1               scran_1.14.1               
## [13] scater_1.14.0               ggplot2_3.3.3              
## [15] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
## [17] DelayedArray_0.12.0         BiocParallel_1.20.0        
## [19] matrixStats_0.57.0          Biobase_2.46.0             
## [21] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [23] IRanges_2.20.0              S4Vectors_0.24.0           
## [25] BiocGenerics_0.32.0         RJSONIO_1.3-1.4            
## [27] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.6               igraph_1.2.6             lazyeval_0.2.2          
##   [4] splines_3.6.1            listenv_0.8.0            scattermore_0.7         
##   [7] digest_0.6.27            htmltools_0.5.1          viridis_0.5.1           
##  [10] magrittr_2.0.1           memoise_1.1.0            tensor_1.5              
##  [13] cluster_2.1.0            ROCR_1.0-11              limma_3.42.0            
##  [16] remotes_2.2.0            globals_0.14.0           askpass_1.1             
##  [19] colorspace_2.0-0         rappdirs_0.3.1           blob_1.2.1              
##  [22] ggrepel_0.9.1            xfun_0.20                dplyr_1.0.3             
##  [25] crayon_1.3.4             RCurl_1.98-1.2           jsonlite_1.7.2          
##  [28] spatstat_1.64-1          spatstat.data_1.7-0      survival_3.2-7          
##  [31] zoo_1.8-8                glue_1.4.2               polyclip_1.10-0         
##  [34] gtable_0.3.0             zlibbioc_1.32.0          XVector_0.26.0          
##  [37] leiden_0.3.6             BiocSingular_1.2.0       future.apply_1.7.0      
##  [40] abind_1.4-5              scales_1.1.1             DBI_1.1.1               
##  [43] edgeR_3.28.0             miniUI_0.1.1.1           viridisLite_0.3.0       
##  [46] xtable_1.8-4             dqrng_0.2.1              bit_4.0.4               
##  [49] rsvd_1.0.3               htmlwidgets_1.5.3        httr_1.4.2              
##  [52] getopt_1.20.3            RColorBrewer_1.1-2       ellipsis_0.3.1          
##  [55] Seurat_3.2.3             ica_1.0-2                farver_2.0.3            
##  [58] pkgconfig_2.0.3          uwot_0.1.10              deldir_0.2-3            
##  [61] locfit_1.5-9.4           labeling_0.4.2           tidyselect_1.1.0        
##  [64] rlang_0.4.10             reshape2_1.4.4           later_1.1.0.1           
##  [67] munsell_0.5.0            tools_3.6.1              generics_0.1.0          
##  [70] RSQLite_2.2.2            ggridges_0.5.3           batchelor_1.2.1         
##  [73] evaluate_0.14            stringr_1.4.0            fastmap_1.0.1           
##  [76] goftest_1.2-2            yaml_2.2.1               knitr_1.30              
##  [79] bit64_4.0.5              fitdistrplus_1.1-3       admisc_0.11             
##  [82] randomForest_4.6-14      purrr_0.3.4              RANN_2.6.1              
##  [85] nlme_3.1-150             pbapply_1.4-3            future_1.21.0           
##  [88] mime_0.9                 formatR_1.7              hdf5r_1.3.3             
##  [91] compiler_3.6.1           beeswarm_0.2.3           plotly_4.9.3            
##  [94] curl_4.3                 png_0.1-7                spatstat.utils_1.20-2   
##  [97] tibble_3.0.5             statmod_1.4.35           stringi_1.5.3           
## [100] RSpectra_0.16-0          lattice_0.20-41          Matrix_1.3-2            
## [103] vctrs_0.3.6              pillar_1.4.7             lifecycle_0.2.0         
## [106] BiocManager_1.30.10      lmtest_0.9-38            RcppAnnoy_0.0.18        
## [109] BiocNeighbors_1.4.0      data.table_1.13.6        bitops_1.0-6            
## [112] irlba_2.3.3              httpuv_1.5.5             patchwork_1.1.1         
## [115] R6_2.5.0                 promises_1.1.1           KernSmooth_2.23-18      
## [118] gridExtra_2.3            vipor_0.4.5              parallelly_1.23.0       
## [121] codetools_0.2-18         MASS_7.3-53              assertthat_0.2.1        
## [124] openssl_1.4.3            withr_2.4.0              sctransform_0.3.2       
## [127] GenomeInfoDbData_1.2.2   mgcv_1.8-33              rpart_4.1-15            
## [130] grid_3.6.1               tidyr_1.1.2              rmarkdown_2.6           
## [133] DelayedMatrixStats_1.8.0 Rtsne_0.15               shiny_1.5.0             
## [136] ggbeeswarm_0.6.0
```

