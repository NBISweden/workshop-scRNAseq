---
title: #INTEG_TITLE:
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

sce <- readRDS("data/3pbmc_qc_dm.rds")
print(names(sce@reducedDims))
```

```
## [1] "PCA"               "tSNE_on_PCA"       "UMAP_on_PCA"      
## [4] "UMAP10_on_PCA"     "UMAP_on_ScaleData" "KNN"              
## [7] "UMAP_on_Graph"
```

We split the combined object into a list, with each dataset as an element. We perform standard preprocessing (log-normalization), and identify variable features individually for each dataset based on a variance stabilizing transformation ("vst").


```r
sce.list <- lapply( unique(sce$sample_id), function(x){
  x <- sce[ , sce$sample_id == x ] })


mypar(1,3)
hvgs_per_dataset <- lapply( sce.list, function(x){
  x <- computeSumFactors(x, sizes=c(20, 40, 60, 80))
  x <- normalize(x)
  var.fit <- trendVar(x, use.spikes=FALSE,method="loess",loess.args=list(span=0.05))
  var.out <- decomposeVar(x, var.fit)
  hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.2),]
  hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
  return(rownames(hvg.out))
})
```

```
## Warning in FUN(...): encountered negative size factor estimates
```

```r
names(hvgs_per_dataset) <- unique(sce$sample_id)

venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn = 1,cexil = 1,lwd=1,col="white",borders = NA)
```

![](scater_03_integration_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

The mutual nearest neighbors (MNN) approach within the scran package utilizes a novel approach to adjust for batch effects. The `fastMNN()` function returns a representation of the data with reduced dimensionality, which can be used in a similar fashion to other lower-dimensional representations such as PCA. In particular, this representation can be used for downstream methods such as clustering. The BNPARAM can be used to specify the specific nearest neighbors method to use from the BiocNeighbors package. Here we make use of the [Annoy library](https://github.com/spotify/annoy) via the `BiocNeighbors::AnnoyParam()` argument. We save the reduced-dimension MNN representation into the reducedDims slot of our sce object.


```r
mnn_out <- fastMNN(sce.list[[1]], sce.list[[2]], sce.list[[3]],
                   subset.row = unique(unlist(hvgs_per_dataset)),
                   k = 20, d = 50, approximate = TRUE,
                   # BPPARAM = BiocParallel::MulticoreParam(4),
                   BNPARAM = BiocNeighbors::AnnoyParam())
```

**NOTE**: `fastMNN()` does not produce a batch-corrected expression matrix. 


```r
mnn_out <- t(mnn_out$corrected)
colnames(mnn_out) <- unlist(lapply(sce.list,function(x){colnames(x)}))
mnn_out <- mnn_out[,colnames(sce)]
rownames(mnn_out) <- paste0("dim",1:50)
reducedDim(sce, "MNN") <- t(mnn_out)
```

We can observe that a new assay slot is now created under the name `MNN`.


```r
names(sce@reducedDims)
```

```
## [1] "PCA"               "tSNE_on_PCA"       "UMAP_on_PCA"      
## [4] "UMAP10_on_PCA"     "UMAP_on_ScaleData" "KNN"              
## [7] "UMAP_on_Graph"     "MNN"
```

Thus, the result from `fastMNN()` should solely be treated as a reduced dimensionality representation, suitable for direct plotting, TSNE/UMAP, clustering, and trajectory analysis that relies on such results.


```r
set.seed(42)
sce <- runTSNE(sce, use_dimred = "MNN", n_dimred = 50, perplexity = 30)
reducedDimNames(sce)[reducedDimNames(sce)=="TSNE"] <- "tSNE_on_MNN"


sce <- runUMAP(sce,use_dimred = "MNN", n_dimred = 50, ncomponents = 2)
reducedDimNames(sce)[reducedDimNames(sce)=="UMAP"] <- "UMAP_on_MNN"
```

We can now plot the un-integrated and the integrated space reduced dimensions.


```r
plot_grid(ncol = 3,
  plotReducedDim(sce,use_dimred = "PCA",colour_by = "sample_id",add_ticks = F, point_size = 0.6)+ ggplot2::ggtitle(label ="PCA"),
  plotReducedDim(sce,use_dimred = "tSNE_on_PCA",colour_by = "sample_id",add_ticks = F, point_size = 0.6)+ ggplot2::ggtitle(label ="tSNE_on_PCA"),
  plotReducedDim(sce,use_dimred = "UMAP_on_PCA",colour_by = "sample_id",add_ticks = F, point_size = 0.6)+ ggplot2::ggtitle(label ="UMAP_on_PCA"),
  
  plotReducedDim(sce,use_dimred = "MNN",colour_by = "sample_id",add_ticks = F, point_size = 0.6)+ ggplot2::ggtitle(label ="MNN"),
  plotReducedDim(sce,use_dimred = "tSNE_on_MNN",colour_by = "sample_id",add_ticks = F, point_size = 0.6)+ ggplot2::ggtitle(label ="tSNE_on_MNN"),
  plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = "sample_id",add_ticks = F, point_size = 0.6)+ ggplot2::ggtitle(label ="UMAP_on_MNN")
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
for(i in c("CD3E","CD4","CD8A","NKG7","GNLY","MS4A1","CD14","LYZ","MS4A7","FCGR3A","CST3","FCER1A")){
  plotlist[[i]] <- plotReducedDim(sce,use_dimred = "UMAP_on_MNN",colour_by = i,by_exprs_values = "logcounts",add_ticks = F, point_size = 0.6) +
  scale_fill_gradientn(colours = colorRampPalette(c("grey90","orange3","firebrick","firebrick","red","red" ))(10)) +
  ggtitle(label = i)+ theme(plot.title = element_text(size=20)) }
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
plot_grid(ncol=3, plotlist = plotlist)
```

![](scater_03_integration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Finally, lets save the integrated data for further analysis.


```r
saveRDS(sce,"data/3pbmc_qc_dr_int.rds")
```

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
##  [1] venn_1.7                    rafalib_1.0.0              
##  [3] cowplot_1.0.0               scran_1.10.1               
##  [5] scater_1.10.1               ggplot2_3.2.1              
##  [7] SingleCellExperiment_1.4.0  SummarizedExperiment_1.12.0
##  [9] DelayedArray_0.8.0          BiocParallel_1.16.6        
## [11] matrixStats_0.55.0          Biobase_2.42.0             
## [13] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1        
## [15] IRanges_2.16.0              S4Vectors_0.20.1           
## [17] BiocGenerics_0.28.0         RJSONIO_1.3-1.2            
## [19] optparse_1.6.4             
## 
## loaded via a namespace (and not attached):
##  [1] viridis_0.5.1            dynamicTreeCut_1.63-1    edgeR_3.24.3            
##  [4] jsonlite_1.6             viridisLite_0.3.0        DelayedMatrixStats_1.4.0
##  [7] assertthat_0.2.1         statmod_1.4.32           GenomeInfoDbData_1.2.0  
## [10] vipor_0.4.5              yaml_2.2.0               pillar_1.4.2            
## [13] lattice_0.20-38          reticulate_1.13          glue_1.3.1              
## [16] limma_3.38.3             digest_0.6.23            RColorBrewer_1.1-2      
## [19] XVector_0.22.0           colorspace_1.4-1         htmltools_0.4.0         
## [22] Matrix_1.2-17            plyr_1.8.4               pkgconfig_2.0.3         
## [25] zlibbioc_1.28.0          purrr_0.3.3              scales_1.0.0            
## [28] RSpectra_0.15-0          HDF5Array_1.10.1         Rtsne_0.15              
## [31] getopt_1.20.3            openssl_1.1              tibble_2.1.3            
## [34] umap_0.2.3.1             withr_2.1.2              lazyeval_0.2.2          
## [37] magrittr_1.5             crayon_1.3.4             evaluate_0.14           
## [40] beeswarm_0.2.3           tools_3.5.1              stringr_1.4.0           
## [43] Rhdf5lib_1.4.3           munsell_0.5.0            locfit_1.5-9.1          
## [46] irlba_2.3.3              compiler_3.5.1           rlang_0.4.2             
## [49] rhdf5_2.26.2             grid_3.5.1               RCurl_1.95-4.12         
## [52] BiocNeighbors_1.0.0      igraph_1.2.4.1           labeling_0.3            
## [55] bitops_1.0-6             rmarkdown_1.17           gtable_0.3.0            
## [58] reshape2_1.4.3           R6_2.4.1                 gridExtra_2.3           
## [61] knitr_1.26               dplyr_0.8.3              stringi_1.4.3           
## [64] ggbeeswarm_0.6.0         Rcpp_1.0.3               tidyselect_0.2.5        
## [67] xfun_0.11
```

