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
  library(Seurat)
  library(cowplot)
  library(ggplot2)
})

alldata <- readRDS("data/3pbmc_qc_dr.rds")
print(names(alldata@reductions))
```

```
## [1] "PCA_on_RNA"        "TSNE_on_RNA"       "UMAP_on_RNA"      
## [4] "UMAP10_on_RNA"     "UMAP_on_ScaleData" "UMAP_on_Graph"
```

We split the combined object into a list, with each dataset as an element. We perform standard preprocessing (log-normalization), and identify variable features individually for each dataset based on a variance stabilizing transformation ("vst").


```r
alldata.list <- SplitObject(alldata, split.by = "orig.ident")
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umap10onrna_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umaponscaledata_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umapongraph_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umap10onrna_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umaponscaledata_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umapongraph_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umap10onrna_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umaponscaledata_'
```

```
## Warning: All object keys must be alphanumeric characters, followed by an
## underscore ('_'), setting key to 'umapongraph_'
```

```r
for (i in 1:length(alldata.list)) {
    alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}

hvgs_per_dataset <- lapply(alldata.list, function(x) { x@assays$RNA@var.features })
venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn = 1,cexil = 1,lwd=1,col="white",frame=F,borders = NA)
```

![](seurat_03_integration_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

We identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.


```r
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30)
```

```
## Computing 2000 integration features
```

```
## Scaling features for provided objects
```

```
## Finding all pairwise anchors
```

```
## Running CCA
```

```
## Merging objects
```

```
## Finding neighborhoods
```

```
## Finding anchors
```

```
## 	Found 2185 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1906 anchors
```

```
## Extracting within-dataset neighbors
```

```
## Running CCA
```

```
## Merging objects
```

```
## Finding neighborhoods
```

```
## Finding anchors
```

```
## 	Found 2279 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1979 anchors
```

```
## Extracting within-dataset neighbors
```

```
## Running CCA
```

```
## Merging objects
```

```
## Finding neighborhoods
```

```
## Finding anchors
```

```
## 	Found 3062 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2681 anchors
```

```
## Extracting within-dataset neighbors
```

We then pass these anchors to the IntegrateData function, which returns a Seurat object.


```r
alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
```

```
## Merging dataset 1 into 3
```

```
## Extracting anchors for merged samples
```

```
## Finding integration vectors
```

```
## Finding integration vector weights
```

```
## Integrating data
```

```
## Merging dataset 2 into 3 1
```

```
## Extracting anchors for merged samples
```

```
## Finding integration vectors
```

```
## Finding integration vector weights
```

```
## Integrating data
```

We can observe that a new assay slot is now created under the name `CCA`.


```r
names(alldata.int@assays)
```

```
## [1] "RNA" "CCA"
```

After running IntegrateData, the Seurat object will contain a new Assay with the integrated (or ‘batch-corrected’) expression matrix. Note that the original (uncorrected values) are still stored in the object in the “RNA” assay, so you can switch back and forth. We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP and TSNE. The integrated datasets cluster by cell type, instead of by technology.


```r
#Run Dimensionality reduction on integrated space
alldata.int <- ScaleData(alldata.int, verbose = FALSE,assay = "CCA")
alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE, assay = "CCA",reduction.name = "PCA_on_CCA")
alldata.int <- RunUMAP(alldata.int, reduction = "PCA_on_CCA", dims = 1:30,reduction.name = "UMAP_on_CCA")
alldata.int <- RunTSNE(alldata.int, reduction = "PCA_on_CCA", dims = 1:30,reduction.name = "TSNE_on_CCA")
```

We can now plot the un-integrated and the integrated space reduced dimensions.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "PCA_on_RNA", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "TSNE_on_RNA", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "UMAP_on_RNA", group.by = "orig.ident"),
  
  DimPlot(alldata.int, reduction = "PCA_on_CCA", group.by = "orig.ident"),
  DimPlot(alldata.int, reduction = "TSNE_on_CCA", group.by = "orig.ident"),
  DimPlot(alldata.int, reduction = "UMAP_on_CCA", group.by = "orig.ident")
)
```

![](seurat_03_integration_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

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
FeaturePlot(alldata.int, reduction = "UMAP_on_CCA",dims = 1:2,features = c("CD3E","CD4","CD8A","NKG7","GNLY","MS4A1","CD14","LYZ","MS4A7","FCGR3A","CST3","FCER1A"),ncol = 4,order = T)
```

![](seurat_03_integration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Finally, lets save the integrated data for further analysis.


```r
saveRDS(alldata.int,"data/3pbmc_qc_dr_int.rds")
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
##  [1] scran_1.10.1                SingleCellExperiment_1.4.0 
##  [3] SummarizedExperiment_1.12.0 DelayedArray_0.8.0         
##  [5] matrixStats_0.55.0          Biobase_2.42.0             
##  [7] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1        
##  [9] IRanges_2.16.0              S4Vectors_0.20.1           
## [11] BiocGenerics_0.28.0         BiocParallel_1.16.6        
## [13] ggplot2_3.2.1               cowplot_1.0.0              
## [15] Matrix_1.2-17               Seurat_3.0.1               
## [17] RJSONIO_1.3-1.2             optparse_1.6.4             
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



