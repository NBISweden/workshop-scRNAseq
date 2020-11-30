---
title: #INTEG_TITLE:
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'November 27, 2020'
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

alldata <- readRDS("data/results/covid_qc_dr.rds")
print(names(alldata@reductions))
```

```
## [1] "pca"               "umap"              "tsne"             
## [4] "UMAP10_on_PCA"     "UMAP_on_ScaleData" "UMAP_on_Graph"
```

We split the combined object into a list, with each dataset as an element. We perform standard preprocessing (log-normalization), and identify variable features individually for each dataset based on a variance stabilizing transformation ("vst").


```r
alldata.list <- SplitObject(alldata, split.by = "orig.ident")
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap10_on_pca_ to umap10onpca_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_scaledata_ to umaponscaledata_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_graph_ to umapongraph_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap10_on_pca_ to umap10onpca_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_scaledata_ to umaponscaledata_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_graph_ to umapongraph_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap10_on_pca_ to umap10onpca_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_scaledata_ to umaponscaledata_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_graph_ to umapongraph_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap10_on_pca_ to umap10onpca_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_scaledata_ to umaponscaledata_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_graph_ to umapongraph_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap10_on_pca_ to umap10onpca_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_scaledata_ to umaponscaledata_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_graph_ to umapongraph_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap10_on_pca_ to umap10onpca_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_scaledata_ to umaponscaledata_
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_on_graph_ to umapongraph_
```

```r
for (i in 1:length(alldata.list)) {
    alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

hvgs_per_dataset <- lapply(alldata.list, function(x) {
    x@assays$RNA@var.features
})
venn::venn(hvgs_per_dataset, opacity = 0.4, zcolor = (scales::hue_pal())(3), cexsn = 1, 
    cexil = 1, lwd = 1, col = "white", frame = F, borders = NA)
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
## 	Found 2005 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1718 anchors
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
## 	Found 2139 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1711 anchors
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
## 	Found 2679 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2227 anchors
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
## 	Found 1823 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1463 anchors
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
## 	Found 2286 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1894 anchors
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
## 	Found 2527 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1972 anchors
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
## 	Found 2104 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1622 anchors
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
## 	Found 2630 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2226 anchors
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
## 	Found 2899 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2171 anchors
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
## 	Found 2920 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2513 anchors
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
## 	Found 2081 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1647 anchors
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
## 	Found 2523 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2166 anchors
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
## 	Found 2821 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2174 anchors
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
## 	Found 2734 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2283 anchors
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
## 	Found 3066 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2773 anchors
```

We then pass these anchors to the IntegrateData function, which returns a Seurat object.


```r
alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
```

```
## Merging dataset 1 into 2
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
## Merging dataset 6 into 5
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
## Merging dataset 3 into 2 1
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
## Merging dataset 4 into 5 6
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
## Merging dataset 2 1 3 into 5 6 4
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
## Warning: Adding a command log without an assay associated with it
```

We can observe that a new assay slot is now created under the name `CCA`.


```r
names(alldata.int@assays)

# by default, Seurat now sets the integrated assay as the default assay, so any
# operation you now perform will be on the ingegrated data.

alldata.int@active.assay
```

```
## [1] "RNA" "CCA"
## [1] "CCA"
```

After running IntegrateData, the Seurat object will contain a new Assay with the integrated (or ‘batch-corrected’) expression matrix. Note that the original (uncorrected values) are still stored in the object in the “RNA” assay, so you can switch back and forth. We can then use this new integrated matrix for downstream analysis and visualization. Here we scale the integrated data, run PCA, and visualize the results with UMAP and TSNE. The integrated datasets cluster by cell type, instead of by technology.


```r
# Run Dimensionality reduction on integrated space
alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)
alldata.int <- RunUMAP(alldata.int, dims = 1:30)
```

```
## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
## This message will be shown once per session
```

```
## 16:59:13 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## 16:59:13 Read 5532 rows and found 30 numeric columns
```

```
## 16:59:13 Using Annoy for neighbor search, n_neighbors = 30
```

```
## 16:59:13 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 16:59:14 Writing NN index file to temp file /var/folders/1s/j9ck5c_162s487xcprlxtmdh0000gp/T//RtmpJo37qJ/file89311ad8ed9d
## 16:59:14 Searching Annoy index using 1 thread, search_k = 3000
## 16:59:15 Annoy recall = 100%
## 16:59:15 Commencing smooth kNN distance calibration using 1 thread
## 16:59:16 Initializing from normalized Laplacian + noise
## 16:59:16 Commencing optimization for 500 epochs, with 254078 positive edges
## 16:59:22 Optimization finished
```

```r
alldata.int <- RunTSNE(alldata.int, dims = 1:30)
```

We can now plot the un-integrated and the integrated space reduced dimensions.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident"),
  
  DimPlot(alldata.int, reduction = "pca", group.by = "orig.ident"),
  DimPlot(alldata.int, reduction = "tsne", group.by = "orig.ident"),
  DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident")
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
FeaturePlot(alldata.int, reduction = "umap", features = c("CD3E", "CD4", "CD8A", 
    "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7", "FCGR3A", "CST3", "FCER1A"), 
    ncol = 4, order = T)
```

![](seurat_03_integration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Finally, lets save the integrated data for further analysis.


```r
saveRDS(alldata.int, "data/results/covid_qc_dr_int.rds")
```


### Session Info
***


```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS/LAPACK: /Users/asbj/miniconda3/envs/scRNAseq2021/lib/libopenblasp-r0.3.12.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] ggplot2_3.3.2   cowplot_1.1.0   Seurat_3.2.2    RJSONIO_1.3-1.4
## [5] optparse_1.6.6 
## 
## loaded via a namespace (and not attached):
##   [1] nlme_3.1-150          matrixStats_0.57.0    RcppAnnoy_0.0.17     
##   [4] RColorBrewer_1.1-2    httr_1.4.2            sctransform_0.3.1    
##   [7] tools_4.0.3           R6_2.5.0              irlba_2.3.3          
##  [10] rpart_4.1-15          KernSmooth_2.23-18    uwot_0.1.9           
##  [13] mgcv_1.8-33           lazyeval_0.2.2        colorspace_2.0-0     
##  [16] withr_2.3.0           tidyselect_1.1.0      gridExtra_2.3        
##  [19] compiler_4.0.3        formatR_1.7           plotly_4.9.2.1       
##  [22] labeling_0.4.2        scales_1.1.1          spatstat.data_1.5-2  
##  [25] lmtest_0.9-38         ggridges_0.5.2        pbapply_1.4-3        
##  [28] goftest_1.2-2         spatstat_1.64-1       stringr_1.4.0        
##  [31] digest_0.6.27         spatstat.utils_1.17-0 rmarkdown_2.5        
##  [34] pkgconfig_2.0.3       htmltools_0.5.0       parallelly_1.21.0    
##  [37] fastmap_1.0.1         htmlwidgets_1.5.2     rlang_0.4.8          
##  [40] shiny_1.5.0           farver_2.0.3          generics_0.1.0       
##  [43] zoo_1.8-8             jsonlite_1.7.1        ica_1.0-2            
##  [46] dplyr_1.0.2           magrittr_2.0.1        patchwork_1.1.0      
##  [49] Matrix_1.2-18         Rcpp_1.0.5            munsell_0.5.0        
##  [52] abind_1.4-5           reticulate_1.18       lifecycle_0.2.0      
##  [55] stringi_1.5.3         yaml_2.2.1            MASS_7.3-53          
##  [58] Rtsne_0.15            plyr_1.8.6            grid_4.0.3           
##  [61] parallel_4.0.3        listenv_0.8.0         promises_1.1.1       
##  [64] ggrepel_0.8.2         venn_1.9              crayon_1.3.4         
##  [67] deldir_0.2-3          miniUI_0.1.1.1        lattice_0.20-41      
##  [70] splines_4.0.3         tensor_1.5            knitr_1.30           
##  [73] pillar_1.4.7          igraph_1.2.6          admisc_0.11          
##  [76] future.apply_1.6.0    reshape2_1.4.4        codetools_0.2-18     
##  [79] leiden_0.3.5          glue_1.4.2            evaluate_0.14        
##  [82] data.table_1.13.2     vctrs_0.3.5           png_0.1-7            
##  [85] httpuv_1.5.4          polyclip_1.10-0       gtable_0.3.0         
##  [88] getopt_1.20.3         RANN_2.6.1            purrr_0.3.4          
##  [91] tidyr_1.1.2           future_1.20.1         xfun_0.19            
##  [94] rsvd_1.0.3            mime_0.9              xtable_1.8-4         
##  [97] RSpectra_0.16-0       later_1.1.0.1         survival_3.2-7       
## [100] viridisLite_0.3.0     tibble_3.0.4          cluster_2.1.0        
## [103] globals_0.14.0        fitdistrplus_1.1-1    ellipsis_0.3.1       
## [106] ROCR_1.0-11
```



