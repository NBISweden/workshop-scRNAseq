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
alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30, 
    reduction = "cca")
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
## 	Retained 1713 anchors
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
## 	Found 2138 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1695 anchors
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
## 	Found 2689 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2240 anchors
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
## 	Found 1820 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1474 anchors
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
## 	Found 2285 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1884 anchors
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
## 	Retained 1984 anchors
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
## 	Found 2110 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1627 anchors
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
## 	Found 2629 anchors
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
## 	Found 2902 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2185 anchors
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
## 	Found 2900 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2481 anchors
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
## 	Found 2075 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1639 anchors
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
## 	Found 2520 anchors
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
## 	Found 2815 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2162 anchors
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
## 	Found 2745 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2294 anchors
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
## 	Found 3058 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2763 anchors
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
## 10:12:06 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 10:12:06 Read 5532 rows and found 30 numeric columns
```

```
## 10:12:06 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 10:12:06 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 10:12:07 Writing NN index file to temp file /tmp/Rtmp0xrzXB/file90310ab489f
## 10:12:07 Searching Annoy index using 1 thread, search_k = 3000
## 10:12:09 Annoy recall = 100%
## 10:12:09 Commencing smooth kNN distance calibration using 1 thread
## 10:12:11 Initializing from normalized Laplacian + noise
## 10:12:11 Commencing optimization for 500 epochs, with 254036 positive edges
## 10:12:21 Optimization finished
```

```r
alldata.int <- RunTSNE(alldata.int, dims = 1:30)
```

We can now plot the un-integrated and the integrated space reduced dimensions.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
  
  DimPlot(alldata.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
  DimPlot(alldata.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
  DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
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
    order = T, slot = "data", combine = T)
```

![](seurat_03_integration_files/figure-html/unnamed-chunk-8-1.png)<!-- -->



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
alldata.harmony <- RunHarmony(alldata, group.by.vars = "orig.ident", reduction = "pca", 
    dims.use = 1:50, assay.use = "RNA")
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
## Harmony 6/10
```

```
## Harmony 7/10
```

```
## Harmony 8/10
```

```
## Harmony 9/10
```

```
## Harmony 10/10
```

```
## Warning: Invalid name supplied, making object name syntactically valid. New
## object name is Seurat..ProjectDim.RNA.harmony; see ?make.names for more details
## on syntax validity
```

```r
# Here we use all PCs computed from Harmony for UMAP calculation
alldata.int[["harmony"]] <- alldata.harmony[["harmony"]]
alldata.int <- RunUMAP(alldata.int, dims = 1:50, reduction = "harmony", reduction.name = "umap_harmony")
```

```
## 10:13:09 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 10:13:09 Read 5532 rows and found 50 numeric columns
```

```
## 10:13:09 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 10:13:09 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 10:13:10 Writing NN index file to temp file /tmp/Rtmp0xrzXB/file9032059d7e1
## 10:13:10 Searching Annoy index using 1 thread, search_k = 3000
## 10:13:12 Annoy recall = 100%
## 10:13:13 Commencing smooth kNN distance calibration using 1 thread
## 10:13:14 Initializing from normalized Laplacian + noise
## 10:13:14 Commencing optimization for 500 epochs, with 253364 positive edges
## 10:13:24 Optimization finished
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap_harmony_'
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_harmony_ to umapharmony_
```



```r
hvgs <- unique(unlist(hvgs_per_dataset))

assaylist <- list()
genelist <- list()
for (i in 1:length(alldata.list)) {
    assaylist[[i]] <- t(as.matrix(GetAssayData(alldata.list[[i]], "data")[hvgs, ]))
    genelist[[i]] <- hvgs
}

lapply(assaylist, dim)
```

```
## [[1]]
## [1]  540 5203
## 
## [[2]]
## [1]  837 5203
## 
## [[3]]
## [1]  976 5203
## 
## [[4]]
## [1] 1022 5203
## 
## [[5]]
## [1] 1129 5203
## 
## [[6]]
## [1] 1028 5203
```




```r
library(reticulate)
```

```
## Warning: package 'reticulate' was built under R version 3.6.3
```

```r
scanorama <- import("scanorama")

integrated.data <- scanorama$integrate(datasets_full = assaylist, genes_list = genelist)

intdimred <- do.call(rbind, integrated.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)
rownames(intdimred) <- colnames(alldata.int)

# Add standard deviations in order to draw Elbow Plots in Seurat
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

alldata.int[["scanorama"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, 
    key = "PC_", assay = "RNA")
```

```
## Warning: Cannot add objects with duplicate keys (offending key: PC_), setting
## key to 'scanorama_'
```

```r
# Here we use all PCs computed from Scanorama for UMAP calculation
alldata.int <- RunUMAP(alldata.int, dims = 1:100, reduction = "scanorama", reduction.name = "umap_scanorama")
```

```
## 10:13:45 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 10:13:45 Read 5532 rows and found 100 numeric columns
```

```
## 10:13:45 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'BiocGenerics'
```

```
## Also defined by 'spam'
```

```
## 10:13:45 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 10:13:46 Writing NN index file to temp file /tmp/Rtmp0xrzXB/file903458276d5
## 10:13:46 Searching Annoy index using 1 thread, search_k = 3000
## 10:13:47 Annoy recall = 100%
## 10:13:48 Commencing smooth kNN distance calibration using 1 thread
## 10:13:49 Initializing from normalized Laplacian + noise
## 10:13:49 Commencing optimization for 500 epochs, with 259588 positive edges
## 10:13:59 Optimization finished
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap_scanorama_'
```

```
## Warning: Keys should be one or more alphanumeric characters followed by an
## underscore, setting key from umap_scanorama_ to umapscanorama_
```



```r
p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident") + ggtitle("UMAP raw_data")
p2 <- DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident") + ggtitle("UMAP CCA")
p3 <- DimPlot(alldata.int, reduction = "umap_harmony", group.by = "orig.ident") + 
    ggtitle("UMAP Harmony")
p4 <- DimPlot(alldata.int, reduction = "umap_scanorama", group.by = "orig.ident") + 
    ggtitle("UMAP Scanorama")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() + 
    NoAxes(), p3 + NoLegend() + NoAxes(), p4 + NoLegend() + NoAxes(), nrow = 2), 
    leg, ncol = 2, widths = c(8, 2))
```

![](seurat_03_integration_files/figure-html/unnamed-chunk-12-1.png)<!-- -->




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
##  [1] parallel  stats4    grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] reticulate_1.18             harmony_1.0                
##  [3] Rcpp_1.0.6                  scran_1.14.1               
##  [5] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
##  [7] DelayedArray_0.12.0         BiocParallel_1.20.0        
##  [9] matrixStats_0.57.0          Biobase_2.46.0             
## [11] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [13] IRanges_2.20.0              S4Vectors_0.24.0           
## [15] BiocGenerics_0.32.0         ggplot2_3.3.3              
## [17] cowplot_1.1.1               KernSmooth_2.23-18         
## [19] fields_11.6                 spam_2.5-1                 
## [21] dotCall64_1.0-0             biomaRt_2.42.1             
## [23] DoubletFinder_2.0.3         Matrix_1.3-2               
## [25] Seurat_3.2.3                RJSONIO_1.3-1.4            
## [27] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] tidyselect_1.1.0         RSQLite_2.2.2            AnnotationDbi_1.48.0    
##   [4] htmlwidgets_1.5.3        Rtsne_0.15               munsell_0.5.0           
##   [7] codetools_0.2-18         ica_1.0-2                statmod_1.4.35          
##  [10] future_1.21.0            miniUI_0.1.1.1           withr_2.4.0             
##  [13] colorspace_2.0-0         knitr_1.30               ROCR_1.0-11             
##  [16] tensor_1.5               listenv_0.8.0            labeling_0.4.2          
##  [19] GenomeInfoDbData_1.2.2   polyclip_1.10-0          bit64_4.0.5             
##  [22] farver_2.0.3             parallelly_1.23.0        vctrs_0.3.6             
##  [25] generics_0.1.0           xfun_0.20                BiocFileCache_1.10.0    
##  [28] R6_2.5.0                 ggbeeswarm_0.6.0         rsvd_1.0.3              
##  [31] locfit_1.5-9.4           hdf5r_1.3.3              bitops_1.0-6            
##  [34] spatstat.utils_1.20-2    assertthat_0.2.1         promises_1.1.1          
##  [37] scales_1.1.1             beeswarm_0.2.3           gtable_0.3.0            
##  [40] globals_0.14.0           goftest_1.2-2            rlang_0.4.10            
##  [43] splines_3.6.1            lazyeval_0.2.2           yaml_2.2.1              
##  [46] reshape2_1.4.4           abind_1.4-5              httpuv_1.5.5            
##  [49] tools_3.6.1              ellipsis_0.3.1           RColorBrewer_1.1-2      
##  [52] ggridges_0.5.3           plyr_1.8.6               progress_1.2.2          
##  [55] zlibbioc_1.32.0          purrr_0.3.4              RCurl_1.98-1.2          
##  [58] prettyunits_1.1.1        rpart_4.1-15             openssl_1.4.3           
##  [61] deldir_0.2-3             pbapply_1.4-3            viridis_0.5.1           
##  [64] zoo_1.8-8                ggrepel_0.9.1            cluster_2.1.0           
##  [67] magrittr_2.0.1           data.table_1.13.6        RSpectra_0.16-0         
##  [70] scattermore_0.7          lmtest_0.9-38            RANN_2.6.1              
##  [73] fitdistrplus_1.1-3       hms_1.0.0                patchwork_1.1.1         
##  [76] mime_0.9                 evaluate_0.14            xtable_1.8-4            
##  [79] XML_3.99-0.3             gridExtra_2.3            compiler_3.6.1          
##  [82] scater_1.14.0            tibble_3.0.5             maps_3.3.0              
##  [85] crayon_1.3.4             htmltools_0.5.1          venn_1.9                
##  [88] mgcv_1.8-33              later_1.1.0.1            tidyr_1.1.2             
##  [91] DBI_1.1.1                formatR_1.7              dbplyr_2.0.0            
##  [94] MASS_7.3-53              rappdirs_0.3.1           getopt_1.20.3           
##  [97] igraph_1.2.6             pkgconfig_2.0.3          plotly_4.9.3            
## [100] vipor_0.4.5              admisc_0.11              dqrng_0.2.1             
## [103] XVector_0.26.0           stringr_1.4.0            digest_0.6.27           
## [106] sctransform_0.3.2        RcppAnnoy_0.0.18         spatstat.data_1.7-0     
## [109] rmarkdown_2.6            leiden_0.3.6             uwot_0.1.10             
## [112] edgeR_3.28.0             DelayedMatrixStats_1.8.0 curl_4.3                
## [115] shiny_1.5.0              lifecycle_0.2.0          nlme_3.1-150            
## [118] jsonlite_1.7.2           BiocNeighbors_1.4.0      viridisLite_0.3.0       
## [121] askpass_1.1              limma_3.42.0             pillar_1.4.7            
## [124] lattice_0.20-41          fastmap_1.0.1            httr_1.4.2              
## [127] survival_3.2-7           glue_1.4.2               remotes_2.2.0           
## [130] spatstat_1.64-1          png_0.1-7                bit_4.0.4               
## [133] stringi_1.5.3            blob_1.2.1               BiocSingular_1.2.0      
## [136] memoise_1.1.0            dplyr_1.0.3              irlba_2.3.3             
## [139] future.apply_1.7.0
```



