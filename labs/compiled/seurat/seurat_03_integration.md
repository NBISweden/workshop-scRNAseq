---
title: #INTEG_TITLE:
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 28, 2021'
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
## 17:36:42 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 17:36:42 Read 5532 rows and found 30 numeric columns
```

```
## 17:36:42 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 17:36:42 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 17:36:43 Writing NN index file to temp file /var/folders/n0/1679kqxs6s1bbdhj59hgpq0rm04rx6/T//Rtmpei6lqo/file6b0b7710496a
## 17:36:43 Searching Annoy index using 1 thread, search_k = 3000
## 17:36:44 Annoy recall = 100%
## 17:36:44 Commencing smooth kNN distance calibration using 1 thread
## 17:36:45 Initializing from normalized Laplacian + noise
## 17:36:46 Commencing optimization for 500 epochs, with 254060 positive edges
## 17:36:53 Optimization finished
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
## 17:37:17 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 17:37:17 Read 5532 rows and found 50 numeric columns
```

```
## 17:37:17 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 17:37:17 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 17:37:18 Writing NN index file to temp file /var/folders/n0/1679kqxs6s1bbdhj59hgpq0rm04rx6/T//Rtmpei6lqo/file6b0b6eaacd00
## 17:37:18 Searching Annoy index using 1 thread, search_k = 3000
## 17:37:19 Annoy recall = 100%
## 17:37:19 Commencing smooth kNN distance calibration using 1 thread
## 17:37:21 Initializing from normalized Laplacian + noise
## 17:37:21 Commencing optimization for 500 epochs, with 253382 positive edges
## 17:37:28 Optimization finished
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
## 17:37:36 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 17:37:36 Read 5532 rows and found 100 numeric columns
```

```
## 17:37:36 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 17:37:36 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 17:37:37 Writing NN index file to temp file /var/folders/n0/1679kqxs6s1bbdhj59hgpq0rm04rx6/T//Rtmpei6lqo/file6b0b6f22ac62
## 17:37:37 Searching Annoy index using 1 thread, search_k = 3000
## 17:37:38 Annoy recall = 100%
## 17:37:38 Commencing smooth kNN distance calibration using 1 thread
## 17:37:40 Initializing from normalized Laplacian + noise
## 17:37:40 Commencing optimization for 500 epochs, with 259588 positive edges
## 17:37:47 Optimization finished
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
##  [1] parallel  stats4    grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] reticulate_1.18             harmony_1.0                
##  [3] Rcpp_1.0.6                  scran_1.18.0               
##  [5] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
##  [7] Biobase_2.50.0              GenomicRanges_1.42.0       
##  [9] GenomeInfoDb_1.26.0         IRanges_2.24.0             
## [11] S4Vectors_0.28.0            BiocGenerics_0.36.0        
## [13] MatrixGenerics_1.2.0        matrixStats_0.57.0         
## [15] ggplot2_3.3.3               cowplot_1.1.1              
## [17] KernSmooth_2.23-18          fields_11.6                
## [19] spam_2.6-0                  dotCall64_1.0-0            
## [21] DoubletFinder_2.0.3         Matrix_1.3-2               
## [23] Seurat_3.2.3                RJSONIO_1.3-1.4            
## [25] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.6                igraph_1.2.6             
##   [3] lazyeval_0.2.2            splines_4.0.3            
##   [5] BiocParallel_1.24.0       listenv_0.8.0            
##   [7] scattermore_0.7           digest_0.6.27            
##   [9] htmltools_0.5.1           magrittr_2.0.1           
##  [11] tensor_1.5                cluster_2.1.0            
##  [13] ROCR_1.0-11               limma_3.46.0             
##  [15] remotes_2.2.0             globals_0.14.0           
##  [17] colorspace_2.0-0          ggrepel_0.9.1            
##  [19] xfun_0.20                 dplyr_1.0.3              
##  [21] crayon_1.3.4              RCurl_1.98-1.2           
##  [23] jsonlite_1.7.2            spatstat_1.64-1          
##  [25] spatstat.data_1.7-0       survival_3.2-7           
##  [27] zoo_1.8-8                 glue_1.4.2               
##  [29] polyclip_1.10-0           gtable_0.3.0             
##  [31] zlibbioc_1.36.0           XVector_0.30.0           
##  [33] leiden_0.3.6              DelayedArray_0.16.0      
##  [35] BiocSingular_1.6.0        future.apply_1.7.0       
##  [37] maps_3.3.0                abind_1.4-5              
##  [39] scales_1.1.1              edgeR_3.32.0             
##  [41] DBI_1.1.1                 miniUI_0.1.1.1           
##  [43] viridisLite_0.3.0         xtable_1.8-4             
##  [45] dqrng_0.2.1               bit_4.0.4                
##  [47] rsvd_1.0.3                htmlwidgets_1.5.3        
##  [49] httr_1.4.2                getopt_1.20.3            
##  [51] RColorBrewer_1.1-2        ellipsis_0.3.1           
##  [53] ica_1.0-2                 scuttle_1.0.0            
##  [55] pkgconfig_2.0.3           farver_2.0.3             
##  [57] uwot_0.1.10               deldir_0.2-9             
##  [59] locfit_1.5-9.4            tidyselect_1.1.0         
##  [61] labeling_0.4.2            rlang_0.4.10             
##  [63] reshape2_1.4.4            later_1.1.0.1            
##  [65] munsell_0.5.0             tools_4.0.3              
##  [67] generics_0.1.0            ggridges_0.5.3           
##  [69] evaluate_0.14             stringr_1.4.0            
##  [71] fastmap_1.0.1             yaml_2.2.1               
##  [73] goftest_1.2-2             knitr_1.30               
##  [75] bit64_4.0.5               fitdistrplus_1.1-3       
##  [77] admisc_0.11               purrr_0.3.4              
##  [79] RANN_2.6.1                sparseMatrixStats_1.2.0  
##  [81] pbapply_1.4-3             future_1.21.0            
##  [83] nlme_3.1-151              mime_0.9                 
##  [85] formatR_1.7               venn_1.9                 
##  [87] hdf5r_1.3.3               compiler_4.0.3           
##  [89] plotly_4.9.3              curl_4.3                 
##  [91] png_0.1-7                 spatstat.utils_1.20-2    
##  [93] statmod_1.4.35            tibble_3.0.5             
##  [95] stringi_1.5.3             RSpectra_0.16-0          
##  [97] bluster_1.0.0             lattice_0.20-41          
##  [99] vctrs_0.3.6               pillar_1.4.7             
## [101] lifecycle_0.2.0           lmtest_0.9-38            
## [103] BiocNeighbors_1.8.0       RcppAnnoy_0.0.18         
## [105] data.table_1.13.6         bitops_1.0-6             
## [107] irlba_2.3.3               httpuv_1.5.5             
## [109] patchwork_1.1.1           R6_2.5.0                 
## [111] promises_1.1.1            gridExtra_2.3            
## [113] parallelly_1.23.0         codetools_0.2-18         
## [115] MASS_7.3-53               assertthat_0.2.1         
## [117] withr_2.4.0               sctransform_0.3.2        
## [119] GenomeInfoDbData_1.2.4    mgcv_1.8-33              
## [121] beachmat_2.6.0            rpart_4.1-15             
## [123] tidyr_1.1.2               DelayedMatrixStats_1.12.0
## [125] rmarkdown_2.6             Rtsne_0.15               
## [127] shiny_1.5.0
```



