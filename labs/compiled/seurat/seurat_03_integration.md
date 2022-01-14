---
title: #INTEG_TITLE:
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
# venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn
# = 1,cexil = 1,lwd=1,col='white',frame=F,borders = NA)

temp <- unique(unlist(hvgs_per_dataset))
overlap <- sapply(hvgs_per_dataset, function(x) {
    temp %in% x
})
pheatmap::pheatmap(t(overlap * 1), cluster_rows = F, color = c("grey90", "grey20"))
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
## 	Retained 1716 anchors
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
## 	Retained 1714 anchors
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
## 	Found 2678 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2223 anchors
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
## 	Found 1824 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1460 anchors
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
## 	Found 2284 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1879 anchors
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
## 	Found 2524 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1979 anchors
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
## 	Retained 1594 anchors
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
## 	Retained 2213 anchors
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
## 	Found 2898 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2173 anchors
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
## 	Found 2919 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2506 anchors
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
## 	Found 2079 anchors
```

```
## Filtering anchors
```

```
## 	Retained 1637 anchors
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
## 	Retained 2161 anchors
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
## 	Found 2819 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2173 anchors
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
## 	Found 2735 anchors
```

```
## Filtering anchors
```

```
## 	Retained 2281 anchors
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
## 	Retained 2769 anchors
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
## 16:13:50 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 16:13:50 Read 5532 rows and found 30 numeric columns
```

```
## 16:13:50 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 16:13:50 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 16:13:51 Writing NN index file to temp file /var/folders/1s/j9ck5c_162s487xcprlxtmdh0000gp/T//RtmpLPXWAN/filed9e1a41c295
## 16:13:51 Searching Annoy index using 1 thread, search_k = 3000
## 16:13:54 Annoy recall = 100%
## 16:13:55 Commencing smooth kNN distance calibration using 1 thread
## 16:13:56 Initializing from normalized Laplacian + noise
## 16:13:57 Commencing optimization for 500 epochs, with 247290 positive edges
## 16:14:10 Optimization finished
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
myfeatures <- c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7",
    "FCGR3A", "CST3", "FCER1A")
plot_list <- list()
for (i in myfeatures) {
    plot_list[[i]] <- FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = i,
        ncol = 3, order = T) + NoLegend() + NoAxes() + NoGrid()
}
plot_grid(ncol = 3, plotlist = plot_list)
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
## Harmony converged after 9 iterations
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
## 16:15:36 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 16:15:36 Read 5532 rows and found 50 numeric columns
```

```
## 16:15:36 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 16:15:36 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 16:15:37 Writing NN index file to temp file /var/folders/1s/j9ck5c_162s487xcprlxtmdh0000gp/T//RtmpLPXWAN/filed9e11df63694
## 16:15:37 Searching Annoy index using 1 thread, search_k = 3000
## 16:15:39 Annoy recall = 100%
## 16:15:39 Commencing smooth kNN distance calibration using 1 thread
## 16:15:41 Initializing from normalized Laplacian + noise
## 16:15:41 Commencing optimization for 500 epochs, with 253252 positive edges
## 16:15:50 Optimization finished
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap_harmony_'
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
## 16:16:08 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 16:16:08 Read 5532 rows and found 100 numeric columns
```

```
## 16:16:08 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 16:16:08 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 16:16:10 Writing NN index file to temp file /var/folders/1s/j9ck5c_162s487xcprlxtmdh0000gp/T//RtmpLPXWAN/filed9e16ccf5ada
## 16:16:10 Searching Annoy index using 1 thread, search_k = 3000
## 16:16:12 Annoy recall = 100%
## 16:16:13 Commencing smooth kNN distance calibration using 1 thread
## 16:16:14 Initializing from normalized Laplacian + noise
## 16:16:15 Commencing optimization for 500 epochs, with 259588 positive edges
## 16:16:28 Optimization finished
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap_scanorama_'
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
##  [1] reticulate_1.22             harmony_1.0                
##  [3] Rcpp_1.0.8                  scran_1.22.0               
##  [5] scuttle_1.4.0               SingleCellExperiment_1.16.0
##  [7] SummarizedExperiment_1.24.0 Biobase_2.54.0             
##  [9] GenomicRanges_1.46.0        GenomeInfoDb_1.30.0        
## [11] IRanges_2.28.0              S4Vectors_0.32.0           
## [13] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
## [15] matrixStats_0.61.0          ggplot2_3.3.5              
## [17] cowplot_1.1.1               KernSmooth_2.23-20         
## [19] fields_13.3                 viridis_0.6.2              
## [21] viridisLite_0.4.0           spam_2.8-0                 
## [23] DoubletFinder_2.0.3         Matrix_1.4-0               
## [25] SeuratObject_4.0.4          Seurat_4.0.6               
## [27] RJSONIO_1.3-1.6             optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                tidyselect_1.1.1         
##   [3] htmlwidgets_1.5.4         grid_4.1.2               
##   [5] BiocParallel_1.28.0       Rtsne_0.15               
##   [7] munsell_0.5.0             ScaledMatrix_1.2.0       
##   [9] codetools_0.2-18          ica_1.0-2                
##  [11] statmod_1.4.36            future_1.23.0            
##  [13] miniUI_0.1.1.1            withr_2.4.3              
##  [15] colorspace_2.0-2          highr_0.9                
##  [17] knitr_1.37                ROCR_1.0-11              
##  [19] tensor_1.5                listenv_0.8.0            
##  [21] labeling_0.4.2            GenomeInfoDbData_1.2.7   
##  [23] polyclip_1.10-0           pheatmap_1.0.12          
##  [25] bit64_4.0.5               farver_2.1.0             
##  [27] rprojroot_2.0.2           parallelly_1.30.0        
##  [29] vctrs_0.3.8               generics_0.1.1           
##  [31] xfun_0.29                 R6_2.5.1                 
##  [33] rsvd_1.0.5                locfit_1.5-9.4           
##  [35] hdf5r_1.3.5               bitops_1.0-7             
##  [37] spatstat.utils_2.3-0      DelayedArray_0.20.0      
##  [39] assertthat_0.2.1          promises_1.2.0.1         
##  [41] scales_1.1.1              gtable_0.3.0             
##  [43] beachmat_2.10.0           globals_0.14.0           
##  [45] processx_3.5.2            goftest_1.2-3            
##  [47] rlang_0.4.12              splines_4.1.2            
##  [49] lazyeval_0.2.2            spatstat.geom_2.3-1      
##  [51] yaml_2.2.1                reshape2_1.4.4           
##  [53] abind_1.4-5               httpuv_1.6.5             
##  [55] tools_4.1.2               ellipsis_0.3.2           
##  [57] spatstat.core_2.3-2       jquerylib_0.1.4          
##  [59] RColorBrewer_1.1-2        ggridges_0.5.3           
##  [61] plyr_1.8.6                sparseMatrixStats_1.6.0  
##  [63] zlibbioc_1.40.0           purrr_0.3.4              
##  [65] RCurl_1.98-1.5            ps_1.6.0                 
##  [67] prettyunits_1.1.1         rpart_4.1-15             
##  [69] deldir_1.0-6              pbapply_1.5-0            
##  [71] zoo_1.8-9                 ggrepel_0.9.1            
##  [73] cluster_2.1.2             here_1.0.1               
##  [75] magrittr_2.0.1            data.table_1.14.2        
##  [77] RSpectra_0.16-0           scattermore_0.7          
##  [79] lmtest_0.9-39             RANN_2.6.1               
##  [81] fitdistrplus_1.1-6        patchwork_1.1.1          
##  [83] mime_0.12                 evaluate_0.14            
##  [85] xtable_1.8-4              gridExtra_2.3            
##  [87] compiler_4.1.2            tibble_3.1.6             
##  [89] maps_3.4.0                crayon_1.4.2             
##  [91] htmltools_0.5.2           mgcv_1.8-38              
##  [93] later_1.2.0               tidyr_1.1.4              
##  [95] DBI_1.1.2                 formatR_1.11             
##  [97] MASS_7.3-55               getopt_1.20.3            
##  [99] cli_3.1.0                 metapod_1.2.0            
## [101] parallel_4.1.2            dotCall64_1.0-1          
## [103] igraph_1.2.11             pkgconfig_2.0.3          
## [105] plotly_4.10.0             spatstat.sparse_2.1-0    
## [107] bslib_0.3.1               dqrng_0.3.0              
## [109] XVector_0.34.0            stringr_1.4.0            
## [111] callr_3.7.0               digest_0.6.29            
## [113] sctransform_0.3.3         RcppAnnoy_0.0.19         
## [115] spatstat.data_2.1-2       rmarkdown_2.11           
## [117] leiden_0.3.9              edgeR_3.36.0             
## [119] uwot_0.1.11               DelayedMatrixStats_1.16.0
## [121] curl_4.3.2                shiny_1.7.1              
## [123] lifecycle_1.0.1           nlme_3.1-155             
## [125] jsonlite_1.7.2            BiocNeighbors_1.12.0     
## [127] limma_3.50.0              fansi_1.0.0              
## [129] pillar_1.6.4              lattice_0.20-45          
## [131] fastmap_1.1.0             httr_1.4.2               
## [133] pkgbuild_1.3.1            survival_3.2-13          
## [135] glue_1.6.0                remotes_2.4.2            
## [137] png_0.1-7                 bluster_1.4.0            
## [139] bit_4.0.4                 stringi_1.7.6            
## [141] sass_0.4.0                BiocSingular_1.10.0      
## [143] dplyr_1.0.7               irlba_2.3.5              
## [145] future.apply_1.8.1
```



