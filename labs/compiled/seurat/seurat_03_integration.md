---
title: #INTEG_TITLE:
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'December 22, 2021'
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
## 14:52:30 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 14:52:30 Read 5532 rows and found 30 numeric columns
```

```
## 14:52:30 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 14:52:30 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 14:52:31 Writing NN index file to temp file /var/folders/n0/1679kqxs6s1bbdhj59hgpq0rm04rx6/T//RtmpLopeeM/file161db19090a57
## 14:52:31 Searching Annoy index using 1 thread, search_k = 3000
## 14:52:32 Annoy recall = 100%
## 14:52:33 Commencing smooth kNN distance calibration using 1 thread
## 14:52:34 Initializing from normalized Laplacian + noise
## 14:52:34 Commencing optimization for 500 epochs, with 247290 positive edges
## 14:52:42 Optimization finished
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
## 14:53:08 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 14:53:08 Read 5532 rows and found 50 numeric columns
```

```
## 14:53:08 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 14:53:08 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 14:53:09 Writing NN index file to temp file /var/folders/n0/1679kqxs6s1bbdhj59hgpq0rm04rx6/T//RtmpLopeeM/file161db18f6481f
## 14:53:09 Searching Annoy index using 1 thread, search_k = 3000
## 14:53:10 Annoy recall = 100%
## 14:53:11 Commencing smooth kNN distance calibration using 1 thread
## 14:53:12 Initializing from normalized Laplacian + noise
## 14:53:12 Commencing optimization for 500 epochs, with 253252 positive edges
## 14:53:20 Optimization finished
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
## 14:53:32 UMAP embedding parameters a = 0.9922 b = 1.112
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 14:53:32 Read 5532 rows and found 100 numeric columns
```

```
## 14:53:32 Using Annoy for neighbor search, n_neighbors = 30
```

```
## Found more than one class "dist" in cache; using the first, from namespace 'spam'
```

```
## Also defined by 'BiocGenerics'
```

```
## 14:53:32 Building Annoy index with metric = cosine, n_trees = 50
```

```
## 0%   10   20   30   40   50   60   70   80   90   100%
```

```
## [----|----|----|----|----|----|----|----|----|----|
```

```
## **************************************************|
## 14:53:32 Writing NN index file to temp file /var/folders/n0/1679kqxs6s1bbdhj59hgpq0rm04rx6/T//RtmpLopeeM/file161db623d7e2e
## 14:53:32 Searching Annoy index using 1 thread, search_k = 3000
## 14:53:34 Annoy recall = 100%
## 14:53:34 Commencing smooth kNN distance calibration using 1 thread
## 14:53:35 Initializing from normalized Laplacian + noise
## 14:53:36 Commencing optimization for 500 epochs, with 259588 positive edges
## 14:53:44 Optimization finished
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
## R version 4.0.5 (2021-03-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS/LAPACK: /Users/paulo.czarnewski/miniconda3/envs/scRNAseq2022/lib/libopenblasp-r0.3.18.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
##  [1] parallel  stats4    grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] reticulate_1.22             harmony_1.0                
##  [3] Rcpp_1.0.7                  scran_1.18.5               
##  [5] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
##  [7] Biobase_2.50.0              GenomicRanges_1.42.0       
##  [9] GenomeInfoDb_1.26.4         IRanges_2.24.1             
## [11] S4Vectors_0.28.1            BiocGenerics_0.36.0        
## [13] MatrixGenerics_1.2.1        matrixStats_0.61.0         
## [15] ggplot2_3.3.5               cowplot_1.1.1              
## [17] KernSmooth_2.23-20          fields_13.3                
## [19] viridis_0.6.2               viridisLite_0.4.0          
## [21] spam_2.7-0                  dotCall64_1.0-1            
## [23] DoubletFinder_2.0.3         Matrix_1.4-0               
## [25] SeuratObject_4.0.4          Seurat_4.0.6               
## [27] RJSONIO_1.3-1.6             optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                tidyselect_1.1.1         
##   [3] htmlwidgets_1.5.4         BiocParallel_1.24.1      
##   [5] Rtsne_0.15                munsell_0.5.0            
##   [7] codetools_0.2-18          ica_1.0-2                
##   [9] statmod_1.4.36            future_1.23.0            
##  [11] miniUI_0.1.1.1            withr_2.4.3              
##  [13] colorspace_2.0-2          highr_0.9                
##  [15] knitr_1.37                ROCR_1.0-11              
##  [17] tensor_1.5                listenv_0.8.0            
##  [19] labeling_0.4.2            GenomeInfoDbData_1.2.4   
##  [21] polyclip_1.10-0           pheatmap_1.0.12          
##  [23] bit64_4.0.5               farver_2.1.0             
##  [25] rprojroot_2.0.2           parallelly_1.30.0        
##  [27] vctrs_0.3.8               generics_0.1.1           
##  [29] xfun_0.29                 R6_2.5.1                 
##  [31] rsvd_1.0.5                locfit_1.5-9.4           
##  [33] hdf5r_1.3.5               bitops_1.0-7             
##  [35] spatstat.utils_2.3-0      DelayedArray_0.16.3      
##  [37] assertthat_0.2.1          promises_1.2.0.1         
##  [39] scales_1.1.1              gtable_0.3.0             
##  [41] beachmat_2.6.4            globals_0.14.0           
##  [43] processx_3.5.2            goftest_1.2-3            
##  [45] rlang_0.4.12              splines_4.0.5            
##  [47] lazyeval_0.2.2            spatstat.geom_2.3-1      
##  [49] yaml_2.2.1                reshape2_1.4.4           
##  [51] abind_1.4-5               httpuv_1.6.4             
##  [53] tools_4.0.5               ellipsis_0.3.2           
##  [55] spatstat.core_2.3-2       jquerylib_0.1.4          
##  [57] RColorBrewer_1.1-2        ggridges_0.5.3           
##  [59] plyr_1.8.6                sparseMatrixStats_1.2.1  
##  [61] zlibbioc_1.36.0           purrr_0.3.4              
##  [63] RCurl_1.98-1.5            ps_1.6.0                 
##  [65] prettyunits_1.1.1         rpart_4.1-15             
##  [67] deldir_1.0-6              pbapply_1.5-0            
##  [69] zoo_1.8-9                 ggrepel_0.9.1            
##  [71] cluster_2.1.2             here_1.0.1               
##  [73] magrittr_2.0.1            data.table_1.14.2        
##  [75] RSpectra_0.16-0           scattermore_0.7          
##  [77] lmtest_0.9-39             RANN_2.6.1               
##  [79] fitdistrplus_1.1-6        patchwork_1.1.1          
##  [81] mime_0.12                 evaluate_0.14            
##  [83] xtable_1.8-4              gridExtra_2.3            
##  [85] compiler_4.0.5            tibble_3.1.6             
##  [87] maps_3.4.0                crayon_1.4.2             
##  [89] htmltools_0.5.2           mgcv_1.8-38              
##  [91] later_1.2.0               tidyr_1.1.4              
##  [93] DBI_1.1.2                 formatR_1.11             
##  [95] MASS_7.3-54               getopt_1.20.3            
##  [97] cli_3.1.0                 igraph_1.2.10            
##  [99] pkgconfig_2.0.3           scuttle_1.0.4            
## [101] plotly_4.10.0             spatstat.sparse_2.1-0    
## [103] bslib_0.3.1               dqrng_0.3.0              
## [105] XVector_0.30.0            stringr_1.4.0            
## [107] callr_3.7.0               digest_0.6.29            
## [109] sctransform_0.3.2         RcppAnnoy_0.0.19         
## [111] spatstat.data_2.1-2       rmarkdown_2.11           
## [113] leiden_0.3.9              edgeR_3.32.1             
## [115] uwot_0.1.11               DelayedMatrixStats_1.12.3
## [117] curl_4.3.2                shiny_1.7.1              
## [119] lifecycle_1.0.1           nlme_3.1-153             
## [121] jsonlite_1.7.2            BiocNeighbors_1.8.2      
## [123] limma_3.46.0              fansi_0.5.0              
## [125] pillar_1.6.4              lattice_0.20-45          
## [127] fastmap_1.1.0             httr_1.4.2               
## [129] pkgbuild_1.3.1            survival_3.2-13          
## [131] glue_1.6.0                remotes_2.4.2            
## [133] png_0.1-7                 bluster_1.0.0            
## [135] bit_4.0.4                 stringi_1.7.6            
## [137] sass_0.4.0                BiocSingular_1.6.0       
## [139] dplyr_1.0.7               irlba_2.3.5              
## [141] future.apply_1.8.1
```



