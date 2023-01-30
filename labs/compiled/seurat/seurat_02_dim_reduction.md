---
title: "Dimensionality reduction"
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 27, 2023'
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
    library(Seurat)
    library(cowplot)
    library(ggplot2)
    library(scran)
})

alldata <- readRDS("data/results/seurat_covid_qc.rds")
```

### Feature selection

Next, we first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.


```r
suppressWarnings(suppressMessages(alldata <- FindVariableFeatures(alldata, selection.method = "vst",
    nfeatures = 2000, verbose = FALSE, assay = "RNA")))
top20 <- head(VariableFeatures(alldata), 20)

LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

### Z-score transformation

Now that the data is prepared, we now proceed with PCA. Since each gene has a different expression level, it means that genes with higher expression values will naturally have higher variation that will be captured by PCA. This means that we need to somehow give each gene a similar weight when performing PCA (see below). The common practice is to center and scale each gene before performing PCA. This exact scaling is called Z-score normalization it is very useful for PCA, clustering and plotting heatmaps. <br>Additionally, we can use regression to remove any unwanted sources of variation from the dataset, such as `cell cycle`, `sequencing depth`, `percent mitocondria`. This is achieved by doing a generalized linear regression using these parameters as covariates in the model. Then the residuals of the model are taken as the "regressed data". Although perhaps not in the best way, batch effect regression can also be done here.


```r
alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"),
    assay = "RNA")
```


## PCA
***

Performing PCA has many useful applications and interpretations, which much depends on the data used. In the case of life sciences, we want to segregate samples based on gene expression patterns in the data.


```r
alldata <- RunPCA(alldata, npcs = 50, verbose = F)
```

We then plot the first principal components.


```r
plot_grid(ncol = 3, DimPlot(alldata, reduction = "pca", group.by = "orig.ident",
    dims = 1:2), DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 3:4),
    DimPlot(alldata, reduction = "pca", group.by = "orig.ident", dims = 5:6))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

To identify which genes (Seurat) or metadata paramters (Scater/Scran) contribute the most to each PC, one can retreive the loading matrix information. Unfortunatelly this is not implemented in Scater/Scran, so you will need to compute PCA using `logcounts`.


```r
VizDimLoadings(alldata, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


We can also plot the amount of variance explained by each PC.


```r
ElbowPlot(alldata, reduction = "pca", ndims = 50)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Based on this plot, we can see that the top 8 PCs retain a lot of information, while other PCs contain pregressivelly less. However, it is still advisable to use more PCs since they might contain informaktion about rare cell types (such as platelets and DCs in this dataset)

## tSNE
***

We can now run [BH-tSNE](https://arxiv.org/abs/1301.3342).


```r
alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30, perplexity = 30, max_iter = 1000,
    theta = 0.5, eta = 200, num_threads = 0)
# see ?Rtsne and ?RunTSNE for more info
```

We can now plot the tSNE colored per dataset. We can clearly see the effect of batches present in the dataset.


```r
plot_grid(ncol = 3, DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


***
## UMAP
***

We can now run [UMAP](https://arxiv.org/abs/1802.03426) for cell embeddings.


```r
alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30, n.components = 2, n.neighbors = 30,
    n.epochs = 200, min.dist = 0.3, learning.rate = 1, spread = 1)
# see ?RunUMAP for more info
```

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

We have now done Variable gene selection, PCA and UMAP with the settings we chose. Test a few different ways of selecting variable genes, number of PCs for UMAP and check how it influences your embedding.
</div>

Another usefullness of UMAP is that it is not limitted by the number of dimensions the data cen be reduced into (unlike tSNE). We can simply reduce the dimentions altering the `n.components` parameter.


```r
# we can add in additional reductions, by defulat they are named 'pca', 'umap',
# 'tsne' etc. But we can specify alternative names with reduction.name

alldata <- RunUMAP(alldata, reduction.name = "UMAP10_on_PCA", reduction = "pca",
    dims = 1:30, n.components = 10, n.neighbors = 30, n.epochs = 200, min.dist = 0.3,
    learning.rate = 1, spread = 1)
# see ?RunUMAP for more info
```

We can now plot the UMAP colored per dataset. Although less distinct as in the tSNE, we still see quite an effect of the different batches in the data.


```r
plot_grid(ncol = 3, DimPlot(alldata, reduction = "umap", group.by = "orig.ident") +
    ggplot2::ggtitle(label = "UMAP_on_PCA"), DimPlot(alldata, reduction = "UMAP10_on_PCA",
    group.by = "orig.ident", dims = 1:2) + ggplot2::ggtitle(label = "UMAP10_on_PCA"),
    DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident", dims = 3:4) +
        ggplot2::ggtitle(label = "UMAP10_on_PCA"))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

We can now plot PCA, UMAP and tSNE side by side for comparison. Here, we can conclude that our dataset contains a batch effect that needs to be corrected before proceeding to clustering and differential gene expression analysis.


```r
plot_grid(ncol = 3, DimPlot(alldata, reduction = "pca", group.by = "orig.ident"),
    DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"), DimPlot(alldata,
        reduction = "umap", group.by = "orig.ident"))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

## Using ScaledData and graphs for DR
***

Althought running a sencond dimmensionality reduction (i.e tSNE or UMAP) on PCA would be a standard approach (because it allows higher computation efficiency), the options are actually limiteless. Below we will show a couple of other common options such as running directly on the scaled data (which was used for PCA) or on a graph built from scaled data. We will show from now on only UMAP, but the same applies for tSNE.

### Using ScaledData for UMAP

To run tSNE or UMAP on the scaled data, one firts needs to select the number of variables to use. This is because including dimentions that do contribute to the separation of your cell types will in the end mask those differences. Another reason for it is because running with all genes/features also will take longer or might be computationally unfeasible. Therefore we will use the scaled data of the highly variable genes.


```r
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_ScaleData", features = alldata@assays$RNA@var.features,
    assay = "RNA", n.components = 2, n.neighbors = 30, n.epochs = 200, min.dist = 0.3,
    learning.rate = 1, spread = 1)
```

### Using a Graph for UMAP

To run tSNE or UMAP on the a graph, we first need to build a graph from the data. In fact, both tSNE and UMAP first build a graph from the data using a specified distance metrix and then optimize the embedding. Since a graph is just a matrix containing distances from cell to cell and as such, you can run either UMAP or tSNE using any other distance metric desired. Euclidean and Correlation are ususally the most commonly used.


```r
# Build Graph
alldata <- FindNeighbors(alldata, reduction = "pca", assay = "RNA", k.param = 20,
    features = alldata@assays$RNA@var.features)

# Run UMAP on a graph
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph", umap.method = "umap-learn",
    graph = "RNA_snn", assay = "RNA")
```

We can now plot the UMAP comparing both on PCA vs ScaledSata vs Graph.


```r
p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_PCA")
p2 <- DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident") +
    ggplot2::ggtitle(label = "UMAP_on_ScaleData")
p3 <- DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_Graph")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() +
    NoAxes(), p3 + NoLegend() + NoAxes(), leg, nrow = 2), ncol = 1, widths = c(1))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

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
myfeatures <- c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7",
    "FCGR3A", "CST3", "FCER1A")
plot_list <- list()
for (i in myfeatures) {
    plot_list[[i]] <- FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = i,
        ncol = 3, order = T) + NoLegend() + NoAxes() + NoGrid()
}
plot_grid(ncol = 3, plotlist = plot_list)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


We can finally save the object for use in future steps.


```r
saveRDS(alldata, "data/results/covid_qc_dr.rds")
```


<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

Select some of your dimensionality reductions and plot some of the QC stats that were calculated in the previous lab. Can you see if some of the separation in your data is driven by quality of the cells?
</div>


### Session Info
***


```r
sessionInfo()
```

```
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Big Sur/Monterey 10.16
## 
## Matrix products: default
## BLAS/LAPACK: /Users/asabjor/miniconda3/envs/scRNAseq2023/lib/libopenblasp-r0.3.21.dylib
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] scran_1.22.1                scuttle_1.4.0              
##  [3] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
##  [5] Biobase_2.54.0              GenomicRanges_1.46.1       
##  [7] GenomeInfoDb_1.30.1         IRanges_2.28.0             
##  [9] S4Vectors_0.32.4            BiocGenerics_0.40.0        
## [11] MatrixGenerics_1.6.0        matrixStats_0.63.0         
## [13] ggplot2_3.4.0               cowplot_1.1.1              
## [15] SeuratObject_4.1.3          Seurat_4.3.0               
## [17] RJSONIO_1.3-1.7             optparse_1.7.3             
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.8                igraph_1.3.5             
##   [3] lazyeval_0.2.2            sp_1.6-0                 
##   [5] splines_4.1.3             BiocParallel_1.28.3      
##   [7] listenv_0.9.0             scattermore_0.8          
##   [9] digest_0.6.31             htmltools_0.5.4          
##  [11] fansi_1.0.4               magrittr_2.0.3           
##  [13] ScaledMatrix_1.2.0        tensor_1.5               
##  [15] cluster_2.1.4             ROCR_1.0-11              
##  [17] limma_3.50.3              globals_0.16.2           
##  [19] spatstat.sparse_3.0-0     colorspace_2.1-0         
##  [21] ggrepel_0.9.2             xfun_0.36                
##  [23] dplyr_1.0.10              RCurl_1.98-1.9           
##  [25] jsonlite_1.8.4            progressr_0.13.0         
##  [27] spatstat.data_3.0-0       survival_3.5-0           
##  [29] zoo_1.8-11                glue_1.6.2               
##  [31] polyclip_1.10-4           gtable_0.3.1             
##  [33] zlibbioc_1.40.0           XVector_0.34.0           
##  [35] leiden_0.4.3              DelayedArray_0.20.0      
##  [37] BiocSingular_1.10.0       future.apply_1.10.0      
##  [39] abind_1.4-5               scales_1.2.1             
##  [41] edgeR_3.36.0              DBI_1.1.3                
##  [43] spatstat.random_3.0-1     miniUI_0.1.1.1           
##  [45] Rcpp_1.0.10               viridisLite_0.4.1        
##  [47] xtable_1.8-4              dqrng_0.3.0              
##  [49] reticulate_1.27           rsvd_1.0.5               
##  [51] metapod_1.2.0             htmlwidgets_1.6.1        
##  [53] httr_1.4.4                getopt_1.20.3            
##  [55] RColorBrewer_1.1-3        ellipsis_0.3.2           
##  [57] ica_1.0-3                 farver_2.1.1             
##  [59] pkgconfig_2.0.3           sass_0.4.5               
##  [61] uwot_0.1.14               deldir_1.0-6             
##  [63] here_1.0.1                locfit_1.5-9.7           
##  [65] utf8_1.2.2                labeling_0.4.2           
##  [67] tidyselect_1.2.0          rlang_1.0.6              
##  [69] reshape2_1.4.4            later_1.3.0              
##  [71] munsell_0.5.0             tools_4.1.3              
##  [73] cachem_1.0.6              cli_3.6.0                
##  [75] generics_0.1.3            ggridges_0.5.4           
##  [77] evaluate_0.20             stringr_1.5.0            
##  [79] fastmap_1.1.0             yaml_2.3.7               
##  [81] goftest_1.2-3             knitr_1.41               
##  [83] fitdistrplus_1.1-8        purrr_1.0.1              
##  [85] RANN_2.6.1                sparseMatrixStats_1.6.0  
##  [87] pbapply_1.7-0             future_1.30.0            
##  [89] nlme_3.1-161              mime_0.12                
##  [91] formatR_1.14              compiler_4.1.3           
##  [93] plotly_4.10.1             png_0.1-8                
##  [95] spatstat.utils_3.0-1      statmod_1.5.0            
##  [97] tibble_3.1.8              bslib_0.4.2              
##  [99] stringi_1.7.12            highr_0.10               
## [101] bluster_1.4.0             lattice_0.20-45          
## [103] Matrix_1.5-3              vctrs_0.5.2              
## [105] pillar_1.8.1              lifecycle_1.0.3          
## [107] spatstat.geom_3.0-5       lmtest_0.9-40            
## [109] jquerylib_0.1.4           BiocNeighbors_1.12.0     
## [111] RcppAnnoy_0.0.20          data.table_1.14.6        
## [113] bitops_1.0-7              irlba_2.3.5.1            
## [115] httpuv_1.6.8              patchwork_1.1.2          
## [117] R6_2.5.1                  promises_1.2.0.1         
## [119] KernSmooth_2.23-20        gridExtra_2.3            
## [121] parallelly_1.34.0         codetools_0.2-18         
## [123] MASS_7.3-58.2             assertthat_0.2.1         
## [125] rprojroot_2.0.3           withr_2.5.0              
## [127] sctransform_0.3.5         GenomeInfoDbData_1.2.7   
## [129] parallel_4.1.3            beachmat_2.10.0          
## [131] grid_4.1.3                tidyr_1.2.1              
## [133] DelayedMatrixStats_1.16.0 rmarkdown_2.20           
## [135] Rtsne_0.16                spatstat.explore_3.0-5   
## [137] shiny_1.7.4
```



