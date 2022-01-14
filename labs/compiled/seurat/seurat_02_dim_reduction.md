---
title: "Dimensionality reduction"
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
alldata <- FindNeighbors(alldata, reduction = "pca", graph.name = "SNN", assay = "RNA",
    k.param = 20, features = alldata@assays$RNA@var.features)

# Run UMAP on a graph
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph", graph = "SNN", assay = "RNA")
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
##  [1] scran_1.22.0                scuttle_1.4.0              
##  [3] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
##  [5] Biobase_2.54.0              GenomicRanges_1.46.0       
##  [7] GenomeInfoDb_1.30.0         IRanges_2.28.0             
##  [9] S4Vectors_0.32.0            BiocGenerics_0.40.0        
## [11] MatrixGenerics_1.6.0        matrixStats_0.61.0         
## [13] ggplot2_3.3.5               cowplot_1.1.1              
## [15] KernSmooth_2.23-20          fields_13.3                
## [17] viridis_0.6.2               viridisLite_0.4.0          
## [19] spam_2.8-0                  DoubletFinder_2.0.3        
## [21] Matrix_1.4-0                SeuratObject_4.0.4         
## [23] Seurat_4.0.6                RJSONIO_1.3-1.6            
## [25] optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                reticulate_1.22          
##   [3] tidyselect_1.1.1          htmlwidgets_1.5.4        
##   [5] grid_4.1.2                BiocParallel_1.28.0      
##   [7] Rtsne_0.15                munsell_0.5.0            
##   [9] ScaledMatrix_1.2.0        codetools_0.2-18         
##  [11] ica_1.0-2                 statmod_1.4.36           
##  [13] future_1.23.0             miniUI_0.1.1.1           
##  [15] withr_2.4.3               colorspace_2.0-2         
##  [17] highr_0.9                 knitr_1.37               
##  [19] ROCR_1.0-11               tensor_1.5               
##  [21] listenv_0.8.0             labeling_0.4.2           
##  [23] GenomeInfoDbData_1.2.7    polyclip_1.10-0          
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
##  [61] Rcpp_1.0.8                plyr_1.8.6               
##  [63] sparseMatrixStats_1.6.0   zlibbioc_1.40.0          
##  [65] purrr_0.3.4               RCurl_1.98-1.5           
##  [67] ps_1.6.0                  prettyunits_1.1.1        
##  [69] rpart_4.1-15              deldir_1.0-6             
##  [71] pbapply_1.5-0             zoo_1.8-9                
##  [73] ggrepel_0.9.1             cluster_2.1.2            
##  [75] here_1.0.1                magrittr_2.0.1           
##  [77] data.table_1.14.2         RSpectra_0.16-0          
##  [79] scattermore_0.7           lmtest_0.9-39            
##  [81] RANN_2.6.1                fitdistrplus_1.1-6       
##  [83] patchwork_1.1.1           mime_0.12                
##  [85] evaluate_0.14             xtable_1.8-4             
##  [87] gridExtra_2.3             compiler_4.1.2           
##  [89] tibble_3.1.6              maps_3.4.0               
##  [91] crayon_1.4.2              htmltools_0.5.2          
##  [93] mgcv_1.8-38               later_1.2.0              
##  [95] tidyr_1.1.4               DBI_1.1.2                
##  [97] formatR_1.11              MASS_7.3-55              
##  [99] getopt_1.20.3             cli_3.1.0                
## [101] metapod_1.2.0             parallel_4.1.2           
## [103] dotCall64_1.0-1           igraph_1.2.11            
## [105] pkgconfig_2.0.3           plotly_4.10.0            
## [107] spatstat.sparse_2.1-0     bslib_0.3.1              
## [109] dqrng_0.3.0               XVector_0.34.0           
## [111] stringr_1.4.0             callr_3.7.0              
## [113] digest_0.6.29             sctransform_0.3.3        
## [115] RcppAnnoy_0.0.19          spatstat.data_2.1-2      
## [117] rmarkdown_2.11            leiden_0.3.9             
## [119] edgeR_3.36.0              uwot_0.1.11              
## [121] DelayedMatrixStats_1.16.0 curl_4.3.2               
## [123] shiny_1.7.1               lifecycle_1.0.1          
## [125] nlme_3.1-155              jsonlite_1.7.2           
## [127] BiocNeighbors_1.12.0      limma_3.50.0             
## [129] fansi_1.0.0               pillar_1.6.4             
## [131] lattice_0.20-45           fastmap_1.1.0            
## [133] httr_1.4.2                pkgbuild_1.3.1           
## [135] survival_3.2-13           glue_1.6.0               
## [137] remotes_2.4.2             png_0.1-7                
## [139] bluster_1.4.0             bit_4.0.4                
## [141] stringi_1.7.6             sass_0.4.0               
## [143] BiocSingular_1.10.0       dplyr_1.0.7              
## [145] irlba_2.3.5               future.apply_1.8.1
```



