---
title: "Dimensionality reduction"
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
### TASK: Test settings

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
FeaturePlot(alldata, reduction = "umap", dims = 1:2, features = myfeatures, ncol = 3, 
    order = T)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
for (i in myfeatures) {
    
}
p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_PCA")
p2 <- DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident") + 
    ggplot2::ggtitle(label = "UMAP_on_ScaleData")
p3 <- DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident") + ggplot2::ggtitle(label = "UMAP_on_Graph")
leg <- get_legend(p1)

gridExtra::grid.arrange(gridExtra::arrangeGrob(p1 + NoLegend() + NoAxes(), p2 + NoLegend() + 
    NoAxes(), p3 + NoLegend() + NoAxes(), leg, nrow = 2), ncol = 1, widths = c(1))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-17-2.png)<!-- -->


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
##  [1] scran_1.18.0                SingleCellExperiment_1.12.0
##  [3] SummarizedExperiment_1.20.0 Biobase_2.50.0             
##  [5] GenomicRanges_1.42.0        GenomeInfoDb_1.26.0        
##  [7] IRanges_2.24.0              S4Vectors_0.28.0           
##  [9] BiocGenerics_0.36.0         MatrixGenerics_1.2.0       
## [11] matrixStats_0.57.0          ggplot2_3.3.3              
## [13] cowplot_1.1.1               KernSmooth_2.23-18         
## [15] fields_11.6                 spam_2.6-0                 
## [17] dotCall64_1.0-0             DoubletFinder_2.0.3        
## [19] Matrix_1.3-2                Seurat_3.2.3               
## [21] RJSONIO_1.3-1.4             optparse_1.6.6             
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
##  [43] Rcpp_1.0.6                viridisLite_0.3.0        
##  [45] xtable_1.8-4              dqrng_0.2.1              
##  [47] reticulate_1.18           bit_4.0.4                
##  [49] rsvd_1.0.3                htmlwidgets_1.5.3        
##  [51] httr_1.4.2                getopt_1.20.3            
##  [53] RColorBrewer_1.1-2        ellipsis_0.3.1           
##  [55] ica_1.0-2                 scuttle_1.0.0            
##  [57] pkgconfig_2.0.3           farver_2.0.3             
##  [59] uwot_0.1.10               deldir_0.2-9             
##  [61] locfit_1.5-9.4            tidyselect_1.1.0         
##  [63] labeling_0.4.2            rlang_0.4.10             
##  [65] reshape2_1.4.4            later_1.1.0.1            
##  [67] munsell_0.5.0             tools_4.0.3              
##  [69] generics_0.1.0            ggridges_0.5.3           
##  [71] evaluate_0.14             stringr_1.4.0            
##  [73] fastmap_1.0.1             yaml_2.2.1               
##  [75] goftest_1.2-2             knitr_1.30               
##  [77] bit64_4.0.5               fitdistrplus_1.1-3       
##  [79] purrr_0.3.4               RANN_2.6.1               
##  [81] sparseMatrixStats_1.2.0   pbapply_1.4-3            
##  [83] future_1.21.0             nlme_3.1-151             
##  [85] mime_0.9                  formatR_1.7              
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



