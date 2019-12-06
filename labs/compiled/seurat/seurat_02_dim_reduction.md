---
title: "Dimensionality reduction"
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

alldata <- readRDS("data/3pbmc_qc.rds")
```

### Feature selection

Next, we first need to define which features/genes are important in our dataset to distinguish cell types. For this purpose, we need to find genes that are highly variable across cells, which in turn will also provide a good separation of the cell clusters.


```r
suppressWarnings(suppressMessages(alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 850,verbose = FALSE,assay = "RNA")))
top20 <- head(VariableFeatures(alldata), 20)

LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)
```

```
## Warning: Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
## Please use `as_label()` or `as_name()` instead.
## This warning is displayed once per session.
```

```
## When using repel, set xnudge and ynudge to 0 for optimal results
```

```
## Warning: Transformation introduced infinite values in continuous x-axis
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

### Z-score transformation

Now that the data is prepared, we now proceed with PCA. Since each gene has a different expression level, it means that genes with higher expression values will naturally have higher variation that will be captured by PCA. This means that we need to somehow give each gene a similar weight when performing PCA (see below). The common practice is to center and scale each gene before performing PCA. This exact scaling is called Z-score normalization it is very useful for PCA, clustering and plotting heatmaps. <br>Additionally, we can use regression to remove any unwanted sources of variation from the dataset, such as `cell cycle`, `sequencing depth`, `percent mitocondria`. This is achieved by doing a generalized linear regression using these parameters as covariates in the model. Then the residuals of the model are taken as the "regressed data". Although perhaps not in the best way, batch effect regression can also be done here.


```r
alldata <- ScaleData(alldata, vars.to.regress = "percent_mito", assay = "RNA")
```

```
## Regressing out percent_mito
```

```
## Centering and scaling data matrix
```


## PCA
***

Performing PCA has many useful applications and interpretations, which much depends on the data used. In the case of life sciences, we want to segregate samples based on gene expression patterns in the data.


```r
alldata <- RunPCA(alldata, npcs = 50, reduction.name = "PCA_on_RNA", assay = "RNA",verbose = F)
```

We then plot the first principal components.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "PCA_on_RNA", group.by = "orig.ident",dims = 1:2),
  DimPlot(alldata, reduction = "PCA_on_RNA", group.by = "orig.ident",dims = 3:4),
  DimPlot(alldata, reduction = "PCA_on_RNA", group.by = "orig.ident",dims = 5:6) )
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

To identify which genes (Seurat) or metadata paramters (Scater/Scran) contribute the most to each PC, one can retreive the loading matrix information. Unfortunatelly this is not implemented in Scater/Scran, so you will need to compute PCA using `logcounts`.


```r
VizDimLoadings(alldata, dims = 1:5, reduction = "PCA_on_RNA",ncol = 5,balanced = T)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


We can also plot the amount of variance explained by each PC.


```r
ElbowPlot(alldata, reduction = "PCA_on_RNA",ndims = 50)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

Based on this plot, we can see that the top 8 PCs retain a lot of information, while other PCs contain pregressivelly less. However, it is still advisable to use more PCs since they might contain informaktion about rare cell types (such as platelets and DCs in this dataset)

## tSNE
***

We can now run [BH-tSNE](https://arxiv.org/abs/1301.3342).


```r
alldata <- RunTSNE(alldata, reduction = "PCA_on_RNA", dims = 1:30, reduction.name = "TSNE_on_RNA",
                   perplexity=30,
                   max_iter=1000,
                   theta=0.5,
                   eta=200,
                   num_threads=0 )
#see ?Rtsne and ?RunTSNE for more info
```

We can now plot the tSNE colored per dataset. We can clearly see the effect of batches present in the dataset.


```r
plot_grid(ncol = 3,DimPlot(alldata, reduction = "TSNE_on_RNA", group.by = "orig.ident"))
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


***
## UMAP
***

We can now run [UMAP](https://arxiv.org/abs/1802.03426) for cell embeddings.


```r
alldata <- RunUMAP(alldata, reduction = "PCA_on_RNA", dims = 1:30,reduction.name = "UMAP_on_RNA",
                   n.components=2,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )
#see ?RunUMAP for more info
```

Another usefullness of UMAP is that it is not limitted by the number of dimensions the data cen be reduced into (unlike tSNE). We can simply reduce the dimentions altering the `n.components` parameter.


```r
alldata <- RunUMAP(alldata, reduction.name = "UMAP10_on_RNA",
                   reduction = "PCA_on_RNA", 
                   dims = 1:30,
                   n.components=10,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap10_on_rna_'
```

```r
#see ?RunUMAP for more info
```

We can now plot the UMAP colored per dataset. Although less distinct as in the tSNE, we still see quite an effect of the different batches in the data.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "UMAP_on_RNA", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_RNA"),
  DimPlot(alldata, reduction = "UMAP10_on_RNA", group.by = "orig.ident",dims = 1:2)+ ggplot2::ggtitle(label ="UMAP10_on_RNA"),
  DimPlot(alldata, reduction = "UMAP10_on_RNA", group.by = "orig.ident",dims = 3:4)+ ggplot2::ggtitle(label ="UMAP10_on_RNA")
)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

We can now plot PCA, UMAP and tSNE side by side for comparison. Here, we can conclude that our dataset contains a batch effect that needs to be corrected before proceeding to clustering and differential gene expression analysis.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "PCA_on_RNA", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "TSNE_on_RNA", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "UMAP_on_RNA", group.by = "orig.ident")
)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

## Using ScaledData and graphs for DR
***

Althought running a sencond dimmensionality reduction (i.e tSNE or UMAP) on PCA would be a standard approach (because it allows higher computation efficiency), the options are actually limiteless. Below we will show a couple of other common options such as running directly on the scaled data (which was used for PCA) or on a graph built from scaled data. We will show from now on only UMAP, but the same applies for tSNE.

### Using ScaledData for UMAP

To run tSNE or UMAP on the scaled data, one firts needs to select the number of variables to use. This is because including dimentions that do contribute to the separation of your cell types will in the end mask those differences. Another reason for it is because running with all genes/features also will take longer or might be computationally unfeasible. Therefore we will use the scaled data of the highly variable genes.


```r
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_ScaleData",
                   features = alldata@assays$RNA@var.features,
                   assay = "RNA",
                   n.components=2,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap_on_scaledata_'
```

### Using a Graph for UMAP

To run tSNE or UMAP on the a graph, we first need to build a graph from the data. In fact, both tSNE and UMAP first build a graph from the data using a specified distance metrix and then optimize the embedding. Since a graph is just a matrix containing distances from cell to cell and as such, you can run either UMAP or tSNE using any other distance metric desired. Euclidean and Correlation are ususally the most commonly used.


```r
#Build Graph
alldata <- FindNeighbors(alldata,
                         reduction = "PCA_on_RNA",
                         graph.name = "SNN",
                         assay = "RNA",
                         k.param = 20,
                         features = alldata@assays$RNA@var.features)
```

```
## Computing nearest neighbor graph
```

```
## Computing SNN
```

```r
#Run UMAP on a graph
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph",
                   graph = "SNN",
                   assay = "RNA" )
```

```
## Warning: Cannot add objects with duplicate keys (offending key: UMAP_), setting
## key to 'umap_on_graph_'
```

We can now plot the UMAP comparing both on PCA vs ScaledSata vs Graph.


```r
plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "UMAP_on_RNA", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_RNA"),
  DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_ScaleData"),
  DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_Graph")
)
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
myfeatures <- c("CD3E","CD4","CD8A","NKG7","GNLY","MS4A1","CD14","LYZ","MS4A7","FCGR3A","CST3","FCER1A")
FeaturePlot(alldata, reduction = "UMAP_on_RNA",dims = 1:2,
            features = myfeatures,ncol = 3,order = T)
```

![](seurat_02_dim_reduction_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


We can finally save the object for use in future steps.


```r
saveRDS(alldata,"data/3pbmc_qc_dr.rds")
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
##  [70] stringi_1.4.3            caTools_1.17.1.2         bibtex_0.4.2            
##  [73] Rdpack_0.11-0            SDMTools_1.1-221.1       rlang_0.4.2             
##  [76] pkgconfig_2.0.3          bitops_1.0-6             evaluate_0.14           
##  [79] lattice_0.20-38          Rhdf5lib_1.4.3           ROCR_1.0-7              
##  [82] purrr_0.3.3              htmlwidgets_1.5.1        labeling_0.3            
##  [85] bit_1.1-14               tidyselect_0.2.5         plyr_1.8.4              
##  [88] magrittr_1.5             R6_2.4.1                 gplots_3.0.1.1          
##  [91] pillar_1.4.2             withr_2.1.2              fitdistrplus_1.0-14     
##  [94] survival_2.44-1.1        RCurl_1.95-4.12          tibble_2.1.3            
##  [97] future.apply_1.3.0       tsne_0.1-3               crayon_1.3.4            
## [100] hdf5r_1.2.0              KernSmooth_2.23-15       plotly_4.9.1            
## [103] rmarkdown_1.17           viridis_0.5.1            locfit_1.5-9.1          
## [106] grid_3.5.1               data.table_1.11.6        metap_1.1               
## [109] digest_0.6.23            tidyr_1.0.0              R.utils_2.9.0           
## [112] munsell_0.5.0            beeswarm_0.2.3           viridisLite_0.3.0       
## [115] vipor_0.4.5
```



