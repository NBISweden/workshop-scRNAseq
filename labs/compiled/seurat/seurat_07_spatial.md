---
title: "Seurat: Spatial Transcriptomics"
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'Januari 27, 2023'
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

# Spatial transcriptomics
***

This tutorial is adapted from the Seurat vignette: https://satijalab.org/seurat/v3.2/spatial_vignette.html

Spatial transcriptomic data with the Visium platform is in many ways similar to scRNAseq data. It contains UMI counts for 5-20 cells instead of single cells, but is still quite sparse in the same way as scRNAseq data is, but with the additional information about spatial location in the tissue. 

Here we will first run quality control in a similar manner to scRNAseq data, then filtering, dimensionality reduction, integration and clustering. Then we will use scRNAseq data from mouse cortex to run LabelTransfer to predict celltypes in the Visium spots. 

We will use two **Visium** spatial transcriptomics dataset of the mouse brain (Sagittal), which are publicly available from the [10x genomics website](https://support.10xgenomics.com/spatial-gene-expression/datasets/). Note, that these dataset have already been filtered for spots that does not overlap with the tissue.

### Load packages


```r
devtools::install_github("satijalab/seurat-data")

suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(SeuratData))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(dplyr))
```


### Load ST data

The package `SeuratData` has some seurat objects for different datasets. Among those are spatial transcriptomics data from mouse brain and kidney. Here we will download and process sections from the mouse brain. 


```r
outdir = "data/spatial/"
dir.create(outdir, showWarnings = F)

# to list available datasets in SeuratData you can run AvailableData()

# first we dowload the dataset
if (!("stxBrain.SeuratData" %in% rownames(InstalledData()))) {
    InstallData("stxBrain")
}


# now we can list what datasets we have downloaded
InstalledData()
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Dataset"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Version"],"name":[2],"type":["pckg_vrs"],"align":["right"]},{"label":["Summary"],"name":[3],"type":["chr"],"align":["left"]},{"label":["species"],"name":[4],"type":["chr"],"align":["left"]},{"label":["system"],"name":[5],"type":["chr"],"align":["left"]},{"label":["ncells"],"name":[6],"type":["dbl"],"align":["right"]},{"label":["tech"],"name":[7],"type":["chr"],"align":["left"]},{"label":["seurat"],"name":[8],"type":["chr"],"align":["left"]},{"label":["default.dataset"],"name":[9],"type":["chr"],"align":["left"]},{"label":["disk.datasets"],"name":[10],"type":["chr"],"align":["left"]},{"label":["other.datasets"],"name":[11],"type":["chr"],"align":["left"]},{"label":["notes"],"name":[12],"type":["chr"],"align":["left"]},{"label":["Installed"],"name":[13],"type":["lgl"],"align":["right"]},{"label":["InstalledVersion"],"name":[14],"type":["pckg_vrs"],"align":["right"]}],"data":[{"1":"stxBrain","2":"<pckg_vrs>","3":"10X Genomics Visium Mouse Brain Dataset","4":"mouse","5":"brain","6":"12167","7":"visium","8":"NA","9":"__NA__","10":"NA","11":"posterior1, posterior2, anterior1, anterior2","12":"One sample split across four datasets as paired anterior/posterior slices","13":"TRUE","14":"<pckg_vrs>","_rn_":"stxBrain.SeuratData"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# now we will load the seurat object for one section
brain1 <- LoadData("stxBrain", type = "anterior1")
brain2 <- LoadData("stxBrain", type = "posterior1")
```

Merge into one seurat object


```r
brain <- merge(brain1, brain2)

brain
```

```
## An object of class Seurat 
## 31053 features across 6049 samples within 1 assay 
## Active assay: Spatial (31053 features, 0 variable features)
##  2 images present: anterior1, posterior1
```

As you can see, now we do not have the assay "RNA", but instead an assay called "Spatial". 


# Quality control
***

Similar to scRNAseq we use statistics on number of counts, number of features and percent mitochondria for quality control. 

Now the counts and feature counts are calculated on the Spatial assay, so they are named  "nCount_Spatial" and "nFeature_Spatial".


```r
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "percent_mito")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")


VlnPlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito",
    "percent_hb"), pt.size = 0.1, ncol = 2) + NoLegend()
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

We can also plot the same data onto the tissue section.


```r
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito",
    "percent_hb"))
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


As you can see, the spots with low number of counts/features and high mitochondrial content are mainly towards the edges of the tissue. It is quite likely that these regions are damaged tissue. You may also see regions within a tissue with low quality if you have tears or folds in your section. 

But remember, for some tissue types, the amount of genes expressed and proportion mitochondria may also be a biological features, so bear in mind what tissue you are working on and what these features mean.

### Filter

Select all spots with less than **25%** mitochondrial reads, less than **20%** hb-reads and at least **500** detected genes. You must judge for yourself based on your knowledge of the tissue what are appropriate filtering criteria for your dataset.



```r
brain <- brain[, brain$nFeature_Spatial > 500 & brain$percent_mito < 25 & brain$percent_hb <
    20]
```

And replot onto tissue section:


```r
SpatialFeaturePlot(brain, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

### Top expressed genes
As for scRNAseq data, we will look at what the top expressed genes are.


```r
C <- brain@assays$Spatial@counts
C@x <- C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

As you can see, the mitochondrial genes are among the top expressed genes. Also the lncRNA gene *Bc1* (brain cytoplasmic RNA 1). Also one hemoglobin gene.

### Filter genes
We will remove *Bc1* gene, hemoglobin genes (blood contamination) and the mitochondrial genes.


```r
dim(brain)
```

```
## [1] 31053  5789
```

```r
# Filter Bl1
brain <- brain[!grepl("Bc1", rownames(brain)), ]

# Filter Mitocondrial
brain <- brain[!grepl("^mt-", rownames(brain)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
brain <- brain[!grepl("^Hb.*-", rownames(brain)), ]

dim(brain)
```

```
## [1] 31031  5789
```

# Analysis
***

For ST data, the Seurat team recommends to use SCTranform for normalization, so we will do that. `SCTransform` will select variable genes and normalize in one step.


```r
brain <- SCTransform(brain, assay = "Spatial", verbose = TRUE)  #, method = 'poisson')
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |====                                                                  |   5%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |=======                                                               |  11%
  |                                                                            
  |=========                                                             |  13%
  |                                                                            
  |===========                                                           |  16%
  |                                                                            
  |=============                                                         |  18%
  |                                                                            
  |===============                                                       |  21%
  |                                                                            
  |=================                                                     |  24%
  |                                                                            
  |==================                                                    |  26%
  |                                                                            
  |====================                                                  |  29%
  |                                                                            
  |======================                                                |  32%
  |                                                                            
  |========================                                              |  34%
  |                                                                            
  |==========================                                            |  37%
  |                                                                            
  |============================                                          |  39%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |===============================                                       |  45%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=======================================                               |  55%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |==========================================                            |  61%
  |                                                                            
  |============================================                          |  63%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |=====================================================                 |  76%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |=========================================================             |  82%
  |                                                                            
  |===========================================================           |  84%
  |                                                                            
  |=============================================================         |  87%
  |                                                                            
  |===============================================================       |  89%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |==================================================================    |  95%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |====                                                                  |   5%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |=======                                                               |  11%
  |                                                                            
  |=========                                                             |  13%
  |                                                                            
  |===========                                                           |  16%
  |                                                                            
  |=============                                                         |  18%
  |                                                                            
  |===============                                                       |  21%
  |                                                                            
  |=================                                                     |  24%
  |                                                                            
  |==================                                                    |  26%
  |                                                                            
  |====================                                                  |  29%
  |                                                                            
  |======================                                                |  32%
  |                                                                            
  |========================                                              |  34%
  |                                                                            
  |==========================                                            |  37%
  |                                                                            
  |============================                                          |  39%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |===============================                                       |  45%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=======================================                               |  55%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |==========================================                            |  61%
  |                                                                            
  |============================================                          |  63%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |=====================================================                 |  76%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |=========================================================             |  82%
  |                                                                            
  |===========================================================           |  84%
  |                                                                            
  |=============================================================         |  87%
  |                                                                            
  |===============================================================       |  89%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |==================================================================    |  95%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |======================================================================| 100%
```


Now we can plot gene expression of individual genes, the gene *Hpca* is a strong hippocampal marker and *Ttr* is a marker of the choroid plexus.


```r
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

If you want to see the tissue better, you can modify point size and transparency of the points.


```r
SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1, alpha = c(0.1, 1))
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


### Dimensionality reduction and clustering
We can then now run dimensionality reduction and clustering using the same workflow as we use for scRNAseq analysis. 

But make sure you run it on the `SCT` assay.


```r
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
```

We can then plot clusters onto umap or onto the tissue section.


```r
DimPlot(brain, reduction = "umap", group.by = c("ident", "orig.ident"))
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
SpatialDimPlot(brain)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

We can also plot each cluster separately


```r
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(brain), facet.highlight = TRUE,
    ncol = 5)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

### Integration

Quite often there are strong batch effects between different ST sections, so it may be a good idea to integrate the data across sections.

We will do a similar integration as in the Data Integration lab, but this time we will use the SCT assay for integration. Therefore we need to run `PrepSCTIntegration` which will compute the sctransform residuals for all genes in both the datasets. 


```r
# create a list of the original data that we loaded to start with
st.list = list(anterior1 = brain1, posterior1 = brain2)

# run SCT on both datasets
st.list = lapply(st.list, SCTransform, assay = "Spatial", method = "poisson")
```

```
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================                                                   |  28%
  |                                                                            
  |=====================                                                 |  31%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |========================================================              |  81%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================                                                   |  28%
  |                                                                            
  |=====================                                                 |  31%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |========================================================              |  81%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================                                                   |  28%
  |                                                                            
  |=====================                                                 |  31%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |========================================================              |  81%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |======================================================================| 100%
## 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |==============                                                        |  19%
  |                                                                            
  |================                                                      |  22%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |===================                                                   |  28%
  |                                                                            
  |=====================                                                 |  31%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |=================================================                     |  69%
  |                                                                            
  |===================================================                   |  72%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |======================================================                |  78%
  |                                                                            
  |========================================================              |  81%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |======================================================================| 100%
```

```r
# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB


st.features = SelectIntegrationFeatures(st.list, nfeatures = 3000, verbose = FALSE)
st.list <- PrepSCTIntegration(object.list = st.list, anchor.features = st.features,
    verbose = FALSE)
```

Now we can perform the actual integration.


```r
int.anchors <- FindIntegrationAnchors(object.list = st.list, normalization.method = "SCT",
    verbose = FALSE, anchor.features = st.features)
brain.integrated <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
    verbose = FALSE)

rm(int.anchors, st.list)
gc()
```

```
##             used   (Mb) gc trigger   (Mb)   max used   (Mb)
## Ncells   3637071  194.3    6301768  336.6    6301768  336.6
## Vcells 590889428 4508.2 1192885632 9101.0 1192237793 9096.1
```

Then we run dimensionality reduction and clustering as before.


```r
brain.integrated <- RunPCA(brain.integrated, verbose = FALSE)
brain.integrated <- FindNeighbors(brain.integrated, dims = 1:30)
brain.integrated <- FindClusters(brain.integrated, verbose = FALSE)
brain.integrated <- RunUMAP(brain.integrated, dims = 1:30)
```


```r
DimPlot(brain.integrated, reduction = "umap", group.by = c("ident", "orig.ident"))
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

```r
SpatialDimPlot(brain.integrated)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

Do you see any differences between the integrated and non-integrated clustering? Judge for yourself, which of the clusterings do you think looks best? 
As a reference, you can compare to brain regions in the [Allen brain atlas](https://mouse.brain-map.org/experiment/thumbnails/100042147?image_type=atlas). 

### Identification of Spatially Variable Features

There are two main workflows to identify molecular features that correlate with spatial location within a tissue. The first is to perform differential expression based on spatially distinct clusters, the other is to find features that have spatial patterning without taking clusters or spatial annotation into account. 

First, we will do differential expression between clusters just as we did for the scRNAseq data before.

```r
# differential expression between cluster 1 and cluster 6
de_markers <- FindMarkers(brain.integrated, ident.1 = 5, ident.2 = 6)

# plot top markers
SpatialFeaturePlot(object = brain.integrated, features = rownames(de_markers)[1:3],
    alpha = c(0.1, 1), ncol = 3)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-18-1.png)<!-- -->



In `FindSpatiallyVariables` the default method in Seurat (method = 'markvariogram), is inspired by the Trendsceek, which models spatial transcriptomics data as a mark point process and computes a 'variogram', which identifies genes whose expression level is dependent on their spatial location. More specifically, this process calculates gamma(r) values measuring the dependence between two spots a certain "r" distance apart. By default, we use an r-value of '5' in these analyses, and only compute these values for variable genes (where variation is calculated independently of spatial location) to save time.


**OBS!** Takes a long time to run, so skip this step for now!


```r
# brain <- FindSpatiallyVariableFeatures(brain, assay = 'SCT', features =
# VariableFeatures(brain)[1:1000], selection.method = 'markvariogram')

# We would get top features from SpatiallyVariableFeatures top.features <-
# head(SpatiallyVariableFeatures(brain, selection.method = 'markvariogram'), 6)
```




# Single cell data

We can use a scRNAseq dataset as a reference to predict the proportion of different celltypes in the Visium spots. 

Keep in mind that it is important to have a reference that contains all the celltypes you expect to find in your spots. Ideally it should be a scRNAseq reference from the exact same tissue. 

We will use a reference scRNAseq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol.


First download the seurat data from: https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1 to folder `data/spatial/` with command:



```bash

FILE="./data/spatial/allen_cortex.rds"

if [ -e $FILE ]
then
   echo "File $FILE is downloaded."
else
   echo "Downloading $FILE"
   mkdir -p data/spatial
   wget  -O data/spatial/allen_cortex.rds https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
fi

```

```
## File ./data/spatial/allen_cortex.rds is downloaded.
```

For speed, and for a more fair comparison of the celltypes, we will subsample all celltypes to a maximum of **200** cells per class (`subclass`).


```r
allen_reference <- readRDS("data/spatial/allen_cortex.rds")

# check number of cells per subclass
table(allen_reference$subclass)
```

```
## 
##      Astro         CR       Endo    L2/3 IT         L4      L5 IT      L5 PT 
##        368          7         94        982       1401        880        544 
##      L6 CT      L6 IT        L6b      Lamp5 Macrophage      Meis2         NP 
##        960       1872        358       1122         51         45        362 
##      Oligo       Peri      Pvalb   Serpinf1        SMC       Sncg        Sst 
##         91         32       1337         27         55        125       1741 
##        Vip       VLMC 
##       1728         67
```

```r
# select 200 cells per subclass, fist set subclass ass active.ident
Idents(allen_reference) <- allen_reference$subclass
allen_reference <- subset(allen_reference, cells = WhichCells(allen_reference, downsample = 200))

# check again number of cells per subclass
table(allen_reference$subclass)
```

```
## 
##      Astro         CR       Endo    L2/3 IT         L4      L5 IT      L5 PT 
##        200          7         94        200        200        200        200 
##      L6 CT      L6 IT        L6b      Lamp5 Macrophage      Meis2         NP 
##        200        200        200        200         51         45        200 
##      Oligo       Peri      Pvalb   Serpinf1        SMC       Sncg        Sst 
##         91         32        200         27         55        125        200 
##        Vip       VLMC 
##        200         67
```

Then run normalization and dimensionality reduction.


```r
# First run SCTransform and PCA
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE, method = "poisson") %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, label = TRUE)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

# Subset ST for cortex
Since the scRNAseq dataset was generated from the mouse cortex, we will subset the visium dataset in order to select mainly the spots part of the cortex. Note that the integration can also be performed on the whole brain slice, but it would give rise to false positive cell type assignments  and therefore it should be interpreted with more care.


```r
# subset for the anterior dataset
cortex <- subset(brain.integrated, orig.ident == "anterior1")

# there seems to be an error in the subsetting, so the posterior1 image is not
# removed, do it manually
cortex@images$posterior1 = NULL

# subset for a specific region
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

# also subset for FC (Frontal cortex) clusters
cortex <- subset(cortex, idents = c(1, 2, 3, 4, 5))

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-22-1.png)<!-- -->



# Deconvolution
Deconvolution is a method to estimate the abundance (or proportion) of different celltypes in a bulkRNAseq dataset using a single cell reference. As the Visium data can be seen as a small bulk, we can both use methods for traditional bulkRNAseq as well as methods especially developed for Visium data. Some methods for deconvolution are DWLS, cell2location, Tangram, Stereoscope, RCTD, SCDC and many more. 

Here we will use SCDC for deconvolution of celltypes in the Visium spots. For more information on the tool please check their website: https://meichendong.github.io/SCDC/articles/SCDC.html
First, make sure the packages you need are installed. All dependencies should be in the conda environment.



```r
inst = installed.packages()

if (!("xbioc" %in% rownames(inst))) {
    remotes::install_github("renozao/xbioc", dependencies = FALSE)
}
if (!("SCDC" %in% rownames(inst))) {
    remotes::install_github("meichendong/SCDC", dependencies = FALSE)
}

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))
```


### Select genes for deconvolution
Most deconvolution methods does a prior gene selection and there are different options that are used:

* Use variable genes in the SC data.
* Use variable genes in both SC and ST data
* DE genes between clusters in the SC data.

In this case we will use top DEGs per cluster, so first we have to perform differential expression on the scRNAseq data.

For SCDC we want to find unique markers per cluster, so we select top 20 DEGs per cluster.
Ideally you should run with a larger set of genes, perhaps 100 genes per cluster to get better results. However, for the sake of speed, we are now selecting only top20 genes and it still takes about 10 minutes to run.
word

```r
allen_reference@active.assay = "RNA"

markers_sc <- FindAllMarkers(allen_reference, only.pos = TRUE, logfc.threshold = 0.1,
    test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
    return.thresh = 0.05, assay = "RNA")

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(cortex), ]


# Select top 20 genes per cluster, select top by first p-value, then absolute
# diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
    1))
markers_sc %>%
    group_by(cluster) %>%
    top_n(-100, p_val) %>%
    top_n(50, pct.diff) %>%
    top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))
```

### Create ExpressionSets

For SCDC both the SC and the ST data need to be in the format of an **ExpressionSet** with the count matrices as `AssayData`. We also subset the matrices for the genes we selected in the previous step.


```r
eset_SC <- ExpressionSet(assayData = as.matrix(allen_reference@assays$RNA@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(allen_reference@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(cortex@assays$Spatial@counts[m_feats,
    ]), phenoData = AnnotatedDataFrame(cortex@meta.data))
```

### Deconvolve
We then run the deconvolution defining the celltype of interest as "subclass" column in the single cell data.


```r
deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "subclass",
    ct.sub = as.character(unique(eset_SC$subclass)))
```

Now we have a matrix with predicted proportions of each celltypes for each visium spot in `prop.est.mvw`:


```r
head(deconvolution_crc$prop.est.mvw)
```

```
##                      Lamp5 Sncg Serpinf1 Vip          Sst      Pvalb       Endo
## AAACAGAGCGACTCCT-1_1     0    0        0   0 0.0019454427 0.00000000 0.00000000
## AAACCGGGTAGGTACC-1_1     0    0        0   0 0.0915908853 0.00000000 0.01228367
## AAAGGGATGTAGCAAG-1_1     0    0        0   0 0.2744276526 0.00000000 0.01113891
## AAATAACCATACGGGA-1_1     0    0        0   0 0.0002310634 0.00000000 0.00000000
## AAATCGTGTACCACAA-1_1     0    0        0   0 0.0696595459 0.00000000 0.02797452
## AAATGATTCGATCAGC-1_1     0    0        0   0 0.0190284766 0.04230543 0.01467511
##                              Peri        L6 CT          L6b        L6 IT CR
## AAACAGAGCGACTCCT-1_1 1.882987e-05 0.0000000000 7.206678e-05 0.000000e+00  0
## AAACCGGGTAGGTACC-1_1 0.000000e+00 0.0000000000 0.000000e+00 0.000000e+00  0
## AAAGGGATGTAGCAAG-1_1 0.000000e+00 0.0000000000 0.000000e+00 2.234847e-04  0
## AAATAACCATACGGGA-1_1 0.000000e+00 0.0000000000 0.000000e+00 1.722518e-04  0
## AAATCGTGTACCACAA-1_1 5.210378e-06 0.0002166623 2.474756e-01 3.174855e-05  0
## AAATGATTCGATCAGC-1_1 2.742467e-02 0.2912369043 0.000000e+00 1.400229e-01  0
##                          L2/3 IT      L5 PT NP        L4        Oligo
## AAACAGAGCGACTCCT-1_1 0.777266726 0.00000000  0 0.0000000 0.0000000000
## AAACCGGGTAGGTACC-1_1 0.082075888 0.00000000  0 0.4551348 0.0000000000
## AAAGGGATGTAGCAAG-1_1 0.000000000 0.00000000  0 0.1956948 0.0000000000
## AAATAACCATACGGGA-1_1 0.945629429 0.00000000  0 0.0000000 0.0000000000
## AAATCGTGTACCACAA-1_1 0.305548439 0.00000000  0 0.0000000 0.2078974324
## AAATGATTCGATCAGC-1_1 0.009593289 0.09912164  0 0.1596236 0.0001027801
##                             L5 IT        Meis2      Astro  Macrophage VLMC SMC
## AAACAGAGCGACTCCT-1_1 0.0002334628 0.000000e+00 0.18895054 0.031512935    0   0
## AAACCGGGTAGGTACC-1_1 0.0000000000 2.639359e-06 0.33870721 0.020204935    0   0
## AAAGGGATGTAGCAAG-1_1 0.3535093520 0.000000e+00 0.15672588 0.008279916    0   0
## AAATAACCATACGGGA-1_1 0.0000000000 0.000000e+00 0.01209326 0.041873991    0   0
## AAATCGTGTACCACAA-1_1 0.0000000000 0.000000e+00 0.11579684 0.025394016    0   0
## AAATGATTCGATCAGC-1_1 0.1003021059 0.000000e+00 0.09106008 0.005503047    0   0
```

Now we take the deconvolution output and add it to the Seurat object as a new assay.


```r
cortex@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))

# Seems to be a bug in SeuratData package that the key is not set and any
# plotting function etc. will throw an error.
if (length(cortex@assays$SCDC@key) == 0) {
    cortex@assays$SCDC@key = "scdc_"
}
```



```r
DefaultAssay(cortex) <- "SCDC"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2,
    crop = TRUE)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-29-1.png)<!-- -->


Based on these prediction scores, we can also predict cell types whose location is spatially restricted. We use the same methods based on marked point processes to define spatially variable features, but use the cell type prediction scores as the "marks" rather than gene expression.


```r
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "SCDC", selection.method = "markvariogram",
    features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

#ST_R13.7:


```r
VlnPlot(cortex, group.by = "seurat_clusters", features = top.clusters, pt.size = 0,
    ncol = 2)
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

Keep in mind that the deconvolution results are just predictions, depending on how well your scRNAseq data covers the celltypes that are present in the ST data and on how parameters, gene selection etc. are tuned you may get different results.

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn:**

Subset for another region that does not contain cortex cells and check what you get from the label transfer. 

Suggested region is the right end of the posterial section that you can select like this:

</div>



```r
# subset for the anterior dataset
subregion <- subset(brain.integrated, orig.ident == "posterior1")

# there seems to be an error in the subsetting, so the posterior1 image is not
# removed, do it manually
subregion@images$anterior1 = NULL

# subset for a specific region
subregion <- subset(subregion, posterior1_imagecol > 400, invert = FALSE)

p1 <- SpatialDimPlot(subregion, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(subregion, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```

![](seurat_07_spatial_files/figure-html/unnamed-chunk-32-1.png)<!-- -->


### Session info


```r
sessionInfo()
```

```
## R version 4.2.2 (2022-10-31)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: macOS Big Sur ... 10.16
## 
## Matrix products: default
## BLAS/LAPACK: /Users/nimra236/opt/anaconda3/envs/scRNAseq2023/lib/libopenblasp-r0.3.21.dylib
## 
## locale:
## [1] sv_SE.UTF-8/sv_SE.UTF-8/sv_SE.UTF-8/C/sv_SE.UTF-8/sv_SE.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] Biobase_2.58.0            BiocGenerics_0.44.0      
##  [3] SCDC_0.0.0.9000           patchwork_1.1.2          
##  [5] ggplot2_3.4.0             SeuratObject_4.1.3       
##  [7] Seurat_4.3.0              stxBrain.SeuratData_0.1.1
##  [9] SeuratData_0.2.2          dplyr_1.0.10             
## [11] Matrix_1.5-3             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2             spatstat.explore_3.0-5 reticulate_1.26       
##   [4] tidyselect_1.2.0       RSQLite_2.2.20         AnnotationDbi_1.60.0  
##   [7] htmlwidgets_1.6.1      grid_4.2.2             Rtsne_0.16            
##  [10] devtools_2.4.5         munsell_0.5.0          codetools_0.2-18      
##  [13] ica_1.0-3              future_1.30.0          miniUI_0.1.1.1        
##  [16] withr_2.5.0            spatstat.random_3.1-3  colorspace_2.1-0      
##  [19] progressr_0.13.0       highr_0.10             knitr_1.42            
##  [22] rstudioapi_0.14        stats4_4.2.2           ROCR_1.0-11           
##  [25] tensor_1.5             listenv_0.9.0          labeling_0.4.2        
##  [28] GenomeInfoDbData_1.2.9 polyclip_1.10-4        bit64_4.0.5           
##  [31] farver_2.1.1           pheatmap_1.0.12        parallelly_1.34.0     
##  [34] vctrs_0.5.2            generics_0.1.3         xfun_0.36             
##  [37] R6_2.5.1               GenomeInfoDb_1.34.1    ggbeeswarm_0.7.1      
##  [40] L1pack_0.41-2          bitops_1.0-7           spatstat.utils_3.0-1  
##  [43] cachem_1.0.6           assertthat_0.2.1       promises_1.2.0.1      
##  [46] scales_1.2.1           beeswarm_0.4.0         gtable_0.3.1          
##  [49] globals_0.16.2         processx_3.8.0         goftest_1.2-3         
##  [52] rlang_1.0.6            splines_4.2.2          lazyeval_0.2.2        
##  [55] checkmate_2.1.0        spatstat.geom_3.0-5    BiocManager_1.30.19   
##  [58] yaml_2.3.7             reshape2_1.4.4         abind_1.4-5           
##  [61] backports_1.4.1        httpuv_1.6.8           tools_4.2.2           
##  [64] usethis_2.1.6          ellipsis_0.3.2         jquerylib_0.1.4       
##  [67] RColorBrewer_1.1-3     sessioninfo_1.2.2      ggridges_0.5.4        
##  [70] Rcpp_1.0.10            plyr_1.8.8             zlibbioc_1.44.0       
##  [73] purrr_1.0.1            RCurl_1.98-1.9         ps_1.7.2              
##  [76] prettyunits_1.1.1      deldir_1.0-6           pbapply_1.7-0         
##  [79] cowplot_1.1.1          urlchecker_1.0.1       S4Vectors_0.36.0      
##  [82] zoo_1.8-11             ggrepel_0.9.2          cluster_2.1.4         
##  [85] fs_1.6.0               magrittr_2.0.3         data.table_1.14.6     
##  [88] scattermore_0.8        lmtest_0.9-40          RANN_2.6.1            
##  [91] fitdistrplus_1.1-8     matrixStats_0.63.0     pkgload_1.3.2         
##  [94] mime_0.12              evaluate_0.20          xtable_1.8-4          
##  [97] IRanges_2.32.0         gridExtra_2.3          compiler_4.2.2        
## [100] tibble_3.1.8           KernSmooth_2.23-20     crayon_1.5.2          
## [103] htmltools_0.5.4        later_1.3.0            tidyr_1.3.0           
## [106] DBI_1.1.3              formatR_1.14           MASS_7.3-58.2         
## [109] rappdirs_0.3.3         cli_3.6.0              parallel_4.2.2        
## [112] igraph_1.3.5           pkgconfig_2.0.3        registry_0.5-1        
## [115] sp_1.6-0               plotly_4.10.1          xbioc_0.1.19          
## [118] spatstat.sparse_3.0-0  vipor_0.4.5            bslib_0.4.2           
## [121] pkgmaker_0.32.7        XVector_0.38.0         fastmatrix_0.4-1245   
## [124] stringr_1.5.0          callr_3.7.3            digest_0.6.31         
## [127] sctransform_0.3.5      RcppAnnoy_0.0.20       spatstat.data_3.0-0   
## [130] Biostrings_2.66.0      rmarkdown_2.20         leiden_0.4.3          
## [133] uwot_0.1.14            curl_4.3.3             shiny_1.7.4           
## [136] lifecycle_1.0.3        nlme_3.1-161           jsonlite_1.8.4        
## [139] viridisLite_0.4.1      limma_3.54.0           fansi_1.0.4           
## [142] pillar_1.8.1           lattice_0.20-45        ggrastr_1.0.1         
## [145] KEGGREST_1.38.0        fastmap_1.1.0          httr_1.4.4            
## [148] pkgbuild_1.4.0         survival_3.5-0         glue_1.6.2            
## [151] remotes_2.4.2          png_0.1-8              bit_4.0.5             
## [154] stringi_1.7.12         sass_0.4.5             profvis_0.3.7         
## [157] nnls_1.4               blob_1.2.3             memoise_2.0.1         
## [160] irlba_2.3.5.1          future.apply_1.10.0
```
