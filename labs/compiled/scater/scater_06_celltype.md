---
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 23, 2021'
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

## Celltype prediction
***

 Celltype prediction can either be performed on indiviudal cells where each cell gets a predicted celltype label, or on the level of clusters. All methods are based on similarity to other datasets, single cell or sorted bulk RNAseq, or uses know marker genes for each celltype.

We will select one sample from the Covid data, `ctrl_13` and predict celltype by cell on that sample.

Some methods will predict a celltype to each cell based on what it is most similar to even if the celltype of that cell is not included in the reference. Other methods include an uncertainty so that cells with low similarity scores will be unclassified.
There are multiple different methods to predict celltypes, here we will just cover a few of those. 

Here we will use a reference PBMC dataset from the `scPred` package which is provided as a Seurat object with counts. And we will test classification based on the `scPred` and `scMap` methods. Finally we will use gene set enrichment predict celltype based on the DEGs of each cluster.

# Load and process data
## Covid-19 data
First, lets load required libraries and the saved object from the clustering step. Subset for one patient.


```r
suppressPackageStartupMessages({
    library(scater)
    library(scran)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(scPred)
    library(scmap)
})
```




```r
# load the data and select 'ctrl_13` sample
alldata <- readRDS("data/results/covid_qc_dr_int_cl.rds")
ctrl.sce <- alldata[, alldata@colData$sample == "ctrl.13"]

# remove all old dimensionality reductions as they will mess up the analysis
# further down
reducedDims(ctrl.sce) <- NULL
```

## Reference data
Then, load the reference dataset with annotated labels. Also, run all steps of the normal analysis pipeline with normalizaiton, variable gene selection, scaling and dimensionality reduction.



```r
reference <- scPred::pbmc_1

reference
```

```
## An object of class Seurat 
## 32838 features across 3500 samples within 1 assay 
## Active assay: RNA (32838 features, 0 variable features)
```

Convert to a SCE object.


```r
ref.sce = Seurat::as.SingleCellExperiment(reference)
```


## Rerun analysis pipeline
Run normalization, feature selection and dimensionality reduction


```r
# Normalize
ref.sce <- computeSumFactors(ref.sce)
ref.sce <- logNormCounts(ref.sce)

# Variable genes
var.out <- modelGeneVar(ref.sce, method = "loess")
hvg.ref <- getTopHVGs(var.out, n = 1000)

# Dim reduction
ref.sce <- runPCA(ref.sce, exprs_values = "logcounts", scale = T, ncomponents = 30, 
    subset_row = hvg.ref)
ref.sce <- runUMAP(ref.sce, dimred = "PCA")
```


```r
plotReducedDim(ref.sce, dimred = "UMAP", colour_by = "cell_type")
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


Run all steps of the analysis for the ctrl sample as well. Use the clustering from the integration lab with resolution 0.3.


```r
# Normalize
ctrl.sce <- computeSumFactors(ctrl.sce)
ctrl.sce <- logNormCounts(ctrl.sce)

# Variable genes
var.out <- modelGeneVar(ctrl.sce, method = "loess")
hvg.ctrl <- getTopHVGs(var.out, n = 1000)

# Dim reduction
ctrl.sce <- runPCA(ctrl.sce, exprs_values = "logcounts", scale = T, ncomponents = 30, 
    subset_row = hvg.ctrl)
ctrl.sce <- runUMAP(ctrl.sce, dimred = "PCA")
```


```r
plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "louvain_SNNk15")
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

# scMap
The scMap package is one method for projecting cells from a scRNA-seq experiment on to the cell-types or individual cells identified in a different experiment. It can be run on different levels, either projecting by cluster or by single cell, here we will try out both.

## scMap cluster
For scmap cell type labels must be stored in the `cell_type1` column of the `colData` slots, and gene ids that are consistent across both datasets must be stored in the `feature_symbol` column of the `rowData` slots.
Then we can select variable features in both datasets.



```r
# add in slot cell_type1
ref.sce@colData$cell_type1 = ref.sce@colData$cell_type
# create a rowData slot with feature_symbol
rd = data.frame(feature_symbol = rownames(ref.sce))
rownames(rd) = rownames(ref.sce)
rowData(ref.sce) = rd




# same for the ctrl dataset create a rowData slot with feature_symbol
rd = data.frame(feature_symbol = rownames(ctrl.sce))
rownames(rd) = rownames(ctrl.sce)
rowData(ctrl.sce) = rd


# select features
counts(ctrl.sce) <- as.matrix(counts(ctrl.sce))
logcounts(ctrl.sce) <- as.matrix(logcounts(ctrl.sce))
ctrl.sce <- selectFeatures(ctrl.sce, suppress_plot = TRUE)

counts(ref.sce) <- as.matrix(counts(ref.sce))
logcounts(ref.sce) <- as.matrix(logcounts(ref.sce))
ref.sce <- selectFeatures(ref.sce, suppress_plot = TRUE)
```

Then we need to index the reference dataset by cluster, default is the clusters in `cell_type1`.


```r
ref.sce <- indexCluster(ref.sce)
```

Now we project the Covid-19 dataset onto that index.


```r
project_cluster <- scmapCluster(projection = ctrl.sce, index_list = list(ref = metadata(ref.sce)$scmap_cluster_index))

# projected labels
table(project_cluster$scmap_cluster_labs)
```

```
## 
##      B cell  CD4 T cell  CD8 T cell         cDC       cMono      ncMono 
##          66         108         133          38         217         144 
##     NK cell         pDC Plasma cell  unassigned 
##         256           2           3         208
```

Then add the predictions to metadata and plot umap.


```r
# add in predictions
ctrl.sce@colData$scmap_cluster <- project_cluster$scmap_cluster_labs

plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "scmap_cluster")
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

## scMap cell
We can instead index the refernce data based on each single cell and project our data onto the closest neighbor in that dataset.


```r
ref.sce <- indexCell(ref.sce)
```


Again we need to index the reference dataset.


```r
project_cell <- scmapCell(projection = ctrl.sce, index_list = list(ref = metadata(ref.sce)$scmap_cell_index))
```

We now get a table with index for the 5 nearest neigbors in the reference dataset for each cell in our dataset.
We will select the celltype of the closest neighbor and assign it to the data.


```r
cell_type_pred <- colData(ref.sce)$cell_type1[project_cell$ref[[1]][1, ]]

table(cell_type_pred)
```

```
## cell_type_pred
##      B cell  CD4 T cell  CD8 T cell         cDC       cMono      ncMono 
##         100         184         289          52         221         151 
##     NK cell         pDC Plasma cell 
##         175           2           1
```


Then add the predictions to metadata and plot umap.


```r
# add in predictions
ctrl.sce@colData$scmap_cell <- cell_type_pred

plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "scmap_cell")
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Plot both:


```r
cowplot::plot_grid(ncol = 2, plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "scmap_cluster"), 
    plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "scmap_cell"))
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-17-1.png)<!-- -->



# scPred
scPred will train a classifier based on all principal components. First, `getFeatureSpace` will create a scPred object stored in the `@misc` slot where it extracts the PCs that best separates the different celltypes. Then `trainModel` will do the actual training for each celltype.

scPred works with Seurat objects, so we will convert both objects to seurat objects. You may see a lot of warnings about renaming things, but as long as you do not see an Error, you should be fine.


```r
suppressPackageStartupMessages(library(Seurat))

reference <- Seurat::as.Seurat(ref.sce)
ctrl <- Seurat::as.Seurat(ctrl.sce)
```

The loadings matrix is lost when converted to Seurat object, and scPred needs that information. So we need to rerun PCA with Seurat and the same hvgs.


```r
VariableFeatures(reference) = hvg.ref
reference <- reference %>% ScaleData(verbose = F) %>% RunPCA(verbose = F)

VariableFeatures(ctrl) = hvg.ctrl
ctrl <- ctrl %>% ScaleData(verbose = F) %>% RunPCA(verbose = F)
```




```r
reference <- getFeatureSpace(reference, "cell_type")
```

```
## ●  Extracting feature space for each cell type...
## DONE!
```

```r
reference <- trainModel(reference)
```

```
## ●  Training models for each cell type...
## DONE!
```



We can then print how well the training worked for the different celltypes by printing the number of PCs used for each, the ROC value and Sensitivity/Specificity. Which celltypes do you think are harder to classify based on this dataset?


```r
get_scpred(reference)
```

```
## 'scPred' object
## ✔  Prediction variable = cell_type 
## ✔  Discriminant features per cell type
## ✔  Training model(s)
## Summary
## 
## |Cell type   |    n| Features|Method    |   ROC|  Sens|  Spec|
## |:-----------|----:|--------:|:---------|-----:|-----:|-----:|
## |B cell      |  280|       50|svmRadial | 1.000| 1.000| 1.000|
## |CD4 T cell  | 1620|       50|svmRadial | 0.994| 0.972| 0.963|
## |CD8 T cell  |  945|       50|svmRadial | 0.973| 0.859| 0.971|
## |cDC         |   26|       50|svmRadial | 0.994| 0.727| 0.999|
## |cMono       |  212|       50|svmRadial | 1.000| 0.957| 0.997|
## |ncMono      |   79|       50|svmRadial | 1.000| 0.962| 0.999|
## |NK cell     |  312|       50|svmRadial | 0.998| 0.926| 0.995|
## |pDC         |   20|       50|svmRadial | 1.000| 0.950| 1.000|
## |Plasma cell |    6|       50|svmRadial | 1.000| 1.000| 1.000|
```

You can optimize parameters for each dataset by chaning parameters and testing different types of models, see more at: https://powellgenomicslab.github.io/scPred/articles/introduction.html. But for now, we will continue with this model.

 Now, lets predict celltypes on our data, where scPred will align the two datasets with Harmony and then perform classification.


```r
ctrl <- scPredict(ctrl, reference)
```

```
## ●  Matching reference with new dataset...
## 	 ─ 1000 features present in reference loadings
## 	 ─ 937 features shared between reference and new dataset
## 	 ─ 93.7% of features in the reference are present in new dataset
## ●  Aligning new data to reference...
## ●  Classifying cells...
## DONE!
```


```r
DimPlot(ctrl, group.by = "scpred_prediction", label = T, repel = T) + NoAxes()
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

Now plot how many	cells of each celltypes	can be found in	each cluster.


```r
ggplot(ctrl@meta.data, aes(x = louvain_SNNk15, fill = scpred_prediction)) + geom_bar() + 
    theme_classic()
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

Add the predictions into the SCE object


```r
ctrl.sce@colData$scpred_prediction = ctrl$scpred_prediction
```

# Compare results

Now we will compare the output of the two methods using the convenient function in scPred `crossTab` that prints the overlap between two metadata slots.


```r
crossTab(ctrl, "scmap_cell", "scpred_prediction")
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]},{"label":["Plasma cell"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"98","2":"4","3":"2","4":"0","5":"1","6":"0","7":"0","8":"0","9":"0","_rn_":"B cell"},{"1":"1","2":"138","3":"60","4":"0","5":"2","6":"1","7":"1","8":"0","9":"1","_rn_":"CD4 T cell"},{"1":"0","2":"28","3":"176","4":"0","5":"1","6":"1","7":"75","8":"0","9":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"6","5":"14","6":"1","7":"0","8":"0","9":"0","_rn_":"cDC"},{"1":"0","2":"0","3":"0","4":"44","5":"179","6":"47","7":"0","8":"0","9":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"2","5":"5","6":"97","7":"0","8":"0","9":"0","_rn_":"ncMono"},{"1":"0","2":"1","3":"32","4":"0","5":"0","6":"0","7":"86","8":"0","9":"0","_rn_":"NK cell"},{"1":"0","2":"0","3":"0","4":"0","5":"1","6":"0","7":"0","8":"1","9":"0","_rn_":"pDC"},{"1":"1","2":"1","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"Plasma cell"},{"1":"0","2":"12","3":"19","4":"0","5":"18","6":"4","7":"13","8":"1","9":"0","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>


# GSEA with celltype markers

Another option, where celltype can be classified on cluster level is to use gene set enrichment among the DEGs with known markers for different celltypes. Similar to how we did functional enrichment for the DEGs in the Differential expression exercise. 
There are some resources for celltype gene sets that can be used. Such as [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/), [PanglaoDB](https://panglaodb.se/) or celltype gene sets at [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp).
We can also look at overlap between DEGs in a reference dataset and the dataset you are analysing. 

## DEG overlap
First, lets extract top DEGs for our Covid-19 dataset and the reference dataset.
When we run differential expression for our dataset, we want to report as many genes as possible, hence we set the cutoffs quite lenient.


```r
# run differential expression in our dataset, using clustering at resolution 0.3
DGE_list <- scran::findMarkers(x = alldata, groups = as.character(alldata@colData$louvain_SNNk15), 
    pval.type = "all", min.prop = 0)
```


```r
# Compute differential gene expression in reference dataset (that has cell
# annotation)
ref_DGE <- scran::findMarkers(x = ref.sce, groups = as.character(ref.sce@colData$cell_type), 
    pval.type = "all", direction = "up")


# Identify the top cell marker genes in reference dataset select top 50 with
# hihgest foldchange among top 100 signifcant genes.
ref_list <- lapply(ref_DGE, function(x) {
    x$logFC <- rowSums(as.matrix(x[, grep("logFC", colnames(x))]))
    x %>% as.data.frame() %>% filter(p.value < 0.01) %>% top_n(-100, p.value) %>% 
        top_n(50, logFC) %>% rownames()
})

unlist(lapply(ref_list, length))
```

```
##      B cell  CD4 T cell  CD8 T cell         cDC       cMono      ncMono 
##          50          50          19          17          50          50 
##     NK cell         pDC Plasma cell 
##          50          50          24
```


Now we can run GSEA for the DEGs from our dataset and check for enrichment of top DEGs in the reference dataset.


```r
suppressPackageStartupMessages(library(fgsea))

# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    x$logFC <- rowSums(as.matrix(x[, grep("logFC", colnames(x))]))
    gene_rank <- setNames(x$logFC, rownames(x))
    fgseaRes <- fgsea(pathways = ref_list, stats = gene_rank, nperm = 10000)
    return(fgseaRes)
})
names(res) <- names(DGE_list)

# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
    x[x$pval < 0.1, ]
})
res <- lapply(res, function(x) {
    x[x$size > 2, ]
})
res <- lapply(res, function(x) {
    x[order(x$NES, decreasing = T), ]
})
res
```

```
## $`1`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:     B cell 0.0012683649 0.0022830568 -0.7787495 -1.424580           11   47
## 2:    NK cell 0.0008446838 0.0019005385 -0.7795612 -1.430296            7   49
## 3:        cDC 0.0045764036 0.0068646054 -0.8585156 -1.432512           40   17
## 4:     ncMono 0.0001055855 0.0003172589 -0.8307518 -1.524218            0   49
## 5:      cMono 0.0001057530 0.0003172589 -0.8452430 -1.543585            0   46
## 6: CD4 T cell 0.0001055186 0.0003172589 -0.9055773 -1.663893            0   50
##                                                 leadingEdge
## 1:                    RPS23,RPS11,CD52,RPL13A,RPL12,FAU,...
## 2:                    ITGB2,NKG7,MYO1F,GNLY,HCST,IFITM1,...
## 3: HLA-DRA,HLA-DRB1,HLA-DPB1,HLA-DPA1,HLA-DQB1,HLA-DRB5,...
## 4:                  S100A4,S100A11,AIF1,CEBPB,SAT1,PSAP,...
## 5:                   JUND,S100A6,TYROBP,LYZ,NFKBIA,FCN1,...
## 6:                 RPL34,RPS29,EEF1A1,RPL30,RPS14,RPL36,...
## 
## $`2`
##        pathway         pval         padj         ES       NES nMoreExtreme size
## 1:       cMono 0.0001638270 0.0005871608  0.9386094  1.899340            0   46
## 2:      ncMono 0.0035737492 0.0053606238  0.8156152  1.667098           21   49
## 3: Plasma cell 0.0704095720 0.0792107685 -0.7007685 -1.385564          305   24
## 4:         cDC 0.0331858407 0.0426675095 -0.7933131 -1.466581          149   17
## 5:  CD8 T cell 0.0011061947 0.0019911504 -0.9009140 -1.665500            4   17
## 6:     NK cell 0.0002600104 0.0005871608 -0.7893777 -1.772505            0   49
## 7:      B cell 0.0002585315 0.0005871608 -0.8666003 -1.927415            0   47
## 8:  CD4 T cell 0.0002609603 0.0005871608 -0.9369213 -2.111868            0   50
##                                                 leadingEdge
## 1:                  S100A8,S100A9,LYZ,S100A12,VCAN,RETN,...
## 2:               S100A11,S100A4,AIF1,SERPINA1,PSAP,MAFB,...
## 3:                  ISG20,PEBP1,RPL36AL,CYCS,FKBP11,MIF,...
## 4: HLA-DPB1,HLA-DPA1,HLA-DRB1,HLA-DRA,HLA-DQB1,HLA-DRB5,...
## 5:                         IL32,CCL5,GZMH,CD3D,CD2,LYAR,...
## 6:                       GNLY,B2M,GZMA,CTSW,IFITM1,NKG7,...
## 7:                    RPS5,RPL23A,CXCR4,CD52,CD37,RPS23,...
## 8:                   RPL14,RPL3,RPL5,RPS4X,EEF1A1,RPS3A,...
## 
## $`3`
##        pathway         pval         padj         ES       NES nMoreExtreme size
## 1:       cMono 0.0001154601 0.0005195705  0.9459979  1.735485            0   46
## 2:      ncMono 0.0001145607 0.0005195705  0.9069427  1.672099            0   49
## 3:         cDC 0.0095923261 0.0123329908  0.8887496  1.464616           71   17
## 4: Plasma cell 0.0650557621 0.0731877323 -0.6653311 -1.425870          139   24
## 5:  CD8 T cell 0.0004006410 0.0012019231 -0.9280132 -1.844796            0   17
## 6:     NK cell 0.0007855460 0.0012028869 -0.7946767 -1.955455            0   49
## 7:      B cell 0.0007616146 0.0012028869 -0.8352545 -2.039008            0   47
## 8:  CD4 T cell 0.0008019246 0.0012028869 -0.8743367 -2.156161            0   50
##                                                leadingEdge
## 1:                  S100A8,S100A9,LYZ,FCN1,VCAN,TYROBP,...
## 2:                AIF1,S100A11,FCER1G,PSAP,SAT1,S100A4,...
## 3: HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-DMA,HLA-DQB1,HLA-DPA1,...
## 4:                  PEBP1,ISG20,FKBP11,JCHAIN,MIF,CYCS,...
## 5:                        CCL5,IL32,GZMH,CD3D,CD2,CD8A,...
## 6:                     NKG7,GNLY,CST7,GZMA,CTSW,FGFBP2,...
## 7:                CXCR4,RPS5,RPL23A,RPL13A,CD52,RPL18A,...
## 8:                   RPL3,RPS4X,RPS3,RPS29,RPS27A,IL7R,...
## 
## $`4`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:     ncMono 0.0001088732 0.0009798585  0.9748472  1.754290            0   49
## 2:        cDC 0.0071063458 0.0127914225  0.8917037  1.461368           56   17
## 3:     B cell 0.0937500000 0.1406250000 -0.5050561 -1.273605           80   47
## 4: CD8 T cell 0.0005047956 0.0022715800 -0.9301044 -1.934617            0   17
## 5: CD4 T cell 0.0012578616 0.0028301887 -0.7732021 -1.973189            0   50
## 6:    NK cell 0.0012239902 0.0028301887 -0.8109740 -2.063669            0   49
##                                               leadingEdge
## 1:                 LST1,AIF1,COTL1,FCGR3A,FCER1G,PSAP,...
## 2: HLA-DPA1,HLA-DRA,HLA-DRB1,HLA-DPB1,HLA-DRB5,MTMR14,...
## 3:       CXCR4,MS4A1,RPL13A,TNFRSF13C,BANK1,LINC00926,...
## 4:                       CCL5,IL32,GZMH,CD3D,CD2,CD8A,...
## 5:                   LDHB,IL7R,RPS3,MGAT4A,RPL3,RPL31,...
## 6:                    NKG7,GNLY,CST7,GZMA,CTSW,FGFBP2,...
## 
## $`5`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:    NK cell 0.0002133561 0.0004800512  0.9811971  2.105527            0   49
## 2: CD8 T cell 0.0105574324 0.0135738417  0.8763643  1.602903           49   17
## 3: CD4 T cell 0.0067924528 0.0101886792 -0.7349319 -1.586954           35   50
## 4:        cDC 0.0028484618 0.0051272313 -0.8998497 -1.635120           14   17
## 5:     B cell 0.0001891432 0.0004800512 -0.8266821 -1.769288            0   47
## 6:     ncMono 0.0001881468 0.0004800512 -0.8449869 -1.822982            0   49
## 7:      cMono 0.0001895375 0.0004800512 -0.8961279 -1.910601            0   46
##                                                leadingEdge
## 1:                     GNLY,NKG7,FGFBP2,CST7,CTSW,PRF1,...
## 2:                   CCL5,GZMH,IL32,LYAR,CD2,LINC01871,...
## 3:                 RPS13,RPS28,TMEM123,IL7R,RPL22,TPT1,...
## 4: HLA-DRA,HLA-DRB1,HLA-DPA1,HLA-DQB1,HLA-DRB5,HLA-DMA,...
## 5:                  CD37,RPS11,RPL18A,CD52,RPL12,CD79B,...
## 6:                  COTL1,FTH1,AIF1,LST1,SAT1,SERPINA1,...
## 7:                  S100A9,S100A8,LYZ,FCN1,TKT,S100A12,...
## 
## $`6`
##        pathway         pval         padj         ES       NES nMoreExtreme size
## 1:       cMono 0.0001850139 0.0009193054  0.8712632  2.007486            0   46
## 2:      ncMono 0.0003677823 0.0011033468  0.8215868  1.912399            1   49
## 3:         cDC 0.0002042901 0.0009193054  0.9619606  1.888309            0   17
## 4:         pDC 0.0053624260 0.0080436391  0.7616438  1.760837           28   47
## 5:      B cell 0.0015237266 0.0027427079 -0.6804223 -1.707087            6   47
## 6: Plasma cell 0.0012150668 0.0027339004 -0.7864297 -1.741016            5   24
##                                                 leadingEdge
## 1:                  LYZ,TYROBP,S100A9,S100A6,FCN1,AP1S2,...
## 2:                   AIF1,COTL1,FCER1G,S100A4,LST1,PSAP,...
## 3: HLA-DRA,HLA-DPA1,HLA-DPB1,HLA-DRB1,HLA-DRB5,HLA-DQB1,...
## 4:              NPC2,PTPRE,PPP1R14B,PLD4,UNC93B1,BCL11A,...
## 5:              CXCR4,CD79B,MS4A1,TNFRSF13C,BIRC3,BANK1,...
## 6:                    ISG20,CYCS,RABAC1,SUB1,DAD1,SPCS2,...
## 
## $`7`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1: CD4 T cell 0.0002018571 0.0006055713  0.9794603  2.094298            0   50
## 2:     B cell 0.0212467428 0.0273172408  0.7117450  1.507914          105   47
## 3:    NK cell 0.0106782677 0.0160174016 -0.7219152 -1.556888           53   49
## 4:        cDC 0.0021471794 0.0038649229 -0.9031516 -1.643987           10   17
## 5:        pDC 0.0003989627 0.0008976661 -0.7998123 -1.707668            1   47
## 6:      cMono 0.0001997204 0.0006055713 -0.9131339 -1.941423            0   46
## 7:     ncMono 0.0001977457 0.0006055713 -0.9440641 -2.035976            0   49
##                                                 leadingEdge
## 1:                    IL7R,LDHB,RPL3,RPS6,RPS29,PIK3IP1,...
## 2:                RPS5,RPL13A,RPL23A,RPL18A,RPS23,CXCR4,...
## 3:                    NKG7,GNLY,ITGB2,MYO1F,CST7,FGFBP2,...
## 4: HLA-DRA,HLA-DRB1,HLA-DPA1,HLA-DPB1,HLA-DQB1,HLA-DRB5,...
## 5:                      PLEK,NPC2,CTSB,PTPRE,PLAC8,IRF8,...
## 6:                 S100A9,S100A8,LYZ,TYROBP,FCN1,S100A6,...
## 7:                 FCER1G,AIF1,PSAP,LST1,IFITM3,S100A11,...
## 
## $`8`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:     B cell 0.0002023063 0.0004551892  0.9639514  2.034410            0   47
## 2: CD4 T cell 0.0002015723 0.0004551892  0.8788376  1.876335            0   50
## 3:        cDC 0.0006088898 0.0009133347  0.9453704  1.696364            2   17
## 4: CD8 T cell 0.0021674877 0.0027867699 -0.8993881 -1.623566           10   17
## 5:      cMono 0.0005931198 0.0009133347 -0.8102584 -1.712948            2   46
## 6:     ncMono 0.0001981768 0.0004551892 -0.9052455 -1.933940            0   49
## 7:    NK cell 0.0001981768 0.0004551892 -0.9092964 -1.942594            0   49
##                                                leadingEdge
## 1:          MS4A1,CD37,CXCR4,TNFRSF13C,BANK1,LINC00926,...
## 2:                   RPS6,RPS29,RPL13,RPL32,RPS3A,RPL3,...
## 3: HLA-DRA,HLA-DQB1,HLA-DPB1,HLA-DRB1,HLA-DPA1,HLA-DMA,...
## 4:                        CCL5,IL32,GZMH,CD3D,CD2,LYAR,...
## 5:                S100A6,S100A9,LYZ,S100A8,TYROBP,FCN1,...
## 6:                S100A4,FCER1G,S100A11,AIF1,LST1,PSAP,...
## 7:                     ITGB2,NKG7,HCST,GNLY,MYO1F,CST7,...
## 
## $`9`
##        pathway         pval         padj         ES       NES nMoreExtreme size
## 1:     NK cell 0.0001712622 0.0005404756  0.9213333  1.918701            0   49
## 2:  CD4 T cell 0.0001706485 0.0005404756  0.9032275  1.888260            0   50
## 3:  CD8 T cell 0.0005313496 0.0009564293  0.9628247  1.708426            2   17
## 4: Plasma cell 0.0415212840 0.0533845080  0.7780907  1.456961          237   24
## 5:         cDC 0.0064279155 0.0096418733 -0.8910223 -1.638995           27   17
## 6:      ncMono 0.0002402114 0.0005404756 -0.9146086 -2.018853            0   49
## 7:       cMono 0.0002400960 0.0005404756 -0.9259084 -2.023623            0   46
##                                             leadingEdge
## 1:                    NKG7,CST7,GZMA,GNLY,GZMM,CTSW,...
## 2:                 RPS29,RPS3,RPL3,IL7R,RPL14,RPS4X,...
## 3:                    CCL5,IL32,GZMH,CD3D,CD8A,LYAR,...
## 4:             RPL36AL,FKBP11,PPIB,ISG20,PEBP1,CYCS,...
## 5: HLA-DRA,HLA-DRB5,HLA-DMA,BASP1,HLA-DQB1,HLA-DRB1,...
## 6:                 FCER1G,AIF1,LST1,COTL1,FTH1,PSAP,...
## 7:            S100A9,S100A8,LYZ,TYROBP,FCN1,S100A12,...
```

Selecing top significant overlap per cluster, we can now rename the clusters according to the predicted labels. OBS! Be aware that if you have some clusters that have bad p-values for all the gene sets, the cluster label will not be very reliable. Also, the gene sets you are using may not cover all the celltypes you have in your dataset and hence predictions may just be the most similar celltype.
Also, some of the clusters have very similar p-values to multiple celltypes, for instance the ncMono and cMono celltypes are equally good for some clusters.


```r
new.cluster.ids <- unlist(lapply(res, function(x) {
    as.data.frame(x)[1, 1]
}))

alldata@colData$ref_gsea <- new.cluster.ids[as.character(alldata@colData$louvain_SNNk15)]

cowplot::plot_grid(ncol = 2, plotReducedDim(alldata, dimred = "UMAP", colour_by = "louvain_SNNk15"), 
    plotReducedDim(alldata, dimred = "UMAP", colour_by = "ref_gsea"))
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

Compare to results with the other celltype prediction methods in the ctrl_13 sample.


```r
ctrl.sce@colData$ref_gsea = alldata@colData$ref_gsea[alldata@colData$sample == "ctrl.13"]

cowplot::plot_grid(ncol = 3, plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "ref_gsea"), 
    plotReducedDim(ctrl.sce, dimred = "UMAP", colour_by = "scmap_cell"), plotReducedDim(ctrl.sce, 
        dimred = "UMAP", colour_by = "scpred_prediction"))
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

## With annotated gene sets
First download celltype gene sets from CellMarker.


```r
# Download gene marker list
if (!dir.exists("data/CellMarker_list/")) {
    dir.create("data/CellMarker_list")
    download.file(url = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt", 
        destfile = "./data/CellMarker_list/Human_cell_markers.txt")
    download.file(url = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt", 
        destfile = "./data/CellMarker_list/Mouse_cell_markers.txt")
    
}
```

Read in the gene lists and do some filering.


```r
# Load the human marker table
markers <- read.delim("data/CellMarker_list/Human_cell_markers.txt")
markers <- markers[markers$speciesType == "Human", ]
markers <- markers[markers$cancerType == "Normal", ]

# Filter by tissue (to reduce computational time and have tissue-specific
# classification) sort(unique(markers$tissueType))
# grep('blood',unique(markers$tissueType),value = T) markers <- markers [
# markers$tissueType %in% c('Blood','Venous blood', 'Serum','Plasma',
# 'Spleen','Bone marrow','Lymph node'), ]


# remove strange characters etc.
celltype_list <- lapply(unique(markers$cellName), function(x) {
    x <- paste(markers$geneSymbol[markers$cellName == x], sep = ",")
    x <- gsub("[[]|[]]| |-", ",", x)
    x <- unlist(strsplit(x, split = ","))
    x <- unique(x[!x %in% c("", "NA", "family")])
    x <- casefold(x, upper = T)
})
names(celltype_list) <- unique(markers$cellName)
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]} )
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) < 100]
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) > 5]
```


```r
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    x$logFC <- rowSums(as.matrix(x[, grep("logFC", colnames(x))]))
    gene_rank <- setNames(x$logFC, rownames(x))
    fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank, nperm = 10000)
    return(fgseaRes)
})
names(res) <- names(DGE_list)


# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {
    x[x$pval < 0.01, ]
})
res <- lapply(res, function(x) {
    x[x$size > 5, ]
})
res <- lapply(res, function(x) {
    x[order(x$NES, decreasing = T), ]
})

# show top 3 for each cluster.
lapply(res, head, 3)
```

```
## $`1`
##                              pathway        pval       padj         ES
## 1: Effector CD8+ memory T (Tem) cell 0.009267841 0.09679745 -0.6935858
## 2:                  Mesenchymal cell 0.005661564 0.09342601 -0.7266267
## 3:               Natural killer cell 0.001642373 0.07073592 -0.7179163
##          NES nMoreExtreme size                             leadingEdge
## 1: -1.317498           89   79 LGALS1,ZEB2,EFHD2,AHNAK,MYO1F,ARL4C,...
## 2: -1.351676           53   59     S100A4,VIM,PTPRC,CD44,ZEB2,CTSC,...
## 3: -1.368422           15   84     PTPRC,NKG7,GNLY,FCGR3A,ID2,CD81,...
## 
## $`2`
##                            pathway         pval       padj        ES      NES
## 1:                      Neutrophil 0.0001534213 0.01518333 0.9056965 1.964839
## 2:          CD1C+_B dendritic cell 0.0001615248 0.01518333 0.9195692 1.891113
## 3: Monocyte derived dendritic cell 0.0032667877 0.09146749 0.9151496 1.608745
##    nMoreExtreme size                                 leadingEdge
## 1:            0   80 S100A8,S100A9,S100A12,MNDA,CD14,S100A11,...
## 2:            0   53     S100A8,S100A9,LYZ,S100A12,VCAN,RETN,...
## 3:           17   17     S100A8,S100A9,CD14,CST3,ITGAM,SIRPA,...
## 
## $`3`
##                   pathway         pval       padj        ES      NES
## 1:             Neutrophil 0.0001084481 0.01079095 0.9155801 1.749205
## 2: CD1C+_B dendritic cell 0.0001147974 0.01079095 0.9054218 1.673682
## 3:           Stromal cell 0.0002388630 0.01496875 0.8928301 1.600967
##    nMoreExtreme size                                 leadingEdge
## 1:            0   80 S100A8,S100A9,S100A12,S100A11,CD14,MNDA,...
## 2:            0   53     S100A8,S100A9,LYZ,FCN1,VCAN,S100A12,...
## 3:            1   38          VIM,CD44,TIMP2,ICAM1,BST1,KLF6,...
## 
## $`4`
##             pathway         pval       padj        ES      NES nMoreExtreme
## 1: Mesenchymal cell 0.0004295533 0.04040694 0.8355991 1.519021            3
## 2:     Stromal cell 0.0011252391 0.07051498 0.8461284 1.496026            9
## 3:    Hemangioblast 0.0004298610 0.04040694 0.9901971 1.484564            2
##    size                           leadingEdge
## 1:   59   COTL1,S100A4,CTSC,ZEB2,HES4,VIM,...
## 2:   38 PECAM1,TIMP1,VIM,CD44,TIMP2,ICAM3,...
## 3:    7                           PECAM1,CD34
## 
## $`5`
##                              pathway         pval        padj        ES
## 1:             CD4+ cytotoxic T cell 0.0002158895 0.008117444 0.9360306
## 2: Effector CD8+ memory T (Tem) cell 0.0002158895 0.008117444 0.8870134
## 3:               Natural killer cell 0.0002156567 0.008117444 0.8375716
##         NES nMoreExtreme size                           leadingEdge
## 1: 2.171691            0   86   GNLY,NKG7,CCL5,GZMB,FGFBP2,CST7,...
## 2: 2.036462            0   79 GNLY,GZMB,FGFBP2,GZMH,KLRD1,SPON2,...
## 3: 1.939301            0   84   GNLY,NKG7,GZMB,GZMA,CD247,KLRD1,...
## 
## $`6`
##                            pathway         pval       padj        ES      NES
## 1:          CD1C+_B dendritic cell 0.0003613369 0.03396567 0.8367127 1.950614
## 2:                     Acinar cell 0.0003601657 0.03396567 0.8250723 1.927913
## 3: Monocyte derived dendritic cell 0.0018718802 0.08797837 0.9131565 1.771746
##    nMoreExtreme size                              leadingEdge
## 1:            1   53 LYZ,S100A9,FCN1,ANXA1,FCER1A,CLEC10A,...
## 2:            1   54    LYZ,LGALS2,EFHD2,NFKBIA,PPIF,YBX3,...
## 3:            8   17 CST3,S100A9,FCER1A,S100A8,CD1C,ITGAX,...
## 
## $`7`
##              pathway         pval        padj        ES      NES nMoreExtreme
## 1: Naive CD8+ T cell 0.0001973554 0.006414193 0.8482016 1.988338            0
## 2: Naive CD4+ T cell 0.0002047083 0.006414193 0.9229411 1.856407            0
## 3:       CD4+ T cell 0.0004111842 0.011043233 0.9157997 1.751179            1
##    size                           leadingEdge
## 1:   91 LDHB,PIK3IP1,NPM1,RPS8,TCF7,NOSIP,...
## 2:   34  IL7R,RPS5,EEF1B2,TCF7,NOSIP,CCR7,...
## 3:   25       IL7R,LTB,CD3E,CD3D,CD3G,CD2,...
## 
## $`8`
##                        pathway        pval       padj         ES       NES
## 1:              Pyramidal cell 0.003343824 0.04190926 -0.9708267 -1.503295
## 2: CD4+CD25+ regulatory T cell 0.001573564 0.02689364 -0.9789523 -1.515877
## 3:      Lake et al.Science.Ex4 0.006465517 0.06386890 -0.9434984 -1.521354
##    nMoreExtreme size              leadingEdge
## 1:           16    6                NRGN,CD3E
## 2:            7    6 CD3E,CD3D,PTPRC,CD3G,CD4
## 3:           32    8                    ANXA1
## 
## $`9`
##                  pathway         pval       padj        ES      NES
## 1: CD4+ cytotoxic T cell 0.0001618647 0.01217222 0.8595021 1.940950
## 2:           CD8+ T cell 0.0001780944 0.01217222 0.9501193 1.726859
## 3:   Natural killer cell 0.0003237294 0.01217222 0.7659932 1.725979
##    nMoreExtreme size                       leadingEdge
## 1:            0   86 CCL5,NKG7,GZMH,CST7,GZMA,GNLY,...
## 2:            0   19 NKG7,CD3D,CD3E,CD3G,CD8A,GZMK,...
## 3:            1   84 NKG7,CD3D,CD3E,GZMA,GNLY,CD3G,...
```

#CT_GSEA8:


```r
new.cluster.ids <- unlist(lapply(res, function(x) {
    as.data.frame(x)[1, 1]
}))
alldata@colData$cellmarker_gsea <- new.cluster.ids[as.character(alldata@colData$louvain_SNNk15)]

cowplot::plot_grid(ncol = 2, plotReducedDim(alldata, dimred = "UMAP", colour_by = "cellmarker_gsea"), 
    plotReducedDim(alldata, dimred = "UMAP", colour_by = "ref_gsea"))
```

![](scater_06_celltype_files/figure-html/unnamed-chunk-35-1.png)<!-- -->




Do you think that the methods overlap well? Where do you see the most inconsistencies?

In this case we do not have any ground truth, and we cannot say wich method performs best. You should keep in mind, that any celltype classification method is just a prediction, and you still need to use your common sense and knowledge of the biological system to judge if the results make sense.

# Save data
Finally, lets save the data with predictions.


```r
saveRDS(ctrl.sce, "data/results/ctrl13_qc_dr_int_cl_celltype.rds")
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] fgsea_1.12.0                Rcpp_1.0.6                 
##  [3] caret_6.0-86                lattice_0.20-41            
##  [5] Seurat_3.2.3                scmap_1.8.0                
##  [7] scPred_1.9.0                rafalib_1.0.0              
##  [9] pheatmap_1.0.12             cowplot_1.1.1              
## [11] dplyr_1.0.3                 scran_1.14.1               
## [13] scater_1.14.0               ggplot2_3.3.3              
## [15] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
## [17] DelayedArray_0.12.0         BiocParallel_1.20.1        
## [19] matrixStats_0.57.0          Biobase_2.46.0             
## [21] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [23] IRanges_2.20.2              S4Vectors_0.24.4           
## [25] BiocGenerics_0.32.0         RJSONIO_1.3-1.4            
## [27] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] reticulate_1.18          tidyselect_1.1.0         htmlwidgets_1.5.3       
##   [4] grid_3.6.1               Rtsne_0.15               pROC_1.17.0.1           
##   [7] munsell_0.5.0            codetools_0.2-18         ica_1.0-2               
##  [10] statmod_1.4.35           future_1.21.0            miniUI_0.1.1.1          
##  [13] withr_2.4.0              colorspace_2.0-0         highr_0.8               
##  [16] knitr_1.30               ROCR_1.0-11              tensor_1.5              
##  [19] listenv_0.8.0            labeling_0.4.2           GenomeInfoDbData_1.2.2  
##  [22] harmony_1.0              polyclip_1.10-0          farver_2.0.3            
##  [25] parallelly_1.23.0        vctrs_0.3.6              generics_0.1.0          
##  [28] ipred_0.9-9              xfun_0.20                randomForest_4.6-14     
##  [31] R6_2.5.0                 ggbeeswarm_0.6.0         rsvd_1.0.3              
##  [34] locfit_1.5-9.4           bitops_1.0-6             spatstat.utils_1.20-2   
##  [37] assertthat_0.2.1         promises_1.1.1           scales_1.1.1            
##  [40] nnet_7.3-14              beeswarm_0.2.3           gtable_0.3.0            
##  [43] globals_0.14.0           goftest_1.2-2            timeDate_3043.102       
##  [46] rlang_0.4.10             splines_3.6.1            lazyeval_0.2.2          
##  [49] ModelMetrics_1.2.2.2     yaml_2.2.1               reshape2_1.4.4          
##  [52] abind_1.4-5              httpuv_1.5.5             tools_3.6.1             
##  [55] lava_1.6.8.1             ellipsis_0.3.1           RColorBrewer_1.1-2      
##  [58] proxy_0.4-24             ggridges_0.5.3           plyr_1.8.6              
##  [61] zlibbioc_1.32.0          purrr_0.3.4              RCurl_1.98-1.2          
##  [64] rpart_4.1-15             deldir_0.2-9             pbapply_1.4-3           
##  [67] viridis_0.5.1            zoo_1.8-8                ggrepel_0.9.1           
##  [70] cluster_2.1.0            magrittr_2.0.1           data.table_1.13.6       
##  [73] RSpectra_0.16-0          scattermore_0.7          lmtest_0.9-38           
##  [76] RANN_2.6.1               fitdistrplus_1.1-3       patchwork_1.1.1         
##  [79] mime_0.9                 evaluate_0.14            xtable_1.8-4            
##  [82] gridExtra_2.3            compiler_3.6.1           tibble_3.0.5            
##  [85] KernSmooth_2.23-18       crayon_1.3.4             htmltools_0.5.1         
##  [88] mgcv_1.8-33              later_1.1.0.1            tidyr_1.1.2             
##  [91] lubridate_1.7.9.2        DBI_1.1.1                formatR_1.7             
##  [94] MASS_7.3-53              Matrix_1.3-2             getopt_1.20.3           
##  [97] cli_2.2.0                gower_0.2.2              igraph_1.2.6            
## [100] pkgconfig_2.0.3          plotly_4.9.3             recipes_0.1.15          
## [103] foreach_1.5.1            vipor_0.4.5              dqrng_0.2.1             
## [106] XVector_0.26.0           prodlim_2019.11.13       stringr_1.4.0           
## [109] digest_0.6.27            sctransform_0.3.2        RcppAnnoy_0.0.18        
## [112] spatstat.data_1.7-0      fastmatch_1.1-0          rmarkdown_2.6           
## [115] leiden_0.3.6             uwot_0.1.10              edgeR_3.28.1            
## [118] DelayedMatrixStats_1.8.0 googleVis_0.6.9          kernlab_0.9-29          
## [121] shiny_1.5.0              lifecycle_0.2.0          nlme_3.1-150            
## [124] jsonlite_1.7.2           BiocNeighbors_1.4.0      fansi_0.4.2             
## [127] viridisLite_0.3.0        limma_3.42.2             pillar_1.4.7            
## [130] fastmap_1.0.1            httr_1.4.2               survival_3.2-7          
## [133] glue_1.4.2               FNN_1.1.3                spatstat_1.64-1         
## [136] png_0.1-7                iterators_1.0.13         class_7.3-17            
## [139] stringi_1.5.3            BiocSingular_1.2.0       irlba_2.3.3             
## [142] e1071_1.7-4              future.apply_1.7.0
```



