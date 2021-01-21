---
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 21, 2021'
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
ctrl.sce <- selectFeatures(ctrl.sce, suppress_plot = TRUE)
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
##          68         107         127          35         200         140 
##     NK cell         pDC Plasma cell  unassigned 
##         254           2           3         200
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
##         101         180         318          39         205         158 
##     NK cell         pDC Plasma cell 
##         132           1           2
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
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]},{"label":["Plasma cell"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"100","2":"0","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"B cell"},{"1":"0","2":"144","3":"52","4":"0","5":"0","6":"2","7":"1","8":"0","9":"0","_rn_":"CD4 T cell"},{"1":"0","2":"25","3":"193","4":"0","5":"0","6":"1","7":"51","8":"0","9":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"12","5":"5","6":"0","7":"0","8":"0","9":"0","_rn_":"cDC"},{"1":"0","2":"0","3":"0","4":"23","5":"180","6":"51","7":"0","8":"0","9":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"0","5":"7","6":"94","7":"0","8":"0","9":"0","_rn_":"ncMono"},{"1":"0","2":"0","3":"51","4":"0","5":"0","6":"0","7":"70","8":"0","9":"0","_rn_":"NK cell"},{"1":"0","2":"0","3":"0","4":"1","5":"0","6":"0","7":"0","8":"1","9":"0","_rn_":"pDC"},{"1":"0","2":"0","3":"0","4":"0","5":"0","6":"1","7":"0","8":"0","9":"1","_rn_":"Plasma cell"},{"1":"1","2":"11","3":"22","4":"3","5":"13","6":"9","7":"10","8":"0","9":"1","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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
ref_list <- lapply(ref_DGE, function(x) x %>% as.data.frame() %>% filter(p.value < 
    0.01) %>% top_n(-100, p.value) %>% top_n(50, summary.logFC) %>% rownames())

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
    gene_rank <- setNames(x$summary.logFC, rownames(x))
    fgseaRes <- fgsea(pathways = ref_list, stats = gene_rank)
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
## Empty data.table (0 rows and 8 cols): pathway,pval,padj,log2err,ES,NES...
## 
## $`2`
##        pathway       pval     padj   log2err         ES       NES size
## 1:       cMono 0.01365441 0.110540 0.3807304  0.8238865  1.615964   46
## 2:  CD4 T cell 0.09663866 0.289916 0.2089550  0.6841648  1.349772   49
## 3: Plasma cell 0.02456446 0.110540 0.3524879 -0.8144882 -1.427707   24
##                                  leadingEdge
## 1: S100A8,RETN,S100A9,S100A12,PLBD1,JUND,...
## 2:    RPL34,RPS14,RPL13,RPL36,RPL9,RPL38,...
## 3:     CYCS,SPCS2,SUB1,RPL36AL,ISG20,MIF,...
## 
## $`3`
##        pathway        pval       padj   log2err         ES       NES size
## 1:  CD8 T cell 0.008473208 0.03812944 0.3807304  0.9079970  1.597666   17
## 2:         cDC 0.036041792 0.10812538 0.3217759  0.8431468  1.483559   17
## 3: Plasma cell 0.097719870 0.18392371 0.2616635  0.7363860  1.365970   24
## 4:  CD4 T cell 0.005047935 0.03812944 0.4070179 -0.7892089 -1.460559   49
##                                             leadingEdge
## 1:                     CCL5,IL32,GZMH,CD3D,CD2,CD8A,...
## 2: HLA-DRA,FCER1A,HLA-DRB5,HLA-DRB1,HLA-DMA,CLEC10A,...
## 3:           DAD1,JCHAIN,MZB1,RPL36AL,TNFRSF17,SUB1,...
## 4:              RPL34,RPL22,RPS14,RPL31,RPL13,RPL36,...
## 
## $`4`
##       pathway         pval         padj   log2err         ES       NES size
## 1:     B cell 4.383478e-07 3.945130e-06 0.6749629  0.9688540  1.736764   46
## 2:        cDC 2.194401e-02 8.206675e-02 0.3524879  0.9086598  1.467294   17
## 3:        pDC 5.539359e-02 1.246356e-01 0.2311267  0.7998316  1.436890   47
## 4:     ncMono 9.385113e-02 1.407767e-01 0.2663507 -0.6733879 -1.339901   49
## 5: CD4 T cell 8.090615e-02 1.407767e-01 0.2878571 -0.6848419 -1.362692   49
## 6:    NK cell 2.735558e-02 8.206675e-02 0.3524879 -0.7648473 -1.521885   49
##                                            leadingEdge
## 1:    MS4A1,TNFRSF13C,LINC00926,BANK1,RALGPS2,CD37,...
## 2: HLA-DQB1,HLA-DRA,HLA-DPB1,HLA-DRB1,HLA-DQA2,HLA-DMA
## 3:              IRF8,TCF4,BCL11A,SPIB,YPEL5,CCDC50,...
## 4:             TIMP1,S100A4,S100A11,COTL1,MT2A,BID,...
## 5:              RPL5,RPS25,RPL34,RPS28,RPL36,RPL13,...
## 6:                  B2M,HCST,NKG7,BIN2,ITGB2,MYO1F,...
## 
## $`5`
##    pathway         pval        padj   log2err       ES      NES size
## 1:  ncMono 2.754733e-10 2.47926e-09 0.8140358 0.972421 1.564507   49
##                                 leadingEdge
## 1: COTL1,CDKN1C,LST1,RHOC,FCGR3A,SMIM25,...
## 
## $`6`
##    pathway         pval         padj   log2err         ES       NES size
## 1: NK cell 0.0000348878 0.0003139902 0.5573322  0.9443137  1.771130   49
## 2:  B cell 0.0864197531 0.1944444444 0.1882041  0.7434475  1.378395   46
## 3:  ncMono 0.0319362877 0.0958088631 0.3217759 -0.7177214 -1.495263   49
## 4:     cDC 0.0251726572 0.0958088631 0.3524879 -0.8438787 -1.517539   17
##                                   leadingEdge
## 1:     GNLY,SPON2,KLRF1,PRF1,FGFBP2,CD247,...
## 2:           CXCR4,CD37,FAU,RPL23A,RPS5,BIRC3
## 3:   COTL1,FCGR3A,FTH1,IFITM3,PSAP,NAP1L1,...
## 4: HLA-DRA,HLA-DPA1,HLA-DMA,HLA-DRB5,HLA-DPB1
## 
## $`7`
##       pathway         pval        padj   log2err        ES      NES size
## 1: CD8 T cell 0.0002959415 0.002663474 0.4984931 0.9781262 1.634839   17
##                           leadingEdge
## 1: CD8A,GZMK,LYAR,KLRG1,CD8B,CD3D,...
## 
## $`8`
##       pathway         pval         padj   log2err         ES       NES size
## 1: CD4 T cell 1.284929e-09 1.156436e-08 0.7881868  0.9604161  1.800741   49
## 2:        pDC 4.405108e-02 1.321532e-01 0.3217759 -0.6987740 -1.466630   47
## 3:        cDC 8.831471e-03 3.974162e-02 0.3807304 -0.8958972 -1.597727   17
##                                         leadingEdge
## 1:           TPT1,LDHB,IL7R,NOSIP,RCAN3,PIK3IP1,...
## 2:         PLEK,C12orf75,PTPRE,PLAC8,PARK7,CHD9,...
## 3: HLA-DRA,HLA-DRB1,HLA-DQB1,HLA-DMA,HLA-DRB5,GPAT3
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
    gene_rank <- setNames(x$summary.logFC, rownames(x))
    fgseaRes <- fgsea(pathways = celltype_list, stats = gene_rank)
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
##                    pathway         pval       padj   log2err        ES      NES
## 1:             Acinar cell 0.0002748268 0.05166744 0.4984931 0.8450115 1.623988
## 2: Hematopoietic stem cell 0.0010512492 0.09881743 0.4550599 0.8331119 1.596825
## 3:  Lake et al.Science.Ex4 0.0017232451 0.10799003 0.4550599 0.9825751 1.455149
##    size                           leadingEdge
## 1:   54 LGALS2,SGK1,SOD2,PPIF,YBX3,NFKBIA,...
## 2:   53 MCL1,SPI1,CD44,ALDH1A1,CD33,ITGA5,...
## 3:    8                        ANXA1,HS3ST3B1
## 
## $`2`
##                         pathway        pval      padj   log2err         ES
## 1:       CD1C+_B dendritic cell 0.002208939 0.1190423 0.4317077  0.8307889
## 2: Hematopoietic precursor cell 0.003047609 0.1190423 0.4317077 -0.9830610
## 3: Red blood cell (erythrocyte) 0.003047609 0.1190423 0.4317077 -0.9831977
##          NES size                                 leadingEdge
## 1:  1.629018   53 S100A8,RETN,S100A9,S100A12,PLBD1,RNASE2,...
## 2: -1.423192    6                                  CD14,PTPRC
## 3: -1.423390    6                                 PTPRC,ITGB3
## 
## $`3`
##                   pathway         pval        padj   log2err        ES      NES
## 1:          Megakaryocyte 3.594429e-05 0.006757526 0.5573322 0.9395377 1.779347
## 2:    Embryonic stem cell 6.318565e-03 0.252244372 0.4070179 0.8635976 1.602579
## 3: Lake et al.Science.Ex4 4.592561e-03 0.252244372 0.4070179 0.9721453 1.520334
##    size                            leadingEdge
## 1:   26       PF4,PPBP,GP9,ITGA2B,FYB1,CD9,...
## 2:   24 CD9,PECAM1,ITGB1,LEFTY2,ITGA6,FUT4,...
## 3:    8                           ANXA1,COL6A1
## 
## $`4`
##                         pathway         pval       padj   log2err        ES
## 1:            Follicular B cell 0.0009562224 0.06950636 0.4772708 0.9597889
## 2: Megakaryocyte erythroid cell 0.0084678667 0.11371135 0.3807304 0.8271617
## 3:           Hematopoietic cell 0.0025880027 0.06950636 0.4317077 0.9522548
##         NES size                          leadingEdge
## 1: 1.590311   22  MS4A1,CD22,FCER2,CD40,CD69,PAX5,...
## 2: 1.579933   83 CD79A,CD83,FCER2,CD69,PTPRC,CD81,...
## 3: 1.559672   19                MS4A1,PTPRC,CD19,CD38
## 
## $`5`
##                         pathway        pval      padj   log2err        ES
## 1: Megakaryocyte erythroid cell 0.001324355 0.1244894 0.4550599 0.8267372
## 2:                Lymphoid cell 0.008528184 0.2266603 0.3807304 0.8424409
## 3:                Hemangioblast 0.002583281 0.1523527 0.4317077 0.9872838
##         NES size                           leadingEdge
## 1: 1.393082   83 FCGR3A,PECAM1,SPN,CD68,ITGAX,CD86,...
## 2: 1.386579   50  FCGR3A,SPN,CD68,ITGAX,IL17RA,CD4,...
## 3: 1.386466    7                           PECAM1,CD34
## 
## $`6`
##                              pathway        pval       padj   log2err        ES
## 1:             CD4+ cytotoxic T cell 0.001263875 0.08907058 0.4550599 0.8213135
## 2: Effector CD8+ memory T (Tem) cell 0.009466795 0.19502407 0.3807304 0.7937377
## 3:                     M2 macrophage 0.007443650 0.19502407 0.4070179 0.9263520
##         NES size                            leadingEdge
## 1: 1.649268   86   GZMB,GNLY,SPON2,PRF1,FGFBP2,CTSW,...
## 2: 1.577288   79 GZMB,GNLY,SPON2,KLRF1,FGFBP2,KLRD1,...
## 3: 1.557114   21                    CXCR4,TNFSF10,STAT6
## 
## $`7`
##                 pathway        pval      padj   log2err        ES      NES size
## 1:        T helper cell 0.002558426 0.1413589 0.4317077 0.8581622 1.633830   56
## 2:   Hematopoietic cell 0.003065630 0.1413589 0.4317077 0.9345773 1.608948   19
## 3: T helper2 (Th2) cell 0.006004844 0.1413589 0.4070179 0.8953885 1.606361   27
##                           leadingEdge
## 1: IL7R,GZMK,CD3D,CD3G,KLRB1,CD3E,...
## 2:    CD8A,CD3D,CD3G,CD3E,PTPRC,ITGAM
## 3: CD3D,CD3G,CXCR4,CD3E,MAF,GATA3,...
## 
## $`8`
##              pathway         pval        padj   log2err        ES      NES size
## 1: Naive CD8+ T cell 8.748345e-05 0.008223444 0.5384341 0.8332577 1.645418   91
## 2: Naive CD4+ T cell 1.869176e-03 0.058567529 0.4550599 0.8811387 1.584745   34
## 3:     M1 macrophage 1.346206e-03 0.050617339 0.4550599 0.9161355 1.573322   24
##                              leadingEdge
## 1: LDHB,NOSIP,RCAN3,PIK3IP1,TCF7,MAL,...
## 2: IL7R,NOSIP,TCF7,MAL,PRKCA,TRABD2A,...
## 3:      IL7R,CCR7,TSPO,SELL,IL2RA,IL17RA
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] fgsea_1.16.0                caret_6.0-86               
##  [3] lattice_0.20-41             Seurat_3.2.3               
##  [5] scmap_1.12.0                scPred_1.9.0               
##  [7] rafalib_1.0.0               pheatmap_1.0.12            
##  [9] cowplot_1.1.1               dplyr_1.0.3                
## [11] scran_1.18.0                scater_1.18.0              
## [13] ggplot2_3.3.3               SingleCellExperiment_1.12.0
## [15] SummarizedExperiment_1.20.0 Biobase_2.50.0             
## [17] GenomicRanges_1.42.0        GenomeInfoDb_1.26.0        
## [19] IRanges_2.24.0              S4Vectors_0.28.0           
## [21] BiocGenerics_0.36.0         MatrixGenerics_1.2.0       
## [23] matrixStats_0.57.0          RJSONIO_1.3-1.4            
## [25] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] reticulate_1.18           tidyselect_1.1.0         
##   [3] htmlwidgets_1.5.3         grid_4.0.3               
##   [5] BiocParallel_1.24.0       Rtsne_0.15               
##   [7] pROC_1.17.0.1             munsell_0.5.0            
##   [9] codetools_0.2-18          ica_1.0-2                
##  [11] statmod_1.4.35            future_1.21.0            
##  [13] miniUI_0.1.1.1            withr_2.4.0              
##  [15] colorspace_2.0-0          highr_0.8                
##  [17] knitr_1.30                ROCR_1.0-11              
##  [19] tensor_1.5                listenv_0.8.0            
##  [21] labeling_0.4.2            GenomeInfoDbData_1.2.4   
##  [23] harmony_1.0               polyclip_1.10-0          
##  [25] farver_2.0.3              parallelly_1.23.0        
##  [27] vctrs_0.3.6               generics_0.1.0           
##  [29] ipred_0.9-9               xfun_0.20                
##  [31] randomForest_4.6-14       R6_2.5.0                 
##  [33] ggbeeswarm_0.6.0          rsvd_1.0.3               
##  [35] locfit_1.5-9.4            bitops_1.0-6             
##  [37] spatstat.utils_1.20-2     DelayedArray_0.16.0      
##  [39] assertthat_0.2.1          promises_1.1.1           
##  [41] scales_1.1.1              nnet_7.3-14              
##  [43] beeswarm_0.2.3            gtable_0.3.0             
##  [45] beachmat_2.6.0            globals_0.14.0           
##  [47] goftest_1.2-2             timeDate_3043.102        
##  [49] rlang_0.4.10              splines_4.0.3            
##  [51] lazyeval_0.2.2            ModelMetrics_1.2.2.2     
##  [53] yaml_2.2.1                reshape2_1.4.4           
##  [55] abind_1.4-5               httpuv_1.5.5             
##  [57] tools_4.0.3               lava_1.6.8.1             
##  [59] ellipsis_0.3.1            RColorBrewer_1.1-2       
##  [61] proxy_0.4-24              ggridges_0.5.3           
##  [63] Rcpp_1.0.6                plyr_1.8.6               
##  [65] sparseMatrixStats_1.2.0   zlibbioc_1.36.0          
##  [67] purrr_0.3.4               RCurl_1.98-1.2           
##  [69] rpart_4.1-15              deldir_0.2-9             
##  [71] pbapply_1.4-3             viridis_0.5.1            
##  [73] zoo_1.8-8                 ggrepel_0.9.1            
##  [75] cluster_2.1.0             magrittr_2.0.1           
##  [77] RSpectra_0.16-0           data.table_1.13.6        
##  [79] scattermore_0.7           lmtest_0.9-38            
##  [81] RANN_2.6.1                fitdistrplus_1.1-3       
##  [83] patchwork_1.1.1           mime_0.9                 
##  [85] evaluate_0.14             xtable_1.8-4             
##  [87] gridExtra_2.3             compiler_4.0.3           
##  [89] tibble_3.0.5              KernSmooth_2.23-18       
##  [91] crayon_1.3.4              htmltools_0.5.1          
##  [93] mgcv_1.8-33               later_1.1.0.1            
##  [95] tidyr_1.1.2               lubridate_1.7.9.2        
##  [97] DBI_1.1.1                 formatR_1.7              
##  [99] MASS_7.3-53               Matrix_1.3-2             
## [101] getopt_1.20.3             cli_2.2.0                
## [103] gower_0.2.2               igraph_1.2.6             
## [105] pkgconfig_2.0.3           plotly_4.9.3             
## [107] scuttle_1.0.0             recipes_0.1.15           
## [109] foreach_1.5.1             vipor_0.4.5              
## [111] dqrng_0.2.1               XVector_0.30.0           
## [113] prodlim_2019.11.13        stringr_1.4.0            
## [115] digest_0.6.27             sctransform_0.3.2        
## [117] RcppAnnoy_0.0.18          spatstat.data_1.7-0      
## [119] fastmatch_1.1-0           rmarkdown_2.6            
## [121] leiden_0.3.6              uwot_0.1.10              
## [123] edgeR_3.32.0              DelayedMatrixStats_1.12.0
## [125] googleVis_0.6.9           kernlab_0.9-29           
## [127] shiny_1.5.0               lifecycle_0.2.0          
## [129] nlme_3.1-151              jsonlite_1.7.2           
## [131] BiocNeighbors_1.8.0       fansi_0.4.2              
## [133] viridisLite_0.3.0         limma_3.46.0             
## [135] pillar_1.4.7              fastmap_1.0.1            
## [137] httr_1.4.2                survival_3.2-7           
## [139] glue_1.4.2                FNN_1.1.3                
## [141] spatstat_1.64-1           png_0.1-7                
## [143] iterators_1.0.13          bluster_1.0.0            
## [145] class_7.3-17              stringi_1.5.3            
## [147] BiocSingular_1.6.0        irlba_2.3.3              
## [149] e1071_1.7-4               future.apply_1.7.0
```



