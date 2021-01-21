---
author: "Åsa Björklund  &  Paulo Czarnewski"
date: 'January 20, 2021'
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
##          68         107         129          36         200         140 
##     NK cell         pDC Plasma cell  unassigned 
##         254           2           3         199
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
##         103         185         282          35         183         181 
##     NK cell         pDC Plasma cell 
##         166           2           1
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
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]},{"label":["Plasma cell"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"100","2":"0","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"B cell"},{"1":"0","2":"151","3":"47","4":"0","5":"2","6":"1","7":"2","8":"0","9":"0","_rn_":"CD4 T cell"},{"1":"1","2":"25","3":"187","4":"0","5":"0","6":"0","7":"57","8":"0","9":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"12","5":"3","6":"7","7":"0","8":"0","9":"0","_rn_":"cDC"},{"1":"0","2":"0","3":"0","4":"19","5":"155","6":"78","7":"0","8":"0","9":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"2","5":"9","6":"90","7":"0","8":"0","9":"0","_rn_":"ncMono"},{"1":"0","2":"0","3":"26","4":"0","5":"0","6":"0","7":"95","8":"0","9":"0","_rn_":"NK cell"},{"1":"1","2":"1","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"Plasma cell"},{"1":"1","2":"8","3":"22","4":"2","5":"14","6":"5","7":"12","8":"2","9":"1","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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
##    pathway        pval       padj   log2err        ES      NES size
## 1:   cMono 0.005125309 0.04612778 0.4070179 0.8629434 1.525658   46
## 2:  ncMono 0.048223350 0.21700508 0.2311267 0.7710674 1.379413   49
##                                       leadingEdge
## 1:             FCN1,VCAN,APLP2,TSPO,LYZ,CSF3R,...
## 2: SERPINA1,PSAP,TNFRSF1B,LILRB2,FCER1G,PILRA,...
## 
## $`10`
##    pathway       pval      padj   log2err         ES       NES size
## 1:  B cell 0.04434692 0.1995612 0.3217759  0.7801183  1.476355   46
## 2:     cDC 0.03583531 0.1995612 0.3217759 -0.9006493 -1.482889   17
##                                         leadingEdge
## 1:            RPS11,FAU,RPL23A,CD37,CXCR4,CD79B,...
## 2: HLA-DPA1,HLA-DRB1,HLA-DPB1,HLA-DQB1,MTMR14,BASP1
## 
## $`2`
##       pathway         pval         padj   log2err         ES       NES size
## 1:     B cell 7.925438e-05 3.566447e-04 0.5384341 -0.8644575 -1.488027   46
## 2: CD4 T cell 4.503842e-07 4.053458e-06 0.6749629 -0.9115921 -1.575981   49
##                                leadingEdge
## 1: RPS23,CD52,RPL18A,RPL12,RPS11,CXCR4,...
## 2:  RPL34,RPL13,RPS14,EEF1A1,RPS6,RPL5,...
## 
## $`3`
##       pathway        pval       padj   log2err         ES       NES size
## 1: CD8 T cell 0.001329283 0.01196354 0.4550599  0.9776990  1.653051   17
## 2:      cMono 0.016990826 0.07645872 0.3524879 -0.7585525 -1.547230   46
##                            leadingEdge
## 1:  CCL5,KLRG1,CD8A,GZMK,LYAR,CD8B,...
## 2: S100A6,TYROBP,JUND,IRS2,AGTRAP,TSPO
## 
## $`4`
##    pathway         pval        padj   log2err       ES      NES size
## 1: NK cell 0.0004670593 0.004203534 0.4984931 0.956979 1.707951   49
##                                leadingEdge
## 1: NKG7,CTSW,GNLY,ABHD17A,KLRF1,FGFBP2,...
## 
## $`5`
##       pathway         pval         padj   log2err         ES       NES size
## 1:     B cell 1.926751e-05 0.0001734076 0.5756103  0.9572993  1.637943   46
## 2: CD4 T cell 2.126767e-02 0.0638030100 0.3524879  0.8485933  1.466702   49
## 3: CD8 T cell 8.108108e-02 0.1824324324 0.2616635 -0.8155308 -1.424347   17
## 4:        cDC 2.143334e-03 0.0096450024 0.4317077 -0.9358663 -1.634516   17
##                                         leadingEdge
## 1: MS4A1,TNFRSF13C,CD37,LINC00926,BANK1,RALGPS2,...
## 2:           RPS25,RPL30,RPL13,RPL34,RPL9,RPL36,...
## 3:                                       CCL5,PATL2
## 4:       HLA-DPB1,HLA-DRB1,HLA-DRA,HLA-DMA,HLA-DQB1
## 
## $`6`
##       pathway         pval         padj   log2err         ES       NES size
## 1:     ncMono 2.922364e-05 0.0002630128 0.5756103  0.9305173  1.593402   49
## 2: CD8 T cell 4.160544e-02 0.1872244751 0.3217759 -0.8340409 -1.432729   17
##                                   leadingEdge
## 1: SERPINA1,LST1,CDKN1C,COTL1,AIF1,IFITM3,...
## 2:                                  IL32,CCL5
## 
## $`7`
##       pathway         pval         padj   log2err         ES       NES size
## 1: CD4 T cell 7.132176e-09 6.418958e-08 0.7614608  0.9801228  1.719471   49
## 2:     ncMono 2.478267e-04 1.115220e-03 0.4984931 -0.8746644 -1.822240   49
##                                 leadingEdge
## 1:  RPS12,RPS27A,RPS25,IL7R,RPL35A,TPT1,...
## 2: TIMP1,S100A4,COTL1,LYN,TNFRSF1B,TPM3,...
## 
## $`8`
##    pathway       pval      padj   log2err         ES       NES size
## 1:  B cell 0.04171735 0.1877281 0.3217759  0.7278306  1.425956   46
## 2:     cDC 0.02085077 0.1876570 0.3524879 -0.8879238 -1.427659   17
##                               leadingEdge
## 1: RPL23A,RPS5,RPL12,RPS23,RPS11,CD52,...
## 2:         HLA-DRA,HLA-DMA,HLA-DRB5,BASP1
## 
## $`9`
##       pathway         pval        padj   log2err         ES       NES size
## 1:        cDC 0.0004222385 0.003800146 0.4984931  0.9760307  1.561240   17
## 2: CD4 T cell 0.0499726012 0.224876705 0.3217759 -0.7033836 -1.514127   49
##                                              leadingEdge
## 1: FCER1A,HLA-DPA1,HLA-DPB1,CLEC10A,HLA-DRB1,HLA-DRA,...
## 2:                RPS12,RPS25,RPL30,RPL3,RPL14,RPL38,...
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
##                   pathway        pval      padj   log2err        ES      NES
## 1: CD1C+_B dendritic cell 0.004553497 0.1272180 0.4070179 0.8530765 1.523119
## 2:       Cancer stem cell 0.006115925 0.1272180 0.4070179 0.9176287 1.510022
## 3:  Endothelial stem cell 0.003305768 0.1242969 0.4317077 0.9427225 1.506469
##    size                                   leadingEdge
## 1:   53             FCN1,VCAN,CD36,LYZ,CSF3R,CD14,...
## 2:   23          CD44,ITGA5,ALCAM,ALDH1A1,ITGB1,CXCR4
## 3:   19 TNFRSF1B,CSF3R,PECAM1,CD14,ITGA5,TNFRSF1A,...
## 
## $`10`
##                           pathway         pval       padj   log2err         ES
## 1:                 Pyramidal cell 0.0003955116 0.07435617 0.4984931 -0.9933840
## 2: 1-cell stage cell (Blastomere) 0.0096146718 0.34287569 0.3807304 -0.8615699
## 3:                    Acinar cell 0.0079885985 0.34287569 0.3807304 -0.8545811
##          NES size             leadingEdge
## 1: -1.430573    6               CD3E,NRGN
## 2: -1.557199   40                     FTL
## 3: -1.586958   54 IL32,EFHD2,LYZ,SLC25A37
## 
## $`2`
##                pathway         pval        padj   log2err        ES      NES
## 1:            Platelet 9.341777e-04 0.029154425 0.4772708 0.8585723 1.783638
## 2:       Megakaryocyte 4.105803e-05 0.003859455 0.5573322 0.9670228 1.782668
## 3: Embryonic stem cell 2.914851e-03 0.054961611 0.4317077 0.9199849 1.672653
##    size                          leadingEdge
## 1:   45 GP9,ITGA2B,CD9,CD151,GP1BA,TGFB1,...
## 2:   26  PF4,PPBP,GP9,ITGA2B,CD9,RASGRP2,...
## 3:   24     CD9,ITGB1,LEFTY2,ITGA6,CCR4,CD59
## 
## $`3`
##               pathway        pval       padj   log2err        ES      NES size
## 1:      T helper cell 0.001721175 0.06471617 0.4550599 0.8918996 1.653010   56
## 2: Hematopoietic cell 0.002606285 0.08166361 0.4317077 0.9656353 1.623621   19
## 3:        CD8+ T cell 0.003503358 0.08260058 0.4317077 0.9632183 1.619557   19
##                           leadingEdge
## 1: GZMK,IL7R,CD3D,KLRB1,CD3G,CD3E,...
## 2:    CD8A,CD3D,CD3G,CD3E,PTPRC,ITGAM
## 3:  NKG7,CD8A,GZMK,CD8B,IL7R,CD3D,...
## 
## $`4`
##                         pathway         pval        padj   log2err        ES
## 1:          Natural killer cell 8.910583e-05 0.008375948 0.5384341 0.9326028
## 2:        CD4+ cytotoxic T cell 2.779759e-05 0.005225946 0.5756103 0.9285372
## 3: Megakaryocyte erythroid cell 1.136372e-03 0.071212624 0.4550599 0.9126643
##         NES size                           leadingEdge
## 1: 1.748509   84 NKG7,GZMB,GNLY,FCGR3A,KLRF1,CD247,...
## 2: 1.743775   86   NKG7,GZMB,CTSW,GNLY,FCGR3A,CCL5,...
## 3: 1.707496   83 GZMB,FCGR3A,KLRD1,CD7,KLRB1,IL2RB,...
## 
## $`5`
##            pathway        pval      padj   log2err        ES      NES size
## 1: Epithelial cell 0.004871627 0.1705421 0.4070179 0.8918637 1.531581   47
## 2:   M2 macrophage 0.006155124 0.1705421 0.4070179 0.9307351 1.487867   21
## 3:     Plasma cell 0.005855237 0.1705421 0.4070179 0.9479495 1.487859   18
##                             leadingEdge
## 1: CD40,VIM,CD24,KLF6,KRT10,ST6GAL1,...
## 2:                     FCER2,CXCR4,CD68
## 3:      MS4A1,CXCR4,CD19,CD81,LY9,ICAM1
## 
## $`6`
##                         pathway         pval       padj   log2err        ES
## 1:               Secretory cell 0.0009738704 0.06102921 0.4772708 0.9096709
## 2: Megakaryocyte erythroid cell 0.0039041750 0.13617946 0.4317077 0.8334086
## 3:                 Myeloid cell 0.0063891867 0.15056130 0.4070179 0.8241628
##         NES size                            leadingEdge
## 1: 1.509771   39 IFITM3,CTSC,LYPD2,CD74,VMO1,LGALS9,...
## 2: 1.467405   83   PECAM1,FCGR3A,CD68,SPN,ITGAX,CD4,...
## 3: 1.440291   77 PECAM1,FCGR3A,CD68,CSF1R,SPN,IL3RA,...
## 
## $`7`
##                     pathway         pval       padj   log2err        ES
## 1:        Naive CD8+ T cell 0.0005009235 0.04708681 0.4772708 0.8744522
## 2:        Naive CD4+ T cell 0.0020459524 0.10383930 0.4317077 0.9315296
## 3: Morula cell (Blastomere) 0.0033479484 0.10913146 0.4317077 0.8511717
##         NES size                            leadingEdge
## 1: 1.592312   91  PIK3IP1,LDHB,NOSIP,TCF7,RCAN3,MAL,...
## 2: 1.586603   34  IL7R,NOSIP,TCF7,MAL,PRKCA,SERINC5,...
## 3: 1.538274   87 RPS27A,LDHB,RPL34,RPL11,RPL22,RPL7,...
## 
## $`8`
##          pathway         pval      padj   log2err        ES     NES size
## 1: Hemangioblast 0.0006668549 0.1253687 0.4772708 0.9939417 1.47638    7
##    leadingEdge
## 1: PECAM1,CDH5
## 
## $`9`
##                                pathway         pval      padj   log2err
## 1:              CD1C+_A dendritic cell 0.0009346057 0.1757059 0.4772708
## 2: Myeloid conventional dendritic cell 0.0053501317 0.2615627 0.4070179
## 3:            PROM1Low progenitor cell 0.0024028279 0.2258658 0.4317077
##           ES      NES size                                 leadingEdge
## 1: 0.9134740 1.594052   31 FCER1A,CLEC10A,CD1C,CLIC2,ADAM8,TMEM273,...
## 2: 0.9190366 1.453189   16    FCER1A,CLEC10A,CD1C,ITGAX,CLEC9A,CD2,...
## 3: 0.9812644 1.408315    7                                  ALCAM,SGCA
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
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS/LAPACK: /Users/asbj/miniconda3/envs/scRNAseq2021/lib/libopenblasp-r0.3.12.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] caret_6.0-86                lattice_0.20-41            
##  [3] Seurat_3.2.3                scmap_1.12.0               
##  [5] scPred_1.9.0                fgsea_1.16.0               
##  [7] msigdbr_7.2.1               enrichR_2.1                
##  [9] dplyr_1.0.3                 igraph_1.2.6               
## [11] pheatmap_1.0.12             reticulate_1.18            
## [13] harmony_1.0                 Rcpp_1.0.5                 
## [15] venn_1.9                    umap_0.2.7.0               
## [17] rafalib_1.0.0               scDblFinder_1.4.0          
## [19] org.Hs.eg.db_3.12.0         AnnotationDbi_1.52.0       
## [21] cowplot_1.1.1               scran_1.18.0               
## [23] scater_1.18.0               ggplot2_3.3.3              
## [25] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
## [27] Biobase_2.50.0              GenomicRanges_1.42.0       
## [29] GenomeInfoDb_1.26.0         IRanges_2.24.0             
## [31] S4Vectors_0.28.0            BiocGenerics_0.36.0        
## [33] MatrixGenerics_1.2.0        matrixStats_0.57.0         
## [35] RJSONIO_1.3-1.4             optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] tidyselect_1.1.0          RSQLite_2.2.1            
##   [3] htmlwidgets_1.5.3         grid_4.0.3               
##   [5] BiocParallel_1.24.0       Rtsne_0.15               
##   [7] pROC_1.16.2               munsell_0.5.0            
##   [9] codetools_0.2-18          ica_1.0-2                
##  [11] statmod_1.4.35            xgboost_1.3.0.1          
##  [13] future_1.21.0             miniUI_0.1.1.1           
##  [15] withr_2.3.0               batchelor_1.6.0          
##  [17] colorspace_2.0-0          highr_0.8                
##  [19] knitr_1.30                ROCR_1.0-11              
##  [21] tensor_1.5                listenv_0.8.0            
##  [23] labeling_0.4.2            GenomeInfoDbData_1.2.4   
##  [25] polyclip_1.10-0           bit64_4.0.5              
##  [27] farver_2.0.3              parallelly_1.23.0        
##  [29] vctrs_0.3.6               generics_0.1.0           
##  [31] ipred_0.9-9               xfun_0.19                
##  [33] randomForest_4.6-14       R6_2.5.0                 
##  [35] ggbeeswarm_0.6.0          rsvd_1.0.3               
##  [37] locfit_1.5-9.4            hdf5r_1.3.3              
##  [39] bitops_1.0-6              spatstat.utils_1.17-0    
##  [41] DelayedArray_0.16.0       assertthat_0.2.1         
##  [43] promises_1.1.1            scales_1.1.1             
##  [45] nnet_7.3-14               beeswarm_0.2.3           
##  [47] gtable_0.3.0              beachmat_2.6.0           
##  [49] globals_0.14.0            goftest_1.2-2            
##  [51] timeDate_3043.102         rlang_0.4.10             
##  [53] splines_4.0.3             lazyeval_0.2.2           
##  [55] ModelMetrics_1.2.2.2      yaml_2.2.1               
##  [57] reshape2_1.4.4            abind_1.4-5              
##  [59] httpuv_1.5.4              lava_1.6.8.1             
##  [61] tools_4.0.3               ellipsis_0.3.1           
##  [63] RColorBrewer_1.1-2        proxy_0.4-24             
##  [65] ggridges_0.5.2            plyr_1.8.6               
##  [67] sparseMatrixStats_1.2.0   zlibbioc_1.36.0          
##  [69] purrr_0.3.4               RCurl_1.98-1.2           
##  [71] rpart_4.1-15              openssl_1.4.3            
##  [73] deldir_0.2-3              pbapply_1.4-3            
##  [75] viridis_0.5.1             zoo_1.8-8                
##  [77] ggrepel_0.9.0             cluster_2.1.0            
##  [79] magrittr_2.0.1            data.table_1.13.6        
##  [81] RSpectra_0.16-0           scattermore_0.7          
##  [83] ResidualMatrix_1.0.0      lmtest_0.9-38            
##  [85] RANN_2.6.1                fitdistrplus_1.1-3       
##  [87] patchwork_1.1.1           mime_0.9                 
##  [89] evaluate_0.14             xtable_1.8-4             
##  [91] gridExtra_2.3             compiler_4.0.3           
##  [93] tibble_3.0.4              KernSmooth_2.23-18       
##  [95] crayon_1.3.4              htmltools_0.5.0          
##  [97] mgcv_1.8-33               later_1.1.0.1            
##  [99] tidyr_1.1.2               lubridate_1.7.9.2        
## [101] DBI_1.1.0                 formatR_1.7              
## [103] MASS_7.3-53               Matrix_1.3-0             
## [105] getopt_1.20.3             cli_2.2.0                
## [107] gower_0.2.2               pkgconfig_2.0.3          
## [109] recipes_0.1.15            plotly_4.9.2.2           
## [111] scuttle_1.0.0             foreach_1.5.1            
## [113] vipor_0.4.5               admisc_0.11              
## [115] dqrng_0.2.1               XVector_0.30.0           
## [117] prodlim_2019.11.13        stringr_1.4.0            
## [119] digest_0.6.27             sctransform_0.3.2        
## [121] RcppAnnoy_0.0.18          spatstat.data_1.7-0      
## [123] fastmatch_1.1-0           rmarkdown_2.6            
## [125] leiden_0.3.6              uwot_0.1.10              
## [127] edgeR_3.32.0              googleVis_0.6.9          
## [129] DelayedMatrixStats_1.12.0 kernlab_0.9-29           
## [131] curl_4.3                  shiny_1.5.0              
## [133] rjson_0.2.20              lifecycle_0.2.0          
## [135] nlme_3.1-151              jsonlite_1.7.2           
## [137] BiocNeighbors_1.8.0       fansi_0.4.1              
## [139] viridisLite_0.3.0         askpass_1.1              
## [141] limma_3.46.0              pillar_1.4.7             
## [143] fastmap_1.0.1             httr_1.4.2               
## [145] survival_3.2-7            glue_1.4.2               
## [147] FNN_1.1.3                 spatstat_1.64-1          
## [149] iterators_1.0.13          png_0.1-7                
## [151] bluster_1.0.0             bit_4.0.4                
## [153] class_7.3-17              stringi_1.5.3            
## [155] blob_1.2.1                BiocSingular_1.6.0       
## [157] memoise_1.1.0             e1071_1.7-4              
## [159] irlba_2.3.3               future.apply_1.7.0
```



