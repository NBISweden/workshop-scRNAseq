---
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
##     B cell CD4 T cell CD8 T cell        cDC      cMono     ncMono    NK cell 
##        101        186        290         83        210        132        171 
##        pDC 
##          2
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
reference <- reference %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F)

VariableFeatures(ctrl) = hvg.ctrl
ctrl <- ctrl %>%
    ScaleData(verbose = F) %>%
    RunPCA(verbose = F)
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
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]}],"data":[{"1":"101","2":"1","3":"1","4":"0","5":"0","6":"1","7":"1","8":"0","_rn_":"B cell"},{"1":"0","2":"148","3":"52","4":"0","5":"3","6":"0","7":"1","8":"0","_rn_":"CD4 T cell"},{"1":"0","2":"27","3":"181","4":"1","5":"0","6":"1","7":"71","8":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"14","5":"6","6":"1","7":"0","8":"0","_rn_":"cDC"},{"1":"0","2":"0","3":"0","4":"57","5":"179","6":"34","7":"0","8":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"6","5":"10","6":"88","7":"0","8":"0","_rn_":"ncMono"},{"1":"0","2":"1","3":"37","4":"0","5":"0","6":"0","7":"81","8":"0","_rn_":"NK cell"},{"1":"0","2":"0","3":"0","4":"1","5":"0","6":"0","7":"0","8":"1","_rn_":"pDC"},{"1":"0","2":"1","3":"1","4":"0","5":"0","6":"0","7":"0","8":"0","_rn_":"Plasma cell"},{"1":"0","2":"8","3":"18","4":"4","5":"12","6":"7","7":"17","8":"1","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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
# run differential expression in our dataset, using clustering at resolution
# 0.3
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
    x %>%
        as.data.frame() %>%
        filter(p.value < 0.01) %>%
        top_n(-100, p.value) %>%
        top_n(50, logFC) %>%
        rownames()
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
## 1:     B cell 0.0002187705 0.0004922336  0.9673032  2.043419            0   47
## 2: CD4 T cell 0.0002178649 0.0004922336  0.8544088  1.825635            0   50
## 3:        cDC 0.0004370629 0.0007867133  0.9501640  1.697292            1   17
## 4: CD8 T cell 0.0035016587 0.0045021326 -0.8965506 -1.589312           18   17
## 5:      cMono 0.0007365126 0.0011047689 -0.8093741 -1.681008            3   47
## 6:     ncMono 0.0001847746 0.0004922336 -0.9077901 -1.894947            0   49
## 7:    NK cell 0.0001847746 0.0004922336 -0.9104637 -1.900528            0   49
##                                                leadingEdge
## 1:          MS4A1,CD37,TNFRSF13C,CXCR4,BANK1,LINC00926,...
## 2:                   RPS29,RPS6,RPL13,RPL32,RPL3,RPL21,...
## 3: HLA-DRA,HLA-DQB1,HLA-DPB1,HLA-DRB1,HLA-DPA1,HLA-DMA,...
## 4:                        CCL5,IL32,GZMH,CD3D,CD2,LYAR,...
## 5:                S100A6,S100A9,LYZ,S100A8,TYROBP,FCN1,...
## 6:                S100A4,FCER1G,S100A11,AIF1,LST1,PSAP,...
## 7:                     ITGB2,HCST,NKG7,GNLY,MYO1F,CST7,...
## 
## $`2`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:      cMono 0.0001840943 0.0006597757  0.9367387  1.948662            0   47
## 2:     ncMono 0.0058479532 0.0087719298  0.7943688  1.662719           31   49
## 3:        cDC 0.0435303514 0.0559675947 -0.7997560 -1.449083          217   17
## 4: CD8 T cell 0.0009984026 0.0017971246 -0.8957949 -1.623096            4   17
## 5:    NK cell 0.0004415011 0.0009933775 -0.8023961 -1.755464            1   49
## 6:     B cell 0.0002188184 0.0006597757 -0.8610299 -1.870011            0   47
## 7: CD4 T cell 0.0002199252 0.0006597757 -0.9379170 -2.063293            0   50
##                                                 leadingEdge
## 1:                  S100A8,S100A9,LYZ,S100A12,VCAN,RETN,...
## 2:             S100A11,S100A4,AIF1,FCER1G,SERPINA1,MAFB,...
## 3: HLA-DPB1,HLA-DPA1,HLA-DRB1,HLA-DRA,HLA-DQB1,HLA-DRB5,...
## 4:                         IL32,CCL5,GZMH,CD3D,CD2,LYAR,...
## 5:                       GNLY,B2M,GZMA,CTSW,IFITM1,NKG7,...
## 6:                    RPL23A,RPS5,CXCR4,CD52,CD37,RPS23,...
## 7:                   RPL14,RPL3,RPL5,EEF1A1,RPS4X,RPS3A,...
## 
## $`3`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1: CD4 T cell 0.0002223705 0.0005003336  0.9813251  2.119593            0   50
## 2:     B cell 0.0411265087 0.0528769398  0.6685773  1.427125          183   47
## 3:    NK cell 0.0057992026 0.0086988039 -0.7482515 -1.574989           31   49
## 4:        cDC 0.0011025358 0.0019845645 -0.9150291 -1.630054            5   17
## 5:        pDC 0.0001808973 0.0005003336 -0.8117985 -1.696967            0   47
## 6:      cMono 0.0001808973 0.0005003336 -0.9179753 -1.918916            0   47
## 7:     ncMono 0.0001812251 0.0005003336 -0.9461093 -1.991459            0   49
##                                                 leadingEdge
## 1:                    IL7R,LDHB,PIK3IP1,RPL3,RPS12,RPS6,...
## 2:                RPS5,RPL13A,RPL23A,RPL18A,RPS23,CXCR4,...
## 3:                    NKG7,GNLY,ITGB2,MYO1F,CST7,FGFBP2,...
## 4: HLA-DRA,HLA-DRB1,HLA-DPA1,HLA-DPB1,HLA-DQB1,HLA-DRB5,...
## 5:                      PLEK,NPC2,PLAC8,CTSB,PTPRE,IRF8,...
## 6:                 S100A9,S100A8,LYZ,TYROBP,FCN1,S100A6,...
## 7:                 FCER1G,PSAP,AIF1,LST1,S100A11,IFITM3,...
## 
## $`4`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:        cDC 0.0095602294 0.0143403442  0.8200401  1.740751            9   17
## 2:     ncMono 0.0047403350 0.0085326030 -0.7384567 -1.372096           44   49
## 3:      cMono 0.0042189642 0.0085326030 -0.7475294 -1.385681           39   47
## 4:    NK cell 0.0015801117 0.0047403350 -0.7575689 -1.407607           14   49
## 5:     B cell 0.0003164223 0.0014239004 -0.7963543 -1.476186            2   47
## 6: CD4 T cell 0.0001052521 0.0009472687 -0.9046845 -1.683403            0   50
##                                      leadingEdge
## 1: HLA-DRB5,HLA-DRA,FCER1A,CLEC10A,CD1C,ENHO,...
## 2:     S100A4,IFITM2,CEBPB,SAT1,S100A11,MT2A,...
## 3:      JUND,S100A6,NFKBIA,TSPO,TYROBP,APLP2,...
## 4:       ITGB2,IFITM1,JAK1,HCST,ABHD17A,CST7,...
## 5:         CD52,RPL13A,CXCR4,RPS11,FAU,RPS23,...
## 6:       RPL34,RPS29,RPS14,RPL36,RPL38,RPL30,...
## 
## $`5`
##       pathway         pval        padj         ES       NES nMoreExtreme size
## 1:     ncMono 0.0001122586 0.001010328  0.9743693  1.770887            0   49
## 2:        cDC 0.0068510858 0.012331954  0.8980110  1.481227           52   17
## 3:     B cell 0.0312779267 0.046916890 -0.5775129 -1.433103           34   47
## 4: CD8 T cell 0.0004413063 0.001985878 -0.9260888 -1.878344            0   17
## 5:    NK cell 0.0009140768 0.002108716 -0.8299196 -2.087201            0   49
## 6: CD4 T cell 0.0009372071 0.002108716 -0.8602064 -2.169614            0   50
##                                               leadingEdge
## 1:               LST1,AIF1,COTL1,FCGR3A,FCER1G,CDKN1C,...
## 2: HLA-DPA1,HLA-DRA,HLA-DPB1,HLA-DRB1,HLA-DRB5,MTMR14,...
## 3:          CXCR4,RPL13A,MS4A1,TNFRSF13C,BANK1,RPL23A,...
## 4:                       CCL5,IL32,GZMH,CD3D,CD2,LYAR,...
## 5:                    NKG7,GNLY,CST7,GZMA,CTSW,FGFBP2,...
## 6:                   LDHB,RPS3,IL7R,RPL31,MGAT4A,RPL3,...
## 
## $`6`
##        pathway         pval         padj         ES       NES nMoreExtreme size
## 1:     NK cell 0.0001903674 0.0003790272  0.9275865  1.977060            0   49
## 2:  CD4 T cell 0.0001901502 0.0003790272  0.8976429  1.918022            0   50
## 3:  CD8 T cell 0.0001963479 0.0003790272  0.9679912  1.729516            0   17
## 4: Plasma cell 0.0192381685 0.0247347881  0.8107690  1.534098           99   24
## 5:         cDC 0.0091668364 0.0137502546 -0.8720490 -1.581833           44   17
## 6:      ncMono 0.0002105706 0.0003790272 -0.9131970 -1.970559            0   49
## 7:       cMono 0.0002092488 0.0003790272 -0.9288154 -1.995460            0   47
##                                             leadingEdge
## 1:                  NKG7,GNLY,CST7,GZMA,GZMM,FGFBP2,...
## 2:                 IL7R,RPS3,RPS29,RPL3,RPS4X,RPL14,...
## 3:                    CCL5,IL32,GZMH,CD3D,CD8A,LYAR,...
## 4:             RPL36AL,FKBP11,PPIB,PEBP1,SUB1,ISG20,...
## 5: HLA-DRA,HLA-DRB5,HLA-DMA,HLA-DQB1,BASP1,HLA-DRB1,...
## 6:             FCER1G,AIF1,LST1,COTL1,PSAP,SERPINA1,...
## 7:            S100A9,S100A8,LYZ,TYROBP,FCN1,S100A12,...
## 
## $`7`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:      cMono 0.0001220107 0.0005490483  0.9470790  1.776072            0   47
## 2:     ncMono 0.0001214034 0.0005490483  0.9020106  1.698596            0   49
## 3:        cDC 0.0082621083 0.0106227106  0.8940976  1.493084           57   17
## 4: CD8 T cell 0.0003353454 0.0008620690 -0.9242759 -1.813704            0   17
## 5:    NK cell 0.0005665722 0.0008620690 -0.8015313 -1.923065            0   49
## 6:     B cell 0.0005537099 0.0008620690 -0.8742060 -2.079186            0   47
## 7: CD4 T cell 0.0005747126 0.0008620690 -0.9053300 -2.181123            0   50
##                                                leadingEdge
## 1:                  S100A8,S100A9,LYZ,FCN1,TYROBP,VCAN,...
## 2:                AIF1,S100A11,FCER1G,PSAP,FTH1,IFITM3,...
## 3: HLA-DRA,HLA-DRB1,HLA-DRB5,HLA-DMA,HLA-DQB1,HLA-DPA1,...
## 4:                        CCL5,IL32,GZMH,CD3D,CD2,CD8A,...
## 5:                     NKG7,GNLY,CST7,GZMA,CTSW,FGFBP2,...
## 6:                CXCR4,RPS5,RPL23A,CD52,RPL13A,RPL18A,...
## 7:                   RPL3,RPS4X,RPS3,RPS29,RPS27A,IL7R,...
## 
## $`8`
##       pathway         pval         padj         ES       NES nMoreExtreme size
## 1:    NK cell 0.0002553626 0.0007403751  0.9800797  2.116033            0   49
## 2: CD8 T cell 0.0112747354 0.0144960884  0.8788425  1.620313           48   17
## 3:        cDC 0.0037128713 0.0055693069 -0.8893133 -1.568723           20   17
## 4: CD4 T cell 0.0022977187 0.0041358936 -0.7884657 -1.654895           13   50
## 5:     B cell 0.0003290556 0.0007403751 -0.8344794 -1.735751            1   47
## 6:     ncMono 0.0001643115 0.0007403751 -0.8622527 -1.805009            0   49
## 7:      cMono 0.0001645278 0.0007403751 -0.8952788 -1.862216            0   47
##                                                leadingEdge
## 1:                     GNLY,NKG7,FGFBP2,CTSW,PRF1,CST7,...
## 2:                   CCL5,GZMH,IL32,LYAR,LINC01871,CD2,...
## 3: HLA-DRA,HLA-DRB1,HLA-DPA1,HLA-DQB1,HLA-DRB5,HLA-DMA,...
## 4:                 RPS13,RPS28,TMEM123,RPL22,IL7R,TPT1,...
## 5:                  CD37,RPS11,RPL18A,CD52,RPL12,RPS23,...
## 6:                  COTL1,AIF1,FTH1,LST1,SAT1,SERPINA1,...
## 7:                  S100A9,S100A8,LYZ,FCN1,TKT,S100A12,...
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
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]}
# )
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
##                        pathway        pval       padj         ES       NES
## 1:              Pyramidal cell 0.007781363 0.06966172 -0.9655176 -1.474180
## 2: CD4+CD25+ regulatory T cell 0.003036629 0.04391433 -0.9776061 -1.492637
## 3:      Lake et al.Science.Ex4 0.007347400 0.06906556 -0.9478770 -1.504341
##    nMoreExtreme size              leadingEdge
## 1:           40    6                CD3E,NRGN
## 2:           15    6 CD3E,CD3D,PTPRC,CD3G,CD4
## 3:           38    8                    ANXA1
## 
## $`2`
##                            pathway         pval       padj        ES      NES
## 1:                      Neutrophil 0.0001759944 0.01460421 0.9028391 2.009983
## 2:          CD1C+_B dendritic cell 0.0001799532 0.01460421 0.9214147 1.941501
## 3: Monocyte derived dendritic cell 0.0036659878 0.13784114 0.9197785 1.624866
##    nMoreExtreme size                                 leadingEdge
## 1:            0   80 S100A8,S100A9,S100A12,MNDA,CD14,S100A11,...
## 2:            0   53     S100A8,S100A9,LYZ,S100A12,VCAN,RETN,...
## 3:           17   17     S100A8,S100A9,CD14,CST3,ITGAM,SIRPA,...
## 
## $`3`
##              pathway         pval       padj        ES      NES nMoreExtreme
## 1: Naive CD8+ T cell 0.0002260398 0.00708258 0.8565713 2.034078            0
## 2: Naive CD4+ T cell 0.0002207018 0.00708258 0.9304061 1.891055            0
## 3:       CD4+ T cell 0.0004468275 0.01200051 0.9209940 1.775692            1
##    size                            leadingEdge
## 1:   91 LDHB,PIK3IP1,TCF7,NPM1,NOSIP,RCAN3,...
## 2:   34    IL7R,RPS5,TCF7,EEF1B2,NOSIP,MAL,...
## 3:   25        IL7R,LTB,CD3E,CD3D,CD3G,CD2,...
## 
## $`4`
##                  pathway        pval       padj         ES       NES
## 1:         Megakaryocyte 0.005561735 0.09676819  0.7665587  1.770018
## 2:   Natural killer cell 0.008100902 0.10770833 -0.6806129 -1.311980
## 3: CD4+ cytotoxic T cell 0.004199529 0.08772349 -0.6942041 -1.340247
##    nMoreExtreme size                             leadingEdge
## 1:            4   26       PPBP,PF4,GP9,ITGA2B,CD9,GP1BA,...
## 2:           78   84    PTPRC,CD69,CD81,FCGR3A,CD3E,IL7R,...
## 3:           40   86 ZEB2,EFHD2,AHNAK,LITAF,FCGR3A,RUNX3,...
## 
## $`5`
##             pathway         pval       padj        ES      NES nMoreExtreme
## 1: Mesenchymal cell 0.0002186031 0.04109739 0.8268579 1.524749            1
## 2:     Stromal cell 0.0013874436 0.06520985 0.8485252 1.516518           11
## 3:    Hemangioblast 0.0004518753 0.04247628 0.9913871 1.475405            2
##    size                           leadingEdge
## 1:   59   COTL1,S100A4,CTSC,HES4,ZEB2,VIM,...
## 2:   38 PECAM1,TIMP1,VIM,CD44,ICAM3,TIMP2,...
## 3:    7                           PECAM1,CD34
## 
## $`6`
##                  pathway         pval        padj        ES      NES
## 1: CD4+ cytotoxic T cell 0.0001877934 0.007934163 0.8626974 1.976047
## 2:   Natural killer cell 0.0001886792 0.007934163 0.7731065 1.765573
## 3:           CD8+ T cell 0.0001981375 0.007934163 0.9491255 1.737430
##    nMoreExtreme size                       leadingEdge
## 1:            0   86 CCL5,NKG7,GZMH,GNLY,CST7,GZMA,...
## 2:            0   84 NKG7,CD3D,GNLY,CD3E,GZMA,CD3G,...
## 3:            0   19 NKG7,CD3D,CD3E,CD3G,CD8A,GZMK,...
## 
## $`7`
##                   pathway         pval       padj        ES      NES
## 1:             Neutrophil 0.0001136105 0.01129129 0.9064678 1.785437
## 2: CD1C+_B dendritic cell 0.0001201201 0.01129129 0.9175506 1.742603
## 3:           Stromal cell 0.0003784057 0.02250419 0.8917884 1.635669
##    nMoreExtreme size                                 leadingEdge
## 1:            0   80 S100A8,S100A9,S100A11,CD14,S100A12,MNDA,...
## 2:            0   53        S100A8,S100A9,LYZ,FCN1,VCAN,CD14,...
## 3:            2   38          VIM,CD44,ICAM1,TIMP2,BST1,KLF6,...
## 
## $`8`
##                              pathway         pval       padj        ES      NES
## 1:             CD4+ cytotoxic T cell 0.0002655337 0.01213099 0.9422406 2.202078
## 2: Effector CD8+ memory T (Tem) cell 0.0002629503 0.01213099 0.8925514 2.061969
## 3:               Natural killer cell 0.0002632965 0.01213099 0.8419787 1.963832
##    nMoreExtreme size                           leadingEdge
## 1:            0   86   GNLY,NKG7,GZMB,CCL5,FGFBP2,CTSW,...
## 2:            0   79 GNLY,GZMB,FGFBP2,KLRD1,GZMH,SPON2,...
## 3:            0   84   GNLY,NKG7,GZMB,GZMA,CD247,KLRD1,...
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
##  [1] caret_6.0-90                lattice_0.20-45            
##  [3] SeuratObject_4.0.4          Seurat_4.0.6               
##  [5] scmap_1.16.0                scPred_1.9.2               
##  [7] fgsea_1.20.0                msigdbr_7.4.1              
##  [9] enrichR_3.0                 dplyr_1.0.7                
## [11] igraph_1.2.11               pheatmap_1.0.12            
## [13] clustree_0.4.4              ggraph_2.0.5               
## [15] reticulate_1.22             harmony_1.0                
## [17] Rcpp_1.0.8                  umap_0.2.7.0               
## [19] rafalib_1.0.0               scDblFinder_1.8.0          
## [21] DoubletFinder_2.0.3         org.Hs.eg.db_3.14.0        
## [23] AnnotationDbi_1.56.1        cowplot_1.1.1              
## [25] scran_1.22.0                scater_1.22.0              
## [27] ggplot2_3.3.5               scuttle_1.4.0              
## [29] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
## [31] Biobase_2.54.0              GenomicRanges_1.46.0       
## [33] GenomeInfoDb_1.30.0         IRanges_2.28.0             
## [35] S4Vectors_0.32.0            BiocGenerics_0.40.0        
## [37] MatrixGenerics_1.6.0        matrixStats_0.61.0         
## [39] RJSONIO_1.3-1.6             optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] scattermore_0.7           ModelMetrics_1.2.2.2     
##   [3] tidyr_1.1.4               bit64_4.0.5              
##   [5] knitr_1.37                irlba_2.3.5              
##   [7] DelayedArray_0.20.0       data.table_1.14.2        
##   [9] rpart_4.1-15              KEGGREST_1.34.0          
##  [11] RCurl_1.98-1.5            generics_0.1.1           
##  [13] ScaledMatrix_1.2.0        RSQLite_2.2.8            
##  [15] RANN_2.6.1                proxy_0.4-26             
##  [17] future_1.23.0             bit_4.0.4                
##  [19] lubridate_1.8.0           spatstat.data_2.1-2      
##  [21] httpuv_1.6.5              assertthat_0.2.1         
##  [23] viridis_0.6.2             gower_0.2.2              
##  [25] xfun_0.29                 jquerylib_0.1.4          
##  [27] babelgene_21.4            evaluate_0.14            
##  [29] promises_1.2.0.1          fansi_1.0.0              
##  [31] DBI_1.1.2                 htmlwidgets_1.5.4        
##  [33] spatstat.geom_2.3-1       purrr_0.3.4              
##  [35] ellipsis_0.3.2            RSpectra_0.16-0          
##  [37] backports_1.4.1           deldir_1.0-6             
##  [39] sparseMatrixStats_1.6.0   vctrs_0.3.8              
##  [41] remotes_2.4.2             here_1.0.1               
##  [43] ROCR_1.0-11               abind_1.4-5              
##  [45] cachem_1.0.6              withr_2.4.3              
##  [47] batchelor_1.10.0          ggforce_0.3.3            
##  [49] googleVis_0.6.11          checkmate_2.0.0          
##  [51] sctransform_0.3.3         getopt_1.20.3            
##  [53] goftest_1.2-3             cluster_2.1.2            
##  [55] lazyeval_0.2.2            crayon_1.4.2             
##  [57] hdf5r_1.3.5               edgeR_3.36.0             
##  [59] recipes_0.1.17            pkgconfig_2.0.3          
##  [61] labeling_0.4.2            tweenr_1.0.2             
##  [63] nlme_3.1-155              vipor_0.4.5              
##  [65] nnet_7.3-17               rlang_0.4.12             
##  [67] globals_0.14.0            lifecycle_1.0.1          
##  [69] miniUI_0.1.1.1            rsvd_1.0.5               
##  [71] randomForest_4.6-14       rprojroot_2.0.2          
##  [73] polyclip_1.10-0           lmtest_0.9-39            
##  [75] Matrix_1.4-0              zoo_1.8-9                
##  [77] beeswarm_0.4.0            ggridges_0.5.3           
##  [79] png_0.1-7                 viridisLite_0.4.0        
##  [81] rjson_0.2.21              bitops_1.0-7             
##  [83] pROC_1.18.0               KernSmooth_2.23-20       
##  [85] Biostrings_2.62.0         blob_1.2.2               
##  [87] DelayedMatrixStats_1.16.0 stringr_1.4.0            
##  [89] parallelly_1.30.0         beachmat_2.10.0          
##  [91] scales_1.1.1              memoise_2.0.1            
##  [93] magrittr_2.0.1            plyr_1.8.6               
##  [95] ica_1.0-2                 zlibbioc_1.40.0          
##  [97] compiler_4.1.2            dqrng_0.3.0              
##  [99] RColorBrewer_1.1-2        fitdistrplus_1.1-6       
## [101] cli_3.1.0                 XVector_0.34.0           
## [103] listenv_0.8.0             patchwork_1.1.1          
## [105] pbapply_1.5-0             formatR_1.11             
## [107] MASS_7.3-55               mgcv_1.8-38              
## [109] tidyselect_1.1.1          stringi_1.7.6            
## [111] highr_0.9                 yaml_2.2.1               
## [113] BiocSingular_1.10.0       askpass_1.1              
## [115] locfit_1.5-9.4            ggrepel_0.9.1            
## [117] grid_4.1.2                sass_0.4.0               
## [119] fastmatch_1.1-3           tools_4.1.2              
## [121] future.apply_1.8.1        parallel_4.1.2           
## [123] bluster_1.4.0             foreach_1.5.1            
## [125] metapod_1.2.0             gridExtra_2.3            
## [127] prodlim_2019.11.13        farver_2.1.0             
## [129] Rtsne_0.15                digest_0.6.29            
## [131] BiocManager_1.30.16       FNN_1.1.3                
## [133] lava_1.6.10               shiny_1.7.1              
## [135] later_1.2.0               RcppAnnoy_0.0.19         
## [137] httr_1.4.2                kernlab_0.9-29           
## [139] colorspace_2.0-2          tensor_1.5               
## [141] splines_4.1.2             uwot_0.1.11              
## [143] statmod_1.4.36            spatstat.utils_2.3-0     
## [145] graphlayouts_0.8.0        xgboost_1.5.0.1          
## [147] plotly_4.10.0             xtable_1.8-4             
## [149] jsonlite_1.7.2            tidygraph_1.2.0          
## [151] timeDate_3043.102         ipred_0.9-12             
## [153] R6_2.5.1                  pillar_1.6.4             
## [155] htmltools_0.5.2           mime_0.12                
## [157] glue_1.6.0                fastmap_1.1.0            
## [159] BiocParallel_1.28.0       BiocNeighbors_1.12.0     
## [161] class_7.3-20              codetools_0.2-18         
## [163] utf8_1.2.2                bslib_0.3.1              
## [165] spatstat.sparse_2.1-0     tibble_3.1.6             
## [167] ResidualMatrix_1.4.0      curl_4.3.2               
## [169] ggbeeswarm_0.6.0          leiden_0.3.9             
## [171] openssl_1.4.6             survival_3.2-13          
## [173] limma_3.50.0              rmarkdown_2.11           
## [175] munsell_0.5.0             e1071_1.7-9              
## [177] GenomeInfoDbData_1.2.7    iterators_1.0.13         
## [179] reshape2_1.4.4            gtable_0.3.0             
## [181] spatstat.core_2.3-2
```



