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

Here we will use a reference PBMC dataset from the `scPred` package which is already a Seurat object with counts. And we will test classification based on label transfer using the function `TransferData` in the Seurat package and the `scPred` method. Finally we will use gene set enrichment predict celltype based on the DEGs of each cluster. 

# Load and process data
## Covid-19 data
First, lets load required libraries and the saved object from the clustering step. Subset for one patient.


```r
suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
    library(scPred)
})
```


```r
# load the data and select 'ctrl_13` sample
alldata <- readRDS("data/results/covid_qc_dr_int_cl.rds")
ctrl = alldata[, alldata$orig.ident == "ctrl_13"]

# set active assay to RNA and remove the CCA assay
ctrl@active.assay = "RNA"
ctrl[["CCA"]] = NULL
ctrl
```

```
## An object of class Seurat 
## 18121 features across 1129 samples within 1 assay 
## Active assay: RNA (18121 features, 0 variable features)
##  6 dimensional reductions calculated: umap, tsne, harmony, umap_harmony, scanorama, umap_scanorama
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


## Rerun analysis pipeline
Here, we will run all the steps that we did in previous labs in one go using the `magittr` package with the pipe-operator `%>%`.


```r
reference <- reference %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30)
```


```r
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


Run all steps of the analysis for the ctrl sample as well. Use the clustering from the integration lab with resolution 0.3.


```r
# Set the identity as louvain with resolution 0.3
ctrl <- SetIdent(ctrl, value = "CCA_snn_res.0.5")

ctrl <- ctrl %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(verbose = F) %>%
    RunUMAP(dims = 1:30)
```


```r
DimPlot(ctrl, label = TRUE, repel = TRUE) + NoAxes()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


# Seurat label transfer
First we will run label transfer using a similar method as in the integration exercise. But, instad of CCA the default for the 'FindTransferAnchors` function is to use "pcaproject", e.g. the query datset is projected onto the PCA of the reference dataset. Then, the labels of the reference data are predicted.


```r
transfer.anchors <- FindTransferAnchors(reference = reference, query = ctrl, dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = reference$cell_type,
    dims = 1:30)
ctrl <- AddMetaData(object = ctrl, metadata = predictions)
```


```r
DimPlot(ctrl, group.by = "predicted.id", label = T, repel = T) + NoAxes()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

Now plot how many cells of each celltypes can be found in each cluster.


```r
ggplot(ctrl@meta.data, aes(x = CCA_snn_res.0.5, fill = predicted.id)) + geom_bar() +
    theme_classic()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

# scPred
scPred will train a classifier based on all principal components. First, `getFeatureSpace` will create a scPred object stored in the `@misc` slot where it extracts the PCs that best separates the different celltypes. Then `trainModel` will do the actual training for each celltype.


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
## maximum number of iterations reached 0.000116588 -0.0001156614DONE!
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
## |B cell      |  280|       50|svmRadial | 1.000| 0.964| 1.000|
## |CD4 T cell  | 1620|       50|svmRadial | 0.997| 0.971| 0.975|
## |CD8 T cell  |  945|       50|svmRadial | 0.985| 0.902| 0.978|
## |cDC         |   26|       50|svmRadial | 0.995| 0.547| 1.000|
## |cMono       |  212|       50|svmRadial | 0.994| 0.958| 0.970|
## |ncMono      |   79|       50|svmRadial | 0.998| 0.582| 1.000|
## |NK cell     |  312|       50|svmRadial | 0.999| 0.936| 0.996|
## |pDC         |   20|       50|svmRadial | 1.000| 0.700| 1.000|
## |Plasma cell |    6|       50|svmRadial | 1.000| 0.800| 1.000|
```

You can optimize parameters for each dataset by chaning parameters and testing different types of models, see more at: https://powellgenomicslab.github.io/scPred/articles/introduction.html. But for now, we will continue with this model.

 Now, lets predict celltypes on our data, where scPred will align the two datasets with Harmony and then perform classification.


```r
ctrl <- scPredict(ctrl, reference)
```

```
## ●  Matching reference with new dataset...
## 	 ─ 2000 features present in reference loadings
## 	 ─ 1774 features shared between reference and new dataset
## 	 ─ 88.7% of features in the reference are present in new dataset
## ●  Aligning new data to reference...
## ●  Classifying cells...
## DONE!
```


```r
DimPlot(ctrl, group.by = "scpred_prediction", label = T, repel = T) + NoAxes()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

Now plot how many	cells of each celltypes	can be found in	each cluster.


```r
ggplot(ctrl@meta.data, aes(x = CCA_snn_res.0.5, fill = scpred_prediction)) + geom_bar() +
    theme_classic()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

# Compare results

Now we will compare the output of the two methods using the convenient function in scPred `crossTab` that prints the overlap between two metadata slots.


```r
crossTab(ctrl, "predicted.id", "scpred_prediction")
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]},{"label":["Plasma cell"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"102","2":"1","3":"2","4":"0","5":"1","6":"0","7":"0","8":"0","9":"0","_rn_":"B cell"},{"1":"0","2":"199","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"CD4 T cell"},{"1":"0","2":"7","3":"226","4":"0","5":"1","6":"0","7":"3","8":"0","9":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"12","5":"7","6":"0","7":"0","8":"0","9":"0","_rn_":"cDC"},{"1":"0","2":"6","3":"2","4":"0","5":"196","6":"4","7":"0","8":"0","9":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"0","5":"8","6":"101","7":"0","8":"0","9":"0","_rn_":"ncMono"},{"1":"0","2":"0","3":"15","4":"0","5":"0","6":"0","7":"152","8":"0","9":"0","_rn_":"NK cell"},{"1":"0","2":"0","3":"0","4":"0","5":"0","6":"0","7":"0","8":"1","9":"0","_rn_":"pDC"},{"1":"0","2":"1","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"2","_rn_":"Plasma cell"},{"1":"0","2":"11","3":"11","4":"0","5":"56","6":"1","7":"1","8":"0","9":"0","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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
alldata <- SetIdent(alldata, value = "CCA_snn_res.0.5")
DGE_table <- FindAllMarkers(alldata, logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
    min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1,
    assay = "RNA")

# split into a list
DGE_list <- split(DGE_table, DGE_table$cluster)

unlist(lapply(DGE_list, nrow))
```

```
##    0    1    2    3    4    5    6    7    8    9 
## 3167 3306 2532 4082 2090 2486 3756 2221 2493 3357
```


```r
# Compute differential gene expression in reference dataset (that has cell
# annotation)
reference <- SetIdent(reference, value = "cell_type")
reference_markers <- FindAllMarkers(reference, min.pct = 0.1, min.diff.pct = 0.2,
    only.pos = T, max.cells.per.ident = 20, return.thresh = 1)

# Identify the top cell marker genes in reference dataset select top 50 with
# hihgest foldchange among top 100 signifcant genes.
reference_markers <- reference_markers[order(reference_markers$avg_log2FC, decreasing = T),
    ]
reference_markers %>%
    group_by(cluster) %>%
    top_n(-100, p_val) %>%
    top_n(50, avg_log2FC) -> top50_cell_selection

# Transform the markers into a list
ref_list = split(top50_cell_selection$gene, top50_cell_selection$cluster)

unlist(lapply(ref_list, length))
```

```
##  CD8 T cell  CD4 T cell       cMono      B cell     NK cell         pDC 
##          30          15          50          50          50          50 
##      ncMono         cDC Plasma cell 
##          50          50          50
```

Now we can run GSEA for the DEGs from our dataset and check for enrichment of top DEGs in the reference dataset.


```r
suppressPackageStartupMessages(library(fgsea))

# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    gene_rank <- setNames(x$avg_log2FC, x$gene)
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
## $`0`
##    pathway         pval        padj        ES      NES nMoreExtreme size
## 1:   cMono 0.0000999900 0.000299970 0.9600537 2.088612            0   48
## 2:  ncMono 0.0000999900 0.000299970 0.8411398 1.826201            0   46
## 3:     cDC 0.0000999900 0.000299970 0.8207395 1.773425            0   42
## 4:     pDC 0.0005010522 0.001127367 0.7678091 1.584659            4   21
## 5:  B cell 0.0089817338 0.016167121 0.7341094 1.478262           88   16
## 6: NK cell 0.0165876777 0.024881517 0.7579090 1.470174          160   11
##                                     leadingEdge
## 1:      S100A8,S100A9,LYZ,S100A12,VCAN,FCN1,...
## 2:     CTSS,TYMP,CST3,S100A11,AIF1,SERPINA1,...
## 3:              LYZ,GRN,TYMP,CST3,AIF1,SPI1,...
## 4:         GRN,MS4A6A,CST3,MPEG1,CTSB,TGFBI,...
## 5: NCF1,LY86,MARCH1,HLA-DRB5,POU2F2,PHACTR1,...
## 6:       TYROBP,FCER1G,SRGN,CCL3,CD63,MYO1F,...
## 
## $`1`
##       pathway         pval        padj        ES      NES nMoreExtreme size
## 1: CD8 T cell 0.0001002406 0.000350842 0.9452465 2.216529            0   29
## 2:    NK cell 0.0001001904 0.000350842 0.8369610 1.976231            0   32
## 3: CD4 T cell 0.0014391675 0.003358057 0.8778732 1.718662           12    7
##                          leadingEdge
## 1: CD8A,CD3D,CCL5,GZMH,CD3G,CD8B,...
## 2: CCL5,NKG7,GZMM,GZMA,CCL4,CST7,...
## 3:  CD3D,CD3G,CD3E,IL7R,PIK3IP1,TCF7
## 
## $`2`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:      B cell 0.0001000000 0.0004067107 0.8980546 2.008264            0   46
## 2:         cDC 0.0001016777 0.0004067107 0.8768017 1.789299            0   14
## 3:         pDC 0.0009058883 0.0024157021 0.7745777 1.623859            8   18
## 4: Plasma cell 0.0631218573 0.1009949716 0.7337163 1.395562          589    8
## 5:      ncMono 0.0871829105 0.1162438807 0.8469751 1.353869          652    3
##                                             leadingEdge
## 1:      CD79A,TCL1A,LINC00926,MS4A1,CD79B,TNFRSF13C,...
## 2: CD74,HLA-DQB1,HLA-DRA,HLA-DPB1,HLA-DRB1,HLA-DQA1,...
## 3:               CD74,TCF4,BCL11A,IRF8,HERPUD1,SPIB,...
## 4:                PLPP5,ISG20,HERPUD1,MZB1,ITM2C,JCHAIN
## 5:                                  HLA-DPA1,POU2F2,LYN
## 
## $`3`
##        pathway        pval         padj        ES      NES nMoreExtreme size
## 1:     NK cell 0.000099990 0.0004018081 0.9370485 2.439191            0   50
## 2:  CD8 T cell 0.000100452 0.0004018081 0.9114297 2.218342            0   23
## 3:         pDC 0.001153161 0.0030750952 0.7964110 1.759809           10   11
## 4:      ncMono 0.019671038 0.0314736608 0.7699121 1.567094          176    7
## 5: Plasma cell 0.016014413 0.0314736608 0.5808952 1.453494          159   30
##                                    leadingEdge
## 1:          SPON2,PRF1,GNLY,GZMB,CD7,CLIC3,...
## 2:         PRF1,GNLY,GZMB,NKG7,CTSW,FGFBP2,...
## 3:   GZMB,PLAC8,C12orf75,RRBP1,ALOX5AP,HSP90B1
## 4:                     FCGR3A,IFITM2,RHOC,HES4
## 5: CD38,FKBP11,SLAMF7,HSP90B1,SDF2L1,PRDM1,...
## 
## $`4`
##       pathway         pval         padj        ES      NES nMoreExtreme size
## 1: CD4 T cell 0.0001015228 0.0007106599 0.9196564 1.807750            0   14
## 2: CD8 T cell 0.0033717642 0.0118011747 0.8774835 1.582642           30    7
##                          leadingEdge
## 1: IL7R,LTB,LDHB,RCAN3,MAL,NOSIP,...
## 2:      IL32,CD3E,CD3D,CD2,CD3G,CD8B
## 
## $`5`
##    pathway        pval       padj        ES      NES nMoreExtreme size
## 1:  B cell 0.009477124 0.04264706 0.8589320 1.582633           86    7
## 2:  ncMono 0.004900000 0.04264706 0.6654054 1.436328           48   39
## 3:     cDC 0.061818546 0.13909173 0.5965618 1.280875          617   35
## 4:   cMono 0.047195280 0.13909173 0.5863420 1.275923          471   45
##                                    leadingEdge
## 1:                   PDLIM1,HLA-DRB5,NCF1,STX7
## 2:     OAZ1,TIMP1,CST3,IFITM3,FKBP1A,COTL1,...
## 3: GAPDH,CST3,FKBP1A,HLA-DRB5,COTL1,FCER1G,...
## 4:      CST3,COTL1,FCER1G,STXBP2,LYZ,AP1S2,...
## 
## $`6`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:     NK cell 0.0001000200 0.0004014049 0.9230703 2.415986            0   45
## 2:  CD8 T cell 0.0001003512 0.0004014049 0.9103731 2.270813            0   26
## 3:      ncMono 0.0056327900 0.0150207734 0.8707280 1.655455           46    5
## 4:         pDC 0.0208666028 0.0417332056 0.7111345 1.559891          195   10
## 5:  CD4 T cell 0.0411220044 0.0548293391 0.8584391 1.450444          301    3
## 6: Plasma cell 0.0394473368 0.0548293391 0.5470747 1.388203          393   31
##                                   leadingEdge
## 1:        FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 2:        FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 3:                         FCGR3A,IFITM2,RHOC
## 4:  GZMB,C12orf75,HSP90B1,ALOX5AP,RRBP1,PLAC8
## 5:                                  CD3E,CD3D
## 6: PRDM1,FKBP11,HSP90B1,PPIB,SPCS2,SDF2L1,...
## 
## $`7`
##       pathway         pval         padj        ES      NES nMoreExtreme size
## 1: CD4 T cell 0.0001016673 0.0006100041 0.9148165 1.965006            0   14
## 2: CD8 T cell 0.0773082560 0.1213035607 0.7913036 1.362316          632    4
##                             leadingEdge
## 1: IL7R,TCF7,TSHZ2,PIK3IP1,LTB,LEF1,...
## 2:                   CD3E,CD3G,CD3D,CD2
## 
## $`8`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:      B cell 0.0000999900 0.0004561119 0.8999636 1.901873            0   45
## 2:         cDC 0.0001013582 0.0004561119 0.8796120 1.722431            0   14
## 3:         pDC 0.0003023279 0.0009069838 0.8278167 1.649379            2   17
## 4: Plasma cell 0.0083770690 0.0150787243 0.7523949 1.491188           82   16
##                                             leadingEdge
## 1:        CD79A,MS4A1,BANK1,CD74,TNFRSF13C,HLA-DQA1,...
## 2: CD74,HLA-DQA1,HLA-DRA,HLA-DPB1,HLA-DQB1,HLA-DPA1,...
## 3:             CD74,JCHAIN,SPIB,HERPUD1,TCF4,CCDC50,...
## 4:                JCHAIN,HERPUD1,ISG20,ITM2C,MZB1,PEBP1
## 
## $`9`
##    pathway        pval         padj        ES      NES nMoreExtreme size
## 1:  ncMono 0.000099990 0.0002666667 0.9586584 2.039615            0   49
## 2:   cMono 0.000100000 0.0002666667 0.8888565 1.857631            0   35
## 3:     cDC 0.000100000 0.0002666667 0.8407213 1.763844            0   37
## 4: NK cell 0.001729048 0.0034580960 0.8044823 1.563721           16   13
## 5:  B cell 0.031923464 0.0510775428 0.6945353 1.377306          316   16
## 6:     pDC 0.053373615 0.0711648204 0.6746244 1.337821          529   16
##                                               leadingEdge
## 1:                CDKN1C,LST1,FCGR3A,AIF1,COTL1,MS4A7,...
## 2:               LST1,AIF1,COTL1,SERPINA1,FCER1G,PSAP,...
## 3:                   LST1,AIF1,COTL1,FCER1G,CST3,SPI1,...
## 4:             FCGR3A,FCER1G,RHOC,TYROBP,IFITM2,MYO1F,...
## 5: HLA-DPA1,HLA-DRB5,POU2F2,HLA-DRA,HLA-DPB1,HLA-DRB1,...
## 6:                    CST3,NPC2,PLD4,MPEG1,VAMP8,CTSB,...
```

Selecing top significant overlap per cluster, we can now rename the clusters according to the predicted labels. OBS! Be aware that if you have some clusters that have bad p-values for all the gene sets, the cluster label will not be very reliable. Also, the gene sets you are using may not cover all the celltypes you have in your dataset and hence predictions may just be the most similar celltype.
Also, some of the clusters have very similar p-values to multiple celltypes, for instance the ncMono and cMono celltypes are equally good for some clusters.


```r
new.cluster.ids <- unlist(lapply(res, function(x) {
    as.data.frame(x)[1, 1]
}))

alldata$ref_gsea <- new.cluster.ids[as.character(alldata@active.ident)]

cowplot::plot_grid(ncol = 2, DimPlot(alldata, label = T, group.by = "CCA_snn_res.0.5") +
    NoAxes(), DimPlot(alldata, label = T, group.by = "ref_gsea") + NoAxes())
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

Compare to results with the other celltype prediction methods in the ctrl_13 sample.


```r
ctrl$ref_gsea = alldata$ref_gsea[alldata$orig.ident == "ctrl_13"]

cowplot::plot_grid(ncol = 3, DimPlot(ctrl, label = T, group.by = "ref_gsea") + NoAxes() +
    ggtitle("GSEA"), DimPlot(ctrl, label = T, group.by = "predicted.id") + NoAxes() +
    ggtitle("LabelTransfer"), DimPlot(ctrl, label = T, group.by = "scpred_prediction") +
    NoAxes() + ggtitle("scPred"))
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

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
    gene_rank <- setNames(x$avg_log2FC, x$gene)
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
## $`0`
##                   pathway         pval        padj        ES      NES
## 1:             Neutrophil 0.0000999900 0.004256203 0.8254733 1.804418
## 2:             Fibroblast 0.0001038098 0.004256203 0.9047595 1.731978
## 3: CD1C+_B dendritic cell 0.0000999900 0.004256203 0.7927680 1.719238
##    nMoreExtreme size                              leadingEdge
## 1:            0   55 S100A8,S100A9,S100A12,CD14,MNDA,G0S2,...
## 2:            0   10       CD14,VIM,CD36,CKAP4,LRP1,ITGAM,...
## 3:            0   47  S100A8,S100A9,LYZ,S100A12,VCAN,FCN1,...
## 
## $`1`
##                         pathway         pval         padj        ES      NES
## 1:                  CD8+ T cell 0.0001032631 0.0009839762 0.9465174 2.011465
## 2:        CD4+ cytotoxic T cell 0.0001000000 0.0009839762 0.7974635 1.958620
## 3: Megakaryocyte erythroid cell 0.0001001904 0.0009839762 0.8422593 1.958329
##    nMoreExtreme size                        leadingEdge
## 1:            0   12  GZMK,CD8A,CD3D,CD3G,CD8B,CD3E,...
## 2:            0   59 CCL5,GZMH,LYAR,KLRG1,NKG7,GZMM,...
## 3:            0   27  CD8A,CD3D,CD3G,CD3E,KLRG1,CD2,...
## 
## $`2`
##                         pathway         pval       padj        ES      NES
## 1:            Follicular B cell 0.0001028383 0.01306047 0.8930077 1.793405
## 2: Megakaryocyte erythroid cell 0.0010586492 0.02223389 0.8744492 1.689118
##    nMoreExtreme size                         leadingEdge
## 1:            0   12 MS4A1,CD69,FCER2,CD22,CD40,PAX5,...
## 2:            9    9           CD79A,CD83,CD69,FCER2,LY9
## 
## $`3`
##                              pathway         pval        padj        ES
## 1:             CD4+ cytotoxic T cell 0.0000999900 0.003173328 0.8553921
## 2: Effector CD8+ memory T (Tem) cell 0.0000999900 0.003173328 0.8167321
## 3:      Megakaryocyte erythroid cell 0.0001004218 0.003173328 0.8537738
##         NES nMoreExtreme size                            leadingEdge
## 1: 2.273160            0   71    SPON2,PRF1,GNLY,PTGDS,GZMB,NKG7,...
## 2: 2.154911            0   64 SPON2,GNLY,GZMB,FGFBP2,KLRF1,KLRD1,...
## 3: 2.089741            0   24    GZMB,CD7,KLRB1,IL2RB,KLRD1,GZMA,...
## 
## $`4`
##             pathway         pval        padj        ES      NES nMoreExtreme
## 1:      CD4+ T cell 0.0001015125 0.005217925 0.8903852 1.752150            0
## 2:      CD8+ T cell 0.0006146281 0.011740331 0.8572271 1.660181            5
## 3: Activated T cell 0.0005244389 0.011740331 0.8791173 1.643812            4
##    size                          leadingEdge
## 1:   14   IL7R,LTB,CD3E,CD5,CD3D,TNFRSF4,...
## 2:   12      IL7R,CD3E,CD5,CD3D,CD2,CD28,...
## 3:    9 CD3E,CD3D,TNFRSF4,CD28,CD3G,CD27,...
## 
## $`5`
##                   pathway         pval       padj        ES      NES
## 1:          Megakaryocyte 0.0002081599 0.03205662 0.9060829 1.766042
## 2: Circulating fetal cell 0.0073896753 0.37933666 0.8154918 1.589471
##    nMoreExtreme size              leadingEdge
## 1:            1   10 PPBP,PF4,GP9,RASGRP2,CD9
## 2:           70   10        PF4,CD9,ACTB,CD68
## 
## $`6`
##                              pathway       pval        padj        ES      NES
## 1:             CD4+ cytotoxic T cell 0.00009999 0.002760671 0.8609307 2.334597
## 2: Effector CD8+ memory T (Tem) cell 0.00009999 0.002760671 0.8301591 2.232280
## 3:               Natural killer cell 0.00010001 0.002760671 0.8346297 2.186784
##    nMoreExtreme size                           leadingEdge
## 1:            0   74   FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 2:            0   65 FGFBP2,GNLY,GZMB,GZMH,KLRF1,KLRD1,...
## 3:            0   44   GNLY,NKG7,GZMB,KLRF1,KLRD1,GZMA,...
## 
## $`7`
##                  pathway         pval        padj        ES      NES
## 1:     Naive CD4+ T cell 0.0001000801 0.006104884 0.8045268 1.861028
## 2:     Naive CD8+ T cell 0.0000999900 0.006104884 0.7563607 1.856823
## 3: Central memory T cell 0.0003349710 0.013622153 0.9104445 1.704946
##    nMoreExtreme size                             leadingEdge
## 1:            0   29    CCR7,IL7R,TCF7,TSHZ2,TRABD2A,MAL,...
## 2:            0   72 CCR7,TCF7,PIK3IP1,LEF1,TRABD2A,LDHB,...
## 3:            2    6                     CCR7,IL7R,CD28,CD27
## 
## $`8`
##              pathway        pval      padj        ES      NES nMoreExtreme size
## 1: Follicular B cell 0.003512397 0.1735709 0.8160205 1.560233           33   11
##                           leadingEdge
## 1: MS4A1,CD24,CD40,CD22,PAX5,CD69,...
## 
## $`9`
##                         pathway         pval       padj        ES      NES
## 1: Megakaryocyte erythroid cell 0.0006069189 0.04321729 0.8281421 1.622351
## 2:                 Myeloid cell 0.0002002804 0.03244542 0.7875254 1.600725
## 3:                Lymphoid cell 0.0016349888 0.05297364 0.8224772 1.575496
##    nMoreExtreme size                            leadingEdge
## 1:            5   15  FCGR3A,PECAM1,CD68,ITGAX,SPN,CD86,...
## 2:            1   23 FCGR3A,CSF1R,PECAM1,CD68,ITGAX,SPN,...
## 3:           15   12              FCGR3A,CD68,ITGAX,SPN,CD4
```


#CT_GSEA8:


```r
new.cluster.ids <- unlist(lapply(res, function(x) {
    as.data.frame(x)[1, 1]
}))
alldata$cellmarker_gsea <- new.cluster.ids[as.character(alldata@active.ident)]

cowplot::plot_grid(ncol = 2, DimPlot(alldata, label = T, group.by = "ref_gsea") +
    NoAxes(), DimPlot(alldata, label = T, group.by = "cellmarker_gsea") + NoAxes())
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

Do you think that the methods overlap well? Where do you see the most inconsistencies?

In this case we do not have any ground truth, and we cannot say wich method performs best. You should keep in mind, that any celltype classification method is just a prediction, and you still need to use your common sense and knowledge of the biological system to judge if the results make sense.

# Save data
Finally, lets save the data with predictions.


```r
saveRDS(ctrl, "data/results/ctrl13_qc_dr_int_cl_celltype.rds")
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
##  [3] scPred_1.9.2                fgsea_1.20.0               
##  [5] msigdbr_7.4.1               enrichR_3.0                
##  [7] dplyr_1.0.7                 rafalib_1.0.0              
##  [9] pheatmap_1.0.12             clustree_0.4.4             
## [11] ggraph_2.0.5                reticulate_1.22            
## [13] harmony_1.0                 Rcpp_1.0.8                 
## [15] scran_1.22.0                scuttle_1.4.0              
## [17] SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
## [19] Biobase_2.54.0              GenomicRanges_1.46.0       
## [21] GenomeInfoDb_1.30.0         IRanges_2.28.0             
## [23] S4Vectors_0.32.0            BiocGenerics_0.40.0        
## [25] MatrixGenerics_1.6.0        matrixStats_0.61.0         
## [27] ggplot2_3.3.5               cowplot_1.1.1              
## [29] KernSmooth_2.23-20          fields_13.3                
## [31] viridis_0.6.2               viridisLite_0.4.0          
## [33] spam_2.8-0                  DoubletFinder_2.0.3        
## [35] Matrix_1.4-0                SeuratObject_4.0.4         
## [37] Seurat_4.0.6                RJSONIO_1.3-1.6            
## [39] optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] scattermore_0.7           ModelMetrics_1.2.2.2     
##   [3] tidyr_1.1.4               bit64_4.0.5              
##   [5] knitr_1.37                irlba_2.3.5              
##   [7] DelayedArray_0.20.0       data.table_1.14.2        
##   [9] rpart_4.1-15              RCurl_1.98-1.5           
##  [11] generics_0.1.1            ScaledMatrix_1.2.0       
##  [13] callr_3.7.0               RANN_2.6.1               
##  [15] future_1.23.0             bit_4.0.4                
##  [17] lubridate_1.8.0           spatstat.data_2.1-2      
##  [19] httpuv_1.6.5              assertthat_0.2.1         
##  [21] gower_0.2.2               xfun_0.29                
##  [23] jquerylib_0.1.4           babelgene_21.4           
##  [25] evaluate_0.14             promises_1.2.0.1         
##  [27] fansi_1.0.0               igraph_1.2.11            
##  [29] DBI_1.1.2                 htmlwidgets_1.5.4        
##  [31] spatstat.geom_2.3-1       purrr_0.3.4              
##  [33] ellipsis_0.3.2            RSpectra_0.16-0          
##  [35] backports_1.4.1           deldir_1.0-6             
##  [37] sparseMatrixStats_1.6.0   vctrs_0.3.8              
##  [39] remotes_2.4.2             here_1.0.1               
##  [41] ROCR_1.0-11               abind_1.4-5              
##  [43] withr_2.4.3               ggforce_0.3.3            
##  [45] checkmate_2.0.0           sctransform_0.3.3        
##  [47] prettyunits_1.1.1         getopt_1.20.3            
##  [49] goftest_1.2-3             cluster_2.1.2            
##  [51] dotCall64_1.0-1           lazyeval_0.2.2           
##  [53] crayon_1.4.2              hdf5r_1.3.5              
##  [55] edgeR_3.36.0              recipes_0.1.17           
##  [57] pkgconfig_2.0.3           labeling_0.4.2           
##  [59] tweenr_1.0.2              vipor_0.4.5              
##  [61] nlme_3.1-155              nnet_7.3-17              
##  [63] rlang_0.4.12              globals_0.14.0           
##  [65] lifecycle_1.0.1           miniUI_0.1.1.1           
##  [67] rsvd_1.0.5                rprojroot_2.0.2          
##  [69] polyclip_1.10-0           lmtest_0.9-39            
##  [71] zoo_1.8-9                 beeswarm_0.4.0           
##  [73] ggridges_0.5.3            processx_3.5.2           
##  [75] png_0.1-7                 rjson_0.2.21             
##  [77] bitops_1.0-7              pROC_1.18.0              
##  [79] DelayedMatrixStats_1.16.0 stringr_1.4.0            
##  [81] parallelly_1.30.0         beachmat_2.10.0          
##  [83] scales_1.1.1              magrittr_2.0.1           
##  [85] plyr_1.8.6                ica_1.0-2                
##  [87] zlibbioc_1.40.0           compiler_4.1.2           
##  [89] dqrng_0.3.0               RColorBrewer_1.1-2       
##  [91] fitdistrplus_1.1-6        cli_3.1.0                
##  [93] XVector_0.34.0            listenv_0.8.0            
##  [95] patchwork_1.1.1           pbapply_1.5-0            
##  [97] ps_1.6.0                  formatR_1.11             
##  [99] MASS_7.3-55               mgcv_1.8-38              
## [101] tidyselect_1.1.1          stringi_1.7.6            
## [103] highr_0.9                 yaml_2.2.1               
## [105] BiocSingular_1.10.0       locfit_1.5-9.4           
## [107] ggrepel_0.9.1             grid_4.1.2               
## [109] sass_0.4.0                fastmatch_1.1-3          
## [111] tools_4.1.2               future.apply_1.8.1       
## [113] parallel_4.1.2            bluster_1.4.0            
## [115] foreach_1.5.1             metapod_1.2.0            
## [117] gridExtra_2.3             prodlim_2019.11.13       
## [119] farver_2.1.0              Rtsne_0.15               
## [121] digest_0.6.29             lava_1.6.10              
## [123] shiny_1.7.1               later_1.2.0              
## [125] RcppAnnoy_0.0.19          httr_1.4.2               
## [127] kernlab_0.9-29            colorspace_2.0-2         
## [129] tensor_1.5                splines_4.1.2            
## [131] uwot_0.1.11               statmod_1.4.36           
## [133] spatstat.utils_2.3-0      graphlayouts_0.8.0       
## [135] plotly_4.10.0             xtable_1.8-4             
## [137] jsonlite_1.7.2            tidygraph_1.2.0          
## [139] timeDate_3043.102         ipred_0.9-12             
## [141] R6_2.5.1                  pillar_1.6.4             
## [143] htmltools_0.5.2           mime_0.12                
## [145] glue_1.6.0                fastmap_1.1.0            
## [147] BiocParallel_1.28.0       BiocNeighbors_1.12.0     
## [149] class_7.3-20              codetools_0.2-18         
## [151] maps_3.4.0                pkgbuild_1.3.1           
## [153] utf8_1.2.2                bslib_0.3.1              
## [155] spatstat.sparse_2.1-0     tibble_3.1.6             
## [157] ggbeeswarm_0.6.0          curl_4.3.2               
## [159] leiden_0.3.9              survival_3.2-13          
## [161] limma_3.50.0              rmarkdown_2.11           
## [163] munsell_0.5.0             GenomeInfoDbData_1.2.7   
## [165] iterators_1.0.13          reshape2_1.4.4           
## [167] gtable_0.3.0              spatstat.core_2.3-2
```

