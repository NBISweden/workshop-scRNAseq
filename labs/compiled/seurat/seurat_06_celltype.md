---
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

## Celltype prediction
***

 Celltype prediction can either be performed on indiviudal cells where each cell gets a predicted celltype label, or on the level of clusters. All methods are based on similarity to other datasets, single cell or sorted bulk RNAseq, or uses known marker genes for each celltype.

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
alldata <- readRDS("data/results//covid_qc_dr_int_cl.rds")
ctrl = alldata[, alldata$orig.ident == "ctrl_13"]

# set active assay to RNA and remove the CCA assay
ctrl@active.assay = "RNA"
ctrl[["CCA"]] = NULL
ctrl
```

```
## An object of class Seurat 
## 18121 features across 1139 samples within 1 assay 
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


Run all steps of the analysis for the ctrl sample as well. Use the clustering from the integration lab with resolution 0.5.


```r
# Set the identity as louvain with resolution 0.5
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
First we will run label transfer using a similar method as in the integration exercise. But, instead of CCA, the default for the 'FindTransferAnchors` function is to use "pcaproject", e.g. the query dataset is projected onto the PCA of the reference dataset. Then, the labels of the reference data are predicted.


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
## maximum number of iterations reached 0.0001222901 -0.0001213054DONE!
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
## |CD4 T cell  | 1620|       50|svmRadial | 0.997| 0.971| 0.974|
## |CD8 T cell  |  945|       50|svmRadial | 0.985| 0.901| 0.978|
## |cDC         |   26|       50|svmRadial | 0.995| 0.547| 1.000|
## |cMono       |  212|       50|svmRadial | 0.995| 0.958| 0.970|
## |ncMono      |   79|       50|svmRadial | 0.998| 0.570| 1.000|
## |NK cell     |  312|       50|svmRadial | 0.999| 0.933| 0.996|
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
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]},{"label":["Plasma cell"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"103","2":"1","3":"2","4":"0","5":"2","6":"0","7":"0","8":"0","9":"0","_rn_":"B cell"},{"1":"0","2":"199","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"CD4 T cell"},{"1":"0","2":"6","3":"227","4":"0","5":"1","6":"0","7":"4","8":"0","9":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"16","5":"8","6":"0","7":"0","8":"0","9":"0","_rn_":"cDC"},{"1":"0","2":"6","3":"2","4":"0","5":"198","6":"4","7":"0","8":"0","9":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"0","5":"8","6":"101","7":"0","8":"0","9":"0","_rn_":"ncMono"},{"1":"0","2":"0","3":"15","4":"0","5":"0","6":"0","7":"150","8":"0","9":"0","_rn_":"NK cell"},{"1":"0","2":"0","3":"0","4":"1","5":"0","6":"0","7":"0","8":"1","9":"0","_rn_":"pDC"},{"1":"0","2":"1","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"2","_rn_":"Plasma cell"},{"1":"0","2":"12","3":"12","4":"0","5":"56","6":"1","7":"0","8":"0","9":"0","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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
# 0.5
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
## 3132 3287 2505 3831 4085 2535 2025 2225 2509 3315
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
## 1:   cMono 0.0000999900 0.000299970 0.9610708 2.081727            0   48
## 2:  ncMono 0.0000999900 0.000299970 0.8399479 1.815004            0   46
## 3:     cDC 0.0000999900 0.000299970 0.8144952 1.754146            0   43
## 4:     pDC 0.0007018951 0.001579264 0.7587093 1.565462            6   21
## 5: NK cell 0.0166666667 0.025000000 0.7579768 1.472597          161   11
## 6:  B cell 0.0133360275 0.024004849 0.7311643 1.468901          131   15
##                                     leadingEdge
## 1:      S100A8,S100A9,LYZ,S100A12,VCAN,FCN1,...
## 2:     CTSS,TYMP,CST3,S100A11,AIF1,SERPINA1,...
## 3:              LYZ,GRN,TYMP,CST3,AIF1,SPI1,...
## 4:         GRN,MS4A6A,CST3,MPEG1,CTSB,TGFBI,...
## 5:       TYROBP,FCER1G,SRGN,CCL3,CD63,MYO1F,...
## 6: NCF1,LY86,MARCH1,HLA-DRB5,POU2F2,PHACTR1,...
## 
## $`1`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:  CD8 T cell 0.0001001201 0.0003504205 0.9473484 2.208290            0   29
## 2:     NK cell 0.0001000901 0.0003504205 0.8349729 1.960430            0   32
## 3:  CD4 T cell 0.0024406479 0.0056948451 0.8721306 1.698773           21    7
## 4: Plasma cell 0.0952047252 0.1332866153 0.5511619 1.280918          950   28
##                                   leadingEdge
## 1:          CD8A,CD3D,GZMH,CCL5,CD3G,CD8B,...
## 2:          CCL5,GZMA,NKG7,GZMM,CST7,CCL4,...
## 3:           CD3D,CD3G,CD3E,IL7R,PIK3IP1,TCF7
## 4: FKBP11,PRDM1,PEBP1,SEC11C,PPIB,SELENOS,...
## 
## $`2`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:      B cell 0.0000999900 0.0004064215 0.8997944 2.004740            0   46
## 2:         cDC 0.0001016054 0.0004064215 0.8786392 1.793060            0   14
## 3:         pDC 0.0009056148 0.0024149728 0.7731768 1.617386            8   18
## 4: Plasma cell 0.0833333333 0.1293493046 0.7147296 1.353972          774    8
## 5:      ncMono 0.0970119784 0.1293493046 0.8397282 1.339850          736    3
##                                             leadingEdge
## 1:      CD79A,TCL1A,LINC00926,MS4A1,CD79B,TNFRSF13C,...
## 2: CD74,HLA-DQB1,HLA-DRA,HLA-DPB1,HLA-DRB1,HLA-DQA1,...
## 3:               CD74,TCF4,BCL11A,IRF8,SPIB,HERPUD1,...
## 4:                PLPP5,ISG20,HERPUD1,MZB1,ITM2C,JCHAIN
## 5:                                  HLA-DPA1,POU2F2,LYN
## 
## $`3`
##       pathway         pval         padj        ES      NES nMoreExtreme size
## 1:    NK cell 0.0000999900 0.0004008819 0.9200147 2.398072            0   46
## 2: CD8 T cell 0.0001002205 0.0004008819 0.9065596 2.259676            0   27
## 3:     ncMono 0.0021726701 0.0057937869 0.8715870 1.728037           18    6
## 4:        pDC 0.0205922938 0.0411845876 0.7149673 1.562477          193   10
## 5: CD4 T cell 0.0464333782 0.0742934051 0.8565954 1.448014          344    3
##                                  leadingEdge
## 1:       FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 2:       FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 3:                        FCGR3A,IFITM2,RHOC
## 4: GZMB,C12orf75,HSP90B1,ALOX5AP,PLAC8,RRBP1
## 5:                                 CD3E,CD3D
## 
## $`4`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:     NK cell 0.0000999900 0.0004016064 0.9404031 2.445011            0   50
## 2:  CD8 T cell 0.0001004016 0.0004016064 0.9169186 2.251773            0   24
## 3:         pDC 0.0017845895 0.0047589055 0.8015204 1.773943           16   11
## 4:      ncMono 0.0116861436 0.0233722871 0.7904335 1.607688          104    7
## 5: Plasma cell 0.0334267414 0.0534827862 0.5614518 1.403215          333   29
##                                   leadingEdge
## 1:         SPON2,PRF1,GNLY,GZMB,CD7,CLIC3,...
## 2:        PRF1,GNLY,GZMB,CTSW,NKG7,FGFBP2,...
## 3:  GZMB,PLAC8,C12orf75,RRBP1,ALOX5AP,HSP90B1
## 4:                    FCGR3A,IFITM2,RHOC,HES4
## 5: CD38,FKBP11,HSP90B1,SDF2L1,PPIB,SLAMF7,...
## 
## $`5`
##    pathway       pval      padj        ES      NES nMoreExtreme size
## 1:  B cell 0.04830273 0.2173623 0.8155993 1.448315          433    6
## 2:  ncMono 0.03160316 0.2173623 0.6353541 1.335436          315   35
## 3:     cDC 0.09320932 0.2796280 0.5948777 1.246390          931   33
##                                    leadingEdge
## 1:                   PDLIM1,HLA-DRB5,STX7,NCF1
## 2:    OAZ1,TIMP1,FKBP1A,CST3,IFITM3,FCER1G,...
## 3: GAPDH,FKBP1A,CST3,HLA-DRB5,FCER1G,COTL1,...
## 
## $`6`
##       pathway         pval         padj        ES      NES nMoreExtreme size
## 1: CD4 T cell 0.0001017501 0.0007122507 0.9291289 1.813191            0   13
## 2: CD8 T cell 0.0022764228 0.0079674797 0.8808142 1.591824           20    7
##                          leadingEdge
## 1: IL7R,LTB,LDHB,RCAN3,MAL,NOSIP,...
## 2:      IL32,CD3D,CD3E,CD2,CD3G,CD8B
## 
## $`7`
##       pathway         pval         padj        ES      NES nMoreExtreme size
## 1: CD4 T cell 0.0001015228 0.0006091371 0.9185889 1.959619            0   14
## 2: CD8 T cell 0.0832825265 0.1665650530 0.7847816 1.345127          682    4
##                             leadingEdge
## 1: IL7R,TCF7,TSHZ2,LTB,PIK3IP1,LEF1,...
## 2:                   CD3E,CD3G,CD3D,CD2
## 
## $`8`
##        pathway         pval         padj        ES      NES nMoreExtreme size
## 1:      B cell 0.0000999900 0.0004584352 0.9150315 1.949425            0   45
## 2:         cDC 0.0001018745 0.0004584352 0.9026093 1.770460            0   14
## 3:         pDC 0.0003044758 0.0009134274 0.8385454 1.656346            2   15
## 4: Plasma cell 0.0018220468 0.0040996052 0.7854880 1.560232           17   16
## 5:      ncMono 0.0078833268 0.0138248848 0.9502938 1.484796           59    3
##                                             leadingEdge
## 1:        CD79A,MS4A1,BANK1,TNFRSF13C,CD74,HLA-DQA1,...
## 2: CD74,HLA-DQA1,HLA-DRA,HLA-DPB1,HLA-DQB1,HLA-DPA1,...
## 3:             CD74,JCHAIN,SPIB,HERPUD1,TCF4,CCDC50,...
## 4:                JCHAIN,HERPUD1,ISG20,ITM2C,MZB1,PEBP1
## 5:                                      HLA-DPA1,POU2F2
## 
## $`9`
##    pathway        pval         padj        ES      NES nMoreExtreme size
## 1:  ncMono 0.000099990 0.0002666933 0.9572949 2.043523            0   49
## 2:   cMono 0.000100010 0.0002666933 0.8934553 1.866721            0   33
## 3:     cDC 0.000100010 0.0002666933 0.8488184 1.783052            0   36
## 4: NK cell 0.007321097 0.0146421943 0.7922434 1.513841           70   11
## 5:  B cell 0.032915994 0.0526655897 0.6970239 1.382703          325   16
## 6:     pDC 0.051746835 0.0689957806 0.6841886 1.350268          510   15
##                                               leadingEdge
## 1:                CDKN1C,LST1,FCGR3A,COTL1,AIF1,MS4A7,...
## 2:               LST1,COTL1,AIF1,SERPINA1,FCER1G,PSAP,...
## 3:                   LST1,COTL1,AIF1,FCER1G,SPI1,CST3,...
## 4:             FCGR3A,RHOC,FCER1G,TYROBP,IFITM2,MYO1F,...
## 5: HLA-DPA1,POU2F2,HLA-DRB5,HLA-DPB1,HLA-DRA,HLA-DRB1,...
## 6:                   CST3,NPC2,PLD4,MPEG1,TGFBI,VAMP8,...
```

Selecing top significant overlap per cluster, we can now rename the clusters according to the predicted labels. OBS! Be aware that if you have some clusters that have non-significant p-values for all the gene sets, the cluster label will not be very reliable. Also, the gene sets you are using may not cover all the celltypes you have in your dataset and hence predictions may just be the most similar celltype.
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
    download.file(url = "http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt",
        destfile = "./data/CellMarker_list/Human_cell_markers.txt")
    download.file(url = "http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt",
        destfile = "./data/CellMarker_list/Mouse_cell_markers.txt")
}
```

Read in the gene lists and do some filtering.


```r
# Load the human marker table
markers <- read.delim("data/CellMarker_list/Human_cell_markers.txt")
markers <- markers[markers$speciesType == "Human", ]
markers <- markers[markers$cancerType == "Normal", ]

# Filter by tissue (to reduce computational time and have tissue-specific
# classification)
sort(unique(markers$tissueType))
```

```
##   [1] "Abdominal adipose tissue"       "Adipose tissue"                
##   [3] "Adrenal gland"                  "Adventitia"                    
##   [5] "Airway epithelium"              "Alveolus"                      
##   [7] "Amniotic fluid"                 "Amniotic membrane"             
##   [9] "Antecubital vein"               "Anterior cruciate ligament"    
##  [11] "Artery"                         "Ascites"                       
##  [13] "Bladder"                        "Blood"                         
##  [15] "Blood vessel"                   "Bone"                          
##  [17] "Bone marrow"                    "Brain"                         
##  [19] "Breast"                         "Bronchoalveolar system"        
##  [21] "Brown adipose tissue"           "Cartilage"                     
##  [23] "Chorionic villus"               "Colon"                         
##  [25] "Colorectum"                     "Cornea"                        
##  [27] "Corneal endothelium"            "Corneal epithelium"            
##  [29] "Corpus luteum"                  "Decidua"                       
##  [31] "Deciduous tooth"                "Dental pulp"                   
##  [33] "Dermis"                         "Dorsolateral prefrontal cortex"
##  [35] "Duodenum"                       "Embryo"                        
##  [37] "Embryoid body"                  "Embryonic brain"               
##  [39] "Embryonic prefrontal cortex"    "Embryonic stem cell"           
##  [41] "Endometrium"                    "Endometrium stroma"            
##  [43] "Epithelium"                     "Esophagus"                     
##  [45] "Eye"                            "Fat pad"                       
##  [47] "Fetal brain"                    "Fetal gonad"                   
##  [49] "Fetal kidney"                   "Fetal liver"                   
##  [51] "Fetal pancreas"                 "Foreskin"                      
##  [53] "Gastric corpus"                 "Gastric epithelium"            
##  [55] "Gastric gland"                  "Gastrointestinal tract"        
##  [57] "Germ"                           "Gingiva"                       
##  [59] "Gonad"                          "Gut"                           
##  [61] "Hair follicle"                  "Heart"                         
##  [63] "Hippocampus"                    "Inferior colliculus"           
##  [65] "Intervertebral disc"            "Intestinal crypt"              
##  [67] "Intestine"                      "Jejunum"                       
##  [69] "Kidney"                         "Lacrimal gland"                
##  [71] "Large intestine"                "Laryngeal squamous epithelium" 
##  [73] "Ligament"                       "Limbal epithelium"             
##  [75] "Liver"                          "Lung"                          
##  [77] "Lymph"                          "Lymph node"                    
##  [79] "Lymphoid tissue"                "Mammary epithelium"            
##  [81] "Meniscus"                       "Midbrain"                      
##  [83] "Molar"                          "Muscle"                        
##  [85] "Myocardium"                     "Myometrium"                    
##  [87] "Nasal concha"                   "Nasal epithelium"              
##  [89] "Nerve"                          "Nucleus pulposus"              
##  [91] "Optic nerve"                    "Oral mucosa"                   
##  [93] "Osteoarthritic cartilage"       "Ovarian cortex"                
##  [95] "Ovarian follicle"               "Ovary"                         
##  [97] "Oviduct"                        "Pancreas"                      
##  [99] "Pancreatic acinar tissue"       "Pancreatic islet"              
## [101] "Periodontal ligament"           "Periosteum"                    
## [103] "Peripheral blood"               "Placenta"                      
## [105] "Plasma"                         "Pluripotent stem cell"         
## [107] "Premolar"                       "Primitive streak"              
## [109] "Prostate"                       "Pyloric gland"                 
## [111] "Rectum"                         "Renal glomerulus"              
## [113] "Retina"                         "Retinal pigment epithelium"    
## [115] "Salivary gland"                 "Scalp"                         
## [117] "Sclerocorneal tissue"           "Seminal plasma"                
## [119] "Serum"                          "Sinonasal mucosa"              
## [121] "Skeletal muscle"                "Skin"                          
## [123] "Small intestinal crypt"         "Small intestine"               
## [125] "Spinal cord"                    "Spleen"                        
## [127] "Splenic red pulp"               "Sputum"                        
## [129] "Stomach"                        "Subcutaneous adipose tissue"   
## [131] "Submandibular gland"            "Sympathetic ganglion"          
## [133] "Synovial fluid"                 "Synovium"                      
## [135] "Tendon"                         "Testis"                        
## [137] "Thymus"                         "Thyroid"                       
## [139] "Tonsil"                         "Tooth"                         
## [141] "Umbilical cord"                 "Umbilical cord blood"          
## [143] "Umbilical vein"                 "Undefined"                     
## [145] "Urine"                          "Uterus"                        
## [147] "Vagina"                         "Venous blood"                  
## [149] "Visceral adipose tissue"        "Vocal fold"                    
## [151] "Whartons jelly"                 "White adipose tissue"
```

```r
grep("blood", unique(markers$tissueType), value = T)
```

```
## [1] "Peripheral blood"     "Umbilical cord blood" "Venous blood"
```

```r
markers <- markers[markers$tissueType %in% c("Blood", "Venous blood", "Serum", "Plasma",
    "Spleen", "Bone marrow", "Lymph node"), ]

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
##                     pathway         pval        padj        ES      NES
## 1:               Neutrophil 0.0002025316 0.003544304 0.8526808 1.704756
## 2:   CD1C+_B dendritic cell 0.0000999900 0.003499650 0.7809054 1.694096
## 3: Mesenchymal stromal cell 0.0023451658 0.020520200 0.8637549 1.602301
##    nMoreExtreme size                             leadingEdge
## 1:            1   15 CD14,CSF3R,JAML,C5AR1,FCGR2A,FCGR1A,...
## 2:            0   49 S100A8,S100A9,LYZ,S100A12,VCAN,FCN1,...
## 3:           21    8              CD14,VIM,ANPEP,CD44,PECAM1
## 
## $`1`
##                     pathway         pval         padj        ES      NES
## 1:                   T cell 0.0001106807 0.0007932573 0.9640244 1.878471
## 2:             Myeloid cell 0.0001023646 0.0007932573 0.8481847 1.842254
## 3: Mesenchymal stromal cell 0.0022722109 0.0080789720 0.8770455 1.661782
##    nMoreExtreme size                        leadingEdge
## 1:            0    7    CD8A,CD3D,CD3G,CD3E,CD2,CD5,...
## 2:            0   14 CD3D,CD3G,CD3E,CD81,PTPRC,CD84,...
## 3:           19    6               CD81,B2M,PTPRC,ITGB1
## 
## $`2`
##    pathway         pval        padj        ES     NES nMoreExtreme size
## 1:  B cell 0.0001033805 0.002791275 0.9041857 1.79463            0   11
##                             leadingEdge
## 1: CD79A,MS4A1,FCER2,PAX5,CD24,CD19,...
## 
## $`3`
##                pathway         pval        padj        ES      NES nMoreExtreme
## 1: Natural killer cell 0.0002176752 0.004026992 0.8926156 1.855310            1
## 2:        Myeloid cell 0.0002053810 0.004026992 0.7931967 1.821039            1
## 3:            Platelet 0.0011805108 0.010919725 0.8310461 1.771412           10
##    size                          leadingEdge
## 1:    8  KLRD1,FCGR3A,CD3E,CD2,NCR1,CD3D,...
## 2:   14 FCGR3A,CD3E,CD81,ZAP70,SPN,IL2RB,...
## 3:    9         CCL5,CD63,SPN,BSG,CD226,CD47
## 
## $`4`
##                        pathway         pval        padj        ES      NES
## 1:      CD1C+_A dendritic cell 0.0002250731 0.002775902 0.9074648 1.848085
## 2:         Natural killer cell 0.0002250731 0.002775902 0.8947885 1.822269
## 3: AXL+SIGLEC6+ dendritic cell 0.0001005328 0.002775902 0.7404960 1.815919
##    nMoreExtreme size                                leadingEdge
## 1:            1    7                    AREG,LPCAT1,ADAM8,NR4A2
## 2:            1    7              KLRD1,FCGR3A,KLRC1,NCAM1,NCR1
## 3:            0   23 PTGDS,CTSW,BHLHE40,TSEN54,CX3CR1,PLAC8,...
## 
## $`5`
## Empty data.table (0 rows and 8 cols): pathway,pval,padj,ES,NES,nMoreExtreme...
## 
## $`6`
##    pathway         pval        padj        ES      NES nMoreExtreme size
## 1:  T cell 0.0002237637 0.005095851 0.9546988 1.678656            1    6
##               leadingEdge
## 1: CD3D,CD3E,CD2,CD5,CD3G
## 
## $`7`
## Empty data.table (0 rows and 8 cols): pathway,pval,padj,ES,NES,nMoreExtreme...
## 
## $`8`
##    pathway        pval       padj        ES      NES nMoreExtreme size
## 1:  B cell 0.001736111 0.05034722 0.8236547 1.612319           16   13
##                                   leadingEdge
## 1: CD79A,MS4A1,POU2AF1,CD24,FCGR2B,POU2F2,...
## 
## $`9`
##         pathway         pval        padj        ES      NES nMoreExtreme size
## 1: Myeloid cell 0.0002006018 0.007422267 0.8070912 1.645310            1   22
## 2:   Neutrophil 0.0031250000 0.038541667 0.8268479 1.565546           29   10
## 3:       B cell 0.0077148756 0.071362599 0.8529374 1.542302           70    7
##                               leadingEdge
## 1: FCGR3A,CSF1R,PECAM1,CD68,ITGAX,SPN,...
## 2:       FCGR3A,PECAM1,ITGAX,C5AR1,FCGR2A
## 3:                LMO2,POU2F2,CD86,FCGR2A
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

In this case we do not have any ground truth, and we cannot say which method performs best. You should keep in mind, that any celltype classification method is just a prediction, and you still need to use your common sense and knowledge of the biological system to judge if the results make sense.

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
##  [1] fgsea_1.24.0       caret_6.0-93       lattice_0.20-45    scPred_1.9.2      
##  [5] rafalib_1.0.0      pheatmap_1.0.12    ggplot2_3.4.0      cowplot_1.1.1     
##  [9] dplyr_1.0.10       SeuratObject_4.1.3 Seurat_4.3.0      
## 
## loaded via a namespace (and not attached):
##   [1] fastmatch_1.1-3        plyr_1.8.8             igraph_1.3.5          
##   [4] lazyeval_0.2.2         sp_1.6-0               splines_4.2.2         
##   [7] BiocParallel_1.32.5    listenv_0.9.0          scattermore_0.8       
##  [10] digest_0.6.31          foreach_1.5.2          htmltools_0.5.4       
##  [13] fansi_1.0.4            magrittr_2.0.3         tensor_1.5            
##  [16] cluster_2.1.4          ROCR_1.0-11            limma_3.54.0          
##  [19] recipes_1.0.4          globals_0.16.2         gower_1.0.1           
##  [22] matrixStats_0.63.0     hardhat_1.2.0          timechange_0.2.0      
##  [25] spatstat.sparse_3.0-0  colorspace_2.1-0       ggrepel_0.9.2         
##  [28] xfun_0.36              crayon_1.5.2           jsonlite_1.8.4        
##  [31] progressr_0.13.0       spatstat.data_3.0-0    survival_3.5-0        
##  [34] zoo_1.8-11             iterators_1.0.14       glue_1.6.2            
##  [37] polyclip_1.10-4        gtable_0.3.1           ipred_0.9-13          
##  [40] leiden_0.4.3           kernlab_0.9-31         future.apply_1.10.0   
##  [43] abind_1.4-5            scales_1.2.1           DBI_1.1.3             
##  [46] spatstat.random_3.1-3  miniUI_0.1.1.1         Rcpp_1.0.10           
##  [49] viridisLite_0.4.1      xtable_1.8-4           reticulate_1.26       
##  [52] stats4_4.2.2           lava_1.7.1             prodlim_2019.11.13    
##  [55] htmlwidgets_1.6.1      httr_1.4.4             RColorBrewer_1.1-3    
##  [58] ellipsis_0.3.2         ica_1.0-3              farver_2.1.1          
##  [61] pkgconfig_2.0.3        nnet_7.3-18            sass_0.4.5            
##  [64] uwot_0.1.14            deldir_1.0-6           utf8_1.2.2            
##  [67] labeling_0.4.2         tidyselect_1.2.0       rlang_1.0.6           
##  [70] reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
##  [73] tools_4.2.2            cachem_1.0.6           cli_3.6.0             
##  [76] generics_0.1.3         ggridges_0.5.4         evaluate_0.20         
##  [79] stringr_1.5.0          fastmap_1.1.0          yaml_2.3.7            
##  [82] goftest_1.2-3          ModelMetrics_1.2.2.2   knitr_1.42            
##  [85] fitdistrplus_1.1-8     purrr_1.0.1            RANN_2.6.1            
##  [88] pbapply_1.7-0          future_1.30.0          nlme_3.1-161          
##  [91] mime_0.12              formatR_1.14           compiler_4.2.2        
##  [94] rstudioapi_0.14        beeswarm_0.4.0         plotly_4.10.1         
##  [97] png_0.1-8              spatstat.utils_3.0-1   tibble_3.1.8          
## [100] bslib_0.4.2            stringi_1.7.12         highr_0.10            
## [103] Matrix_1.5-3           vctrs_0.5.2            pillar_1.8.1          
## [106] lifecycle_1.0.3        spatstat.geom_3.0-5    lmtest_0.9-40         
## [109] jquerylib_0.1.4        RcppAnnoy_0.0.20       data.table_1.14.6     
## [112] irlba_2.3.5.1          httpuv_1.6.8           patchwork_1.1.2       
## [115] R6_2.5.1               promises_1.2.0.1       KernSmooth_2.23-20    
## [118] gridExtra_2.3          vipor_0.4.5            parallelly_1.34.0     
## [121] codetools_0.2-18       MASS_7.3-58.2          assertthat_0.2.1      
## [124] withr_2.5.0            sctransform_0.3.5      harmony_0.1.1         
## [127] parallel_4.2.2         grid_4.2.2             rpart_4.1.19          
## [130] timeDate_4022.108      tidyr_1.3.0            class_7.3-21          
## [133] rmarkdown_2.20         Rtsne_0.16             pROC_1.18.0           
## [136] spatstat.explore_3.0-5 shiny_1.7.4            lubridate_1.9.1       
## [139] ggbeeswarm_0.7.1
```

