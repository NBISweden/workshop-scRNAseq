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

Here we will use a reference PBMC dataset from the `scPred` package which is already a Seurat object with counts. And we will test classification based on label transfer using the function `TransferData` in the Seurat package and the `scPred` method. Finally we will use gene set enrichment predict celltype based on the DEGs of each cluster. 

# Load and process data
## Covid-19 data
First, lets load required libraries and the saved object from the clustering step. Subset for one patient.


```r
suppressPackageStartupMessages({
    library(Seurat)
    library(venn)
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
reference <- reference %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% 
    RunPCA(verbose = F) %>% RunUMAP(dims = 1:30)
```


```r
DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()
```

![](seurat_06_celltype_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


Run all steps of the analysis for the ctrl sample as well. Use the clustering from the integration lab with resolution 0.3.


```r
# Set the identity as louvain with resolution 0.3
ctrl <- SetIdent(ctrl, value = "CCA_snn_res.0.5")

ctrl <- ctrl %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = F) %>% 
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
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["B cell"],"name":[1],"type":["int"],"align":["right"]},{"label":["CD4 T cell"],"name":[2],"type":["int"],"align":["right"]},{"label":["CD8 T cell"],"name":[3],"type":["int"],"align":["right"]},{"label":["cDC"],"name":[4],"type":["int"],"align":["right"]},{"label":["cMono"],"name":[5],"type":["int"],"align":["right"]},{"label":["ncMono"],"name":[6],"type":["int"],"align":["right"]},{"label":["NK cell"],"name":[7],"type":["int"],"align":["right"]},{"label":["pDC"],"name":[8],"type":["int"],"align":["right"]},{"label":["Plasma cell"],"name":[9],"type":["int"],"align":["right"]}],"data":[{"1":"102","2":"1","3":"2","4":"0","5":"1","6":"0","7":"0","8":"0","9":"0","_rn_":"B cell"},{"1":"0","2":"198","3":"1","4":"0","5":"0","6":"0","7":"0","8":"0","9":"0","_rn_":"CD4 T cell"},{"1":"0","2":"7","3":"221","4":"0","5":"1","6":"0","7":"8","8":"0","9":"0","_rn_":"CD8 T cell"},{"1":"0","2":"0","3":"0","4":"11","5":"8","6":"0","7":"0","8":"0","9":"0","_rn_":"cDC"},{"1":"0","2":"4","3":"2","4":"0","5":"198","6":"4","7":"0","8":"0","9":"0","_rn_":"cMono"},{"1":"0","2":"0","3":"0","4":"0","5":"10","6":"99","7":"0","8":"0","9":"0","_rn_":"ncMono"},{"1":"0","2":"0","3":"17","4":"0","5":"0","6":"0","7":"150","8":"0","9":"0","_rn_":"NK cell"},{"1":"0","2":"0","3":"0","4":"0","5":"0","6":"0","7":"0","8":"1","9":"0","_rn_":"pDC"},{"1":"0","2":"1","3":"0","4":"0","5":"0","6":"0","7":"0","8":"0","9":"2","_rn_":"Plasma cell"},{"1":"0","2":"11","3":"11","4":"0","5":"54","6":"3","7":"1","8":"0","9":"0","_rn_":"unassigned"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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
alldata <- SetIdent(alldata, value = "CCA_snn_res.0.5")
DGE_table <- FindAllMarkers(alldata, logfc.threshold = 0, test.use = "wilcox", min.pct = 0, 
    min.diff.pct = 0, only.pos = TRUE, max.cells.per.ident = 20, return.thresh = 1, 
    assay = "RNA")
# split into a list
DGE_list <- split(DGE_table, DGE_table$cluster)

unlist(lapply(DGE_list, nrow))
```

```
##    0    1    2    3    4    5    6    7    8    9   10 
## 3846 4338 5326 4564 4195 6209 3570 5852 3945 3842 4177
```


```r
# Compute differential gene expression in reference dataset (that has cell
# annotation)
reference <- SetIdent(reference, value = "cell_type")
reference_markers <- FindAllMarkers(reference, min.pct = 0.1, min.diff.pct = 0.2, 
    only.pos = T, max.cells.per.ident = 20, return.thresh = 1)

# Identify the top cell marker genes in reference dataset select top 50 with
# hihgest foldchange among top 100 signifcant genes.
reference_markers <- reference_markers[order(reference_markers$avg_logFC, decreasing = T), 
    ]
top50_cell_selection <- reference_markers %>% group_by(cluster) %>% top_n(-100, p_val) %>% 
    top_n(50, avg_logFC)

# Transform the markers into a list
ref_list = split(top50_cell_selection$gene, top50_cell_selection$cluster)

unlist(lapply(ref_list, length))
```

```
##  CD8 T cell  CD4 T cell       cMono      B cell     NK cell         pDC 
##          30          14          50          50          50          50 
##      ncMono         cDC Plasma cell 
##          50          50          50
```

Now we can run GSEA for the DEGs from our dataset and check for enrichment of top DEGs in the reference dataset.


```r
suppressPackageStartupMessages(library(fgsea))

# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    gene_rank <- setNames(x$avg_logFC, x$gene)
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
## $`0`
##    pathway         pval         padj   log2err        ES      NES size
## 1:   cMono 1.000000e-10 4.500000e-10        NA 0.9636931 1.929409   48
## 2:  ncMono 1.000000e-10 4.500000e-10        NA 0.8554222 1.706870   46
## 3:     cDC 4.709497e-10 1.412849e-09 0.8012156 0.8370169 1.665622   44
## 4:     pDC 4.356212e-04 9.801478e-04 0.4984931 0.7923584 1.525090   22
## 5:  B cell 4.591618e-03 8.264912e-03 0.4070179 0.7815194 1.453460   15
## 6: NK cell 1.146738e-02 1.720107e-02 0.3807304 0.7876769 1.429233   11
##                                     leadingEdge
## 1:      S100A8,S100A9,LYZ,S100A12,VCAN,FCN1,...
## 2:     CTSS,TYMP,CST3,S100A11,AIF1,SERPINA1,...
## 3:              LYZ,GRN,TYMP,CST3,AIF1,CPVL,...
## 4:         GRN,MS4A6A,CST3,MPEG1,CTSB,TGFBI,...
## 5: NCF1,LY86,MARCH1,HLA-DRB5,POU2F2,PHACTR1,...
## 6:       TYROBP,FCER1G,SRGN,CCL3,CD63,MYO1F,...
## 
## $`1`
##        pathway         pval         padj   log2err        ES      NES size
## 1:      B cell 1.000000e-10 8.000000e-10        NA 0.9102004 1.797025   48
## 2:         cDC 1.376804e-05 5.507216e-05 0.5933255 0.9111434 1.660869   14
## 3:         pDC 1.172339e-04 3.126238e-04 0.5384341 0.8240613 1.549381   22
## 4: Plasma cell 3.597122e-02 5.076142e-02 0.2413400 0.7829680 1.394993   11
## 5:      ncMono 3.807107e-02 5.076142e-02 0.2616635 0.9100346 1.386870    3
##                                             leadingEdge
## 1:      CD79A,TCL1A,LINC00926,MS4A1,CD79B,TNFRSF13C,...
## 2: CD74,HLA-DQB1,HLA-DRA,HLA-DPB1,HLA-DRB1,HLA-DQA1,...
## 3:               CD74,TCF4,BCL11A,IRF8,HERPUD1,SPIB,...
## 4:             PLPP5,ISG20,HERPUD1,MZB1,ITM2C,IGLL5,...
## 5:                                  HLA-DPA1,POU2F2,LYN
## 
## $`2`
##        pathway         pval         padj   log2err        ES      NES size
## 1:  CD8 T cell 0.0000000001 3.500000e-10        NA 0.9570021 2.061123   29
## 2:     NK cell 0.0000000001 3.500000e-10        NA 0.8805921 1.913106   32
## 3:  CD4 T cell 0.0010133621 2.364511e-03 0.4550599 0.9341152 1.632062    5
## 4:         pDC 0.0540254237 9.454449e-02 0.1978220 0.7536490 1.434682    8
## 5: Plasma cell 0.0839160839 1.174825e-01 0.1511488 0.5753573 1.254476   33
##                                   leadingEdge
## 1:          GZMH,CD8A,CD3D,CD3G,CD8B,CCL5,...
## 2:        CCL5,NKG7,GZMA,FGFBP2,CCL4,GZMM,...
## 3:                             CD3G,CD3E,IL7R
## 4:              C12orf75,GZMB,SELENOS,HSP90B1
## 5: FKBP11,PRDM1,PPIB,SEC11C,SELENOS,PEBP1,...
## 
## $`3`
##        pathway         pval         padj   log2err        ES      NES size
## 1:  CD8 T cell 1.000000e-10 6.000000e-10        NA 0.9644968 1.949766   25
## 2:     NK cell 8.975122e-08 2.692537e-07 0.7049757 0.8758142 1.777905   27
## 3:  CD4 T cell 6.021170e-03 1.204234e-02 0.4070179 0.8849819 1.592248    7
## 4: Plasma cell 7.792208e-02 9.350649e-02 0.1574029 0.6186609 1.261402   29
##                                 leadingEdge
## 1:       DUSP2,CCL5,CD3D,LYAR,CD8A,CD3E,...
## 2:       CCL5,KLRB1,GZMM,CMC1,CST7,GZMA,...
## 3:         CD3E,CD3G,IL7R,PIK3IP1,TCF7,LDHB
## 4: FKBP11,PEBP1,PRDM1,SEC11C,PPIB,LMAN1,...
## 
## $`4`
##    pathway         pval         padj   log2err        ES      NES size
## 1:   cMono 5.080943e-06 4.572849e-05 0.6105269 0.7757720 1.569012   45
## 2:  B cell 3.411135e-03 7.675054e-03 0.4317077 0.9283909 1.550573    5
## 3:  ncMono 1.281016e-04 5.764570e-04 0.5188481 0.7619578 1.532447   39
## 4:     cDC 1.113701e-03 3.341103e-03 0.4550599 0.7414548 1.484590   35
## 5:     pDC 2.270230e-02 4.086413e-02 0.3524879 0.6849178 1.362086   29
##                                 leadingEdge
## 1:   CST3,FCER1G,COTL1,LYZ,STXBP2,AP1S2,...
## 2:              PDLIM1,JCHAIN,HLA-DRB5,NCF1
## 3: OAZ1,TIMP1,CST3,FKBP1A,IFITM3,FCER1G,...
## 4: GAPDH,CST3,FCER1G,COTL1,LYZ,HLA-DRB5,...
## 5:      JCHAIN,PTCRA,CST3,TXN,MZB1,CTSB,...
## 
## $`5`
##        pathway         pval         padj   log2err        ES      NES size
## 1:     NK cell 0.0000000001 0.0000000004        NA 0.9474356 2.148061   50
## 2:  CD8 T cell 0.0000000001 0.0000000004        NA 0.9370127 2.002828   23
## 3:         pDC 0.0021256123 0.0042548172 0.4317077 0.8336225 1.634093   10
## 4: Plasma cell 0.0021274086 0.0042548172 0.4317077 0.6796597 1.484603   29
## 5:      ncMono 0.0198845883 0.0318153412 0.3524879 0.7984732 1.481087    7
##                                  leadingEdge
## 1:        SPON2,GNLY,PRF1,GZMB,CD7,CLIC3,...
## 2:       GNLY,PRF1,GZMB,NKG7,FGFBP2,CTSW,...
## 3: GZMB,C12orf75,RRBP1,PLAC8,ALOX5AP,HSP90B1
## 4:  CD38,FKBP11,SLAMF7,SDF2L1,PRDM1,PPIB,...
## 5:                   FCGR3A,IFITM2,RHOC,HES4
## 
## $`6`
##       pathway         pval         padj   log2err        ES      NES size
## 1: CD4 T cell 6.891751e-06 0.0000413505 0.6105269 0.9309759 1.673875   13
## 2: CD8 T cell 3.046951e-03 0.0091408534 0.4317077 0.9122051 1.569613    7
##                          leadingEdge
## 1: IL7R,LTB,LDHB,RCAN3,MAL,NOSIP,...
## 2:      IL32,CD3E,CD3D,CD2,CD3G,CD8B
## 
## $`7`
##        pathway         pval         padj   log2err        ES      NES size
## 1:     NK cell 0.0000000001 0.0000000004        NA 0.9376360 2.152328   46
## 2:  CD8 T cell 0.0000000001 0.0000000004        NA 0.9521736 2.109402   26
## 3:      ncMono 0.0056013963 0.0130494800 0.4070179 0.8767417 1.615374    6
## 4:         pDC 0.0065247400 0.0130494800 0.4070179 0.7813233 1.581389   11
## 5: Plasma cell 0.0121930024 0.0195088038 0.3807304 0.6362332 1.416387   29
##                                   leadingEdge
## 1:        FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 2:        FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 3:                         FCGR3A,IFITM2,RHOC
## 4:  GZMB,C12orf75,HSP90B1,ALOX5AP,RRBP1,PLAC8
## 5: PRDM1,FKBP11,HSP90B1,PPIB,SPCS2,SDF2L1,...
## 
## $`8`
##        pathway         pval         padj   log2err        ES      NES size
## 1:      B cell 1.000000e-10 0.0000000009        NA 0.9155446 1.749016   46
## 2:         cDC 4.757819e-05 0.0002141018 0.5573322 0.9086034 1.635419   14
## 3:         pDC 1.311549e-04 0.0003934648 0.5188481 0.8596737 1.571001   19
## 4: Plasma cell 1.678928e-02 0.0215862198 0.3524879 0.7859744 1.428745   16
##                                             leadingEdge
## 1:        CD79A,MS4A1,BANK1,CD74,TNFRSF13C,HLA-DQA1,...
## 2: CD74,HLA-DQA1,HLA-DRA,HLA-DPB1,HLA-DQB1,HLA-DPA1,...
## 3:             CD74,JCHAIN,SPIB,HERPUD1,TCF4,CCDC50,...
## 4:            JCHAIN,HERPUD1,ISG20,ITM2C,PEBP1,MZB1,...
## 
## $`9`
##       pathway         pval         padj   log2err        ES      NES size
## 1: CD4 T cell 8.334409e-10 5.000645e-09 0.8012156 0.9582136 1.866999   13
## 2: CD8 T cell 6.820950e-02 1.364190e-01 0.1882041 0.8347074 1.377604    4
##                             leadingEdge
## 1: IL7R,TCF7,PIK3IP1,TSHZ2,LTB,LEF1,...
## 2:                   CD3E,CD3G,CD3D,CD2
## 
## $`10`
##    pathway         pval         padj   log2err        ES      NES size
## 1:  ncMono 1.000000e-10 4.000000e-10        NA 0.9647690 1.893443   49
## 2:   cMono 1.000000e-10 4.000000e-10        NA 0.9030410 1.749942   35
## 3:     cDC 2.557038e-08 6.818768e-08 0.7337620 0.8487409 1.651333   38
## 4: NK cell 2.053339e-03 4.106679e-03 0.4317077 0.8397246 1.519935   13
## 5:     pDC 4.116466e-02 6.586345e-02 0.2220560 0.7217437 1.329270   16
## 6:  B cell 9.538153e-02 1.100529e-01 0.1412251 0.6843492 1.260398   16
##                                               leadingEdge
## 1:                CDKN1C,LST1,FCGR3A,AIF1,COTL1,MS4A7,...
## 2:               LST1,AIF1,COTL1,SERPINA1,FCER1G,PSAP,...
## 3:                   LST1,AIF1,COTL1,FCER1G,CST3,SPI1,...
## 4:             FCGR3A,FCER1G,RHOC,TYROBP,IFITM2,MYO1F,...
## 5:                   CST3,NPC2,PLD4,MPEG1,VAMP8,TGFBI,...
## 6: HLA-DPA1,POU2F2,HLA-DRB5,HLA-DRA,HLA-DPB1,HLA-DRB1,...
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
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]} )
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) < 100]
celltype_list <- celltype_list[unlist(lapply(celltype_list, length)) > 5]
```


```r
# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x) {
    gene_rank <- setNames(x$avg_logFC, x$gene)
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
## $`0`
##                   pathway         pval         padj   log2err        ES
## 1:             Neutrophil 1.000000e-10 1.750000e-08        NA 0.8394617
## 2:             Fibroblast 1.732754e-05 7.580799e-04 0.5756103 0.9063763
## 3: CD1C+_B dendritic cell 7.943459e-09 6.950526e-07 0.7477397 0.8059971
##         NES size                              leadingEdge
## 1: 1.688406   60 S100A8,S100A9,S100A12,CD14,MNDA,G0S2,...
## 2: 1.660297   12        CD14,VIM,CD36,CKAP4,LRP1,CD44,...
## 3: 1.610161   50  S100A8,S100A9,LYZ,S100A12,VCAN,FCN1,...
## 
## $`1`
##              pathway         pval        padj   log2err        ES      NES size
## 1: Follicular B cell 1.688543e-05 0.002667899 0.5756103 0.9107084 1.665647   13
## 2:            Neuron 2.859483e-03 0.058778210 0.4317077 0.9001931 1.577511    8
## 3:      Myeloid cell 3.096126e-03 0.058778210 0.4317077 0.8975807 1.572933    8
##                            leadingEdge
## 1: MS4A1,CD69,FCER2,CD22,CD40,PAX5,...
## 2:       PAX5,FOS,CD24,BCL2,PNOC,CD200
## 3:         CD69,FCER2,CD40,CD24,FCGR2B
## 
## $`2`
##                         pathway         pval         padj   log2err        ES
## 1: Megakaryocyte erythroid cell 2.384940e-10 9.360891e-09 0.8140358 0.8835926
## 2:        CD4+ cytotoxic T cell 1.000000e-10 6.306263e-09        NA 0.8381357
## 3:                  CD8+ T cell 3.750808e-09 1.177754e-07 0.7614608 0.9750913
##         NES size                          leadingEdge
## 1: 1.924322   28    CD8A,CD3D,CD3G,CD2,CD3E,KLRG1,...
## 2: 1.921561   70 GZMH,CCL5,NKG7,KLRG1,GZMA,FGFBP2,...
## 3: 1.907747   10     CD8A,CD3D,CD3G,CD8B,CD2,CD3E,...
## 
## $`3`
##                         pathway         pval         padj   log2err        ES
## 1:                  CD8+ T cell 1.323184e-08 1.161774e-06 0.7477397 0.9655347
## 2: Megakaryocyte erythroid cell 3.261017e-08 1.663119e-06 0.7195128 0.8618991
## 3:                  CD4+ T cell 2.901051e-05 6.340868e-04 0.5756103 0.9124337
##         NES size                         leadingEdge
## 1: 1.845359   13   GZMK,CD3D,CD8A,CD3E,CD3G,CD8B,...
## 2: 1.750688   32 CD3D,CD8A,KLRB1,CD3E,KLRG1,CD3G,...
## 3: 1.740878   12     CD3D,CD3E,CD3G,IL7R,CD2,CD7,...
## 
## $`4`
##                   pathway         pval       padj   log2err        ES      NES
## 1:          Megakaryocyte 0.0001335701 0.02190549 0.5188481 0.8568478 1.655419
## 2:               Platelet 0.0020575246 0.16871702 0.4317077 0.7814344 1.536939
## 3: Circulating fetal cell 0.0068398620 0.30026142 0.4070179 0.8143232 1.536824
##    size                         leadingEdge
## 1:   18 PPBP,PF4,GP9,ITGA2B,CD9,RASGRP2,...
## 2:   23 GP9,ITGA2B,CD9,CD151,CD63,ICAM2,...
## 3:   12                   PF4,CD9,ACTB,CD68
## 
## $`5`
##                              pathway  pval    padj log2err        ES      NES
## 1:             CD4+ cytotoxic T cell 1e-10 5.5e-09      NA 0.8883960 2.067957
## 2: Effector CD8+ memory T (Tem) cell 1e-10 5.5e-09      NA 0.8569550 1.982192
## 3:               Natural killer cell 1e-10 5.5e-09      NA 0.8443081 1.912697
##    size                            leadingEdge
## 1:   74    SPON2,GNLY,PRF1,PTGDS,GZMB,NKG7,...
## 2:   63 SPON2,GNLY,GZMB,FGFBP2,KLRF1,KLRD1,...
## 3:   46     GNLY,GZMB,CD7,CD247,NKG7,KLRF1,...
## 
## $`6`
##             pathway         pval         padj   log2err        ES      NES size
## 1:      CD4+ T cell 3.522691e-06 0.0002659632 0.6272567 0.9097284 1.669248   16
## 2: Cytotoxic T cell 5.772752e-06 0.0002905619 0.6105269 0.9832666 1.656278    6
## 3:      CD8+ T cell 7.120350e-04 0.0134396600 0.4772708 0.8712324 1.575798   13
##                           leadingEdge
## 1: IL7R,LTB,CD3E,CD3D,CD5,TNFRSF4,...
## 2:            IL7R,CD3E,CD3D,CD5,CD3G
## 3:    IL7R,CD3E,CD3D,CD5,CD2,CD3G,...
## 
## $`7`
##                              pathway  pval         padj log2err        ES
## 1:             CD4+ cytotoxic T cell 1e-10 5.433333e-09      NA 0.9067957
## 2: Effector CD8+ memory T (Tem) cell 1e-10 5.433333e-09      NA 0.8711152
## 3:               Natural killer cell 1e-10 5.433333e-09      NA 0.8808439
##         NES size                           leadingEdge
## 1: 2.152209   76   FGFBP2,GNLY,NKG7,CST7,GZMB,CTSW,...
## 2: 2.057007   68 FGFBP2,GNLY,GZMB,GZMH,KLRF1,KLRD1,...
## 3: 2.037627   43   GNLY,NKG7,GZMB,KLRF1,KLRD1,GZMA,...
## 
## $`8`
##              pathway        pval      padj   log2err        ES      NES size
## 1: Follicular B cell 0.002072313 0.1074369 0.4317077 0.8693442 1.539359   11
## 2:            Neuron 0.004662605 0.1219958 0.4070179 0.9118174 1.514592    6
## 3:     Memory B cell 0.009400735 0.1219958 0.3807304 0.8909277 1.506049    7
##                           leadingEdge
## 1: MS4A1,CD24,CD40,CD22,PAX5,EBF1,...
## 2:                CD24,PAX5,PNOC,BCL2
## 3:          SPIB,CD19,PAX5,CD27,TLR10
## 
## $`9`
##              pathway         pval         padj   log2err        ES      NES
## 1: Naive CD8+ T cell 1.000000e-10 1.350000e-08        NA 0.8212612 1.830721
## 2: Naive CD4+ T cell 2.669809e-09 1.802121e-07 0.7749390 0.8480762 1.795505
## 3:      Naive T cell 4.985826e-04 1.682716e-02 0.4772708 0.8966874 1.658994
##    size                             leadingEdge
## 1:   80 CCR7,TCF7,PIK3IP1,LEF1,TRABD2A,LDHB,...
## 2:   30    CCR7,IL7R,TCF7,TSHZ2,TRABD2A,MAL,...
## 3:    8           CCR7,IL7R,CD27,CD3E,CD3G,CD3D
## 
## $`10`
##                         pathway         pval       padj   log2err        ES
## 1: Megakaryocyte erythroid cell 0.0004182119 0.03533891 0.4984931 0.8313014
## 2:                Lymphoid cell 0.0027635424 0.07783978 0.4317077 0.8343364
## 3:                 Myeloid cell 0.0002939231 0.03533891 0.4984931 0.7931802
##         NES size                            leadingEdge
## 1: 1.542582   18  FCGR3A,PECAM1,CD68,ITGAX,SPN,CD86,...
## 2: 1.508032   13        FCGR3A,CD68,ITGAX,SPN,CD4,ITGA4
## 3: 1.499494   25 FCGR3A,CSF1R,PECAM1,CD68,ITGAX,SPN,...
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] fgsea_1.16.0    caret_6.0-86    lattice_0.20-41 scPred_1.9.0   
##  [5] rafalib_1.0.0   pheatmap_1.0.12 ggplot2_3.3.3   cowplot_1.1.1  
##  [9] dplyr_1.0.3     venn_1.9        Seurat_3.2.3    RJSONIO_1.3-1.4
## [13] optparse_1.6.6 
## 
## loaded via a namespace (and not attached):
##   [1] fastmatch_1.1-0       plyr_1.8.6            igraph_1.2.6         
##   [4] lazyeval_0.2.2        splines_4.0.3         BiocParallel_1.24.0  
##   [7] listenv_0.8.0         scattermore_0.7       digest_0.6.27        
##  [10] foreach_1.5.1         htmltools_0.5.1       fansi_0.4.2          
##  [13] magrittr_2.0.1        tensor_1.5            cluster_2.1.0        
##  [16] ROCR_1.0-11           limma_3.46.0          recipes_0.1.15       
##  [19] globals_0.14.0        gower_0.2.2           matrixStats_0.57.0   
##  [22] colorspace_2.0-0      ggrepel_0.9.1         xfun_0.20            
##  [25] crayon_1.3.4          jsonlite_1.7.2        spatstat_1.64-1      
##  [28] spatstat.data_1.7-0   survival_3.2-7        zoo_1.8-8            
##  [31] iterators_1.0.13      glue_1.4.2            polyclip_1.10-0      
##  [34] gtable_0.3.0          ipred_0.9-9           leiden_0.3.6         
##  [37] kernlab_0.9-29        future.apply_1.7.0    abind_1.4-5          
##  [40] scales_1.1.1          DBI_1.1.1             miniUI_0.1.1.1       
##  [43] Rcpp_1.0.6            viridisLite_0.3.0     xtable_1.8-4         
##  [46] reticulate_1.18       rsvd_1.0.3            stats4_4.0.3         
##  [49] lava_1.6.8.1          prodlim_2019.11.13    htmlwidgets_1.5.3    
##  [52] httr_1.4.2            getopt_1.20.3         RColorBrewer_1.1-2   
##  [55] ellipsis_0.3.1        ica_1.0-2             pkgconfig_2.0.3      
##  [58] farver_2.0.3          nnet_7.3-14           uwot_0.1.10          
##  [61] deldir_0.2-9          tidyselect_1.1.0      labeling_0.4.2       
##  [64] rlang_0.4.10          reshape2_1.4.4        later_1.1.0.1        
##  [67] munsell_0.5.0         tools_4.0.3           cli_2.2.0            
##  [70] generics_0.1.0        ggridges_0.5.3        evaluate_0.14        
##  [73] stringr_1.4.0         fastmap_1.0.1         yaml_2.2.1           
##  [76] goftest_1.2-2         ModelMetrics_1.2.2.2  knitr_1.30           
##  [79] fitdistrplus_1.1-3    admisc_0.11           purrr_0.3.4          
##  [82] RANN_2.6.1            pbapply_1.4-3         future_1.21.0        
##  [85] nlme_3.1-151          mime_0.9              formatR_1.7          
##  [88] compiler_4.0.3        beeswarm_0.2.3        plotly_4.9.3         
##  [91] png_0.1-7             spatstat.utils_1.20-2 tibble_3.0.5         
##  [94] stringi_1.5.3         highr_0.8             RSpectra_0.16-0      
##  [97] Matrix_1.3-2          vctrs_0.3.6           pillar_1.4.7         
## [100] lifecycle_0.2.0       lmtest_0.9-38         RcppAnnoy_0.0.18     
## [103] data.table_1.13.6     irlba_2.3.3           httpuv_1.5.5         
## [106] patchwork_1.1.1       R6_2.5.0              promises_1.1.1       
## [109] KernSmooth_2.23-18    gridExtra_2.3         vipor_0.4.5          
## [112] parallelly_1.23.0     codetools_0.2-18      MASS_7.3-53          
## [115] assertthat_0.2.1      withr_2.4.0           sctransform_0.3.2    
## [118] harmony_1.0           mgcv_1.8-33           parallel_4.0.3       
## [121] grid_4.0.3            rpart_4.1-15          timeDate_3043.102    
## [124] tidyr_1.1.2           class_7.3-17          rmarkdown_2.6        
## [127] Rtsne_0.15            pROC_1.17.0.1         shiny_1.5.0          
## [130] lubridate_1.7.9.2     ggbeeswarm_0.6.0
```

