---
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


<style>
h1, .h1, h2, .h2, h3, .h3, h4, .h4 { margin-top: 50px }
p.caption {font-size: 0.9em;font-style: italic;color: grey;margin-right: 10%;margin-left: 10%;text-align: justify}
</style>

# Differential gene expression

In this tutorial we will cover about Differetial gene expression, which comprises an extensive range of topics and methods. In single cell, differential expresison can have multiple functionalities such as of identifying marker genes for cell populations, as well as differentially regulated genes across conditions (healthy vs control). We will also exercise on how to account the batch information in your test.

We can first load the data from the clustering session. Moreover, we can already decide which clustering resolution to use. First let's define using the `louvain_2` clustering to identifying differentially expressed genes.  


```r
suppressPackageStartupMessages({
    library(Seurat)
    library(venn)
    library(dplyr)
    library(cowplot)
    library(ggplot2)
    library(pheatmap)
    library(rafalib)
})

alldata <- readRDS("data/3pbmc_qc_dr_int_cl.rds")

# Set the identity as louvain_2 clustering
print(alldata@active.ident[1:10])
```

```
## p3.1k_AAACCCAAGTGGTCAG-1 p3.1k_AAAGGTATCAACTACG-1 p3.1k_AAAGTCCAGCGTGTCC-1 p3.1k_AACACACTCAAGAGTA-1 p3.1k_AACACACTCGACGAGA-1 p3.1k_AACAGGGCAGGAGGTT-1 p3.1k_AACAGGGCAGTGTATC-1 
##                        1                        3                        3                        1                        2                        2                        2 
## p3.1k_AACAGGGTCAGAATAG-1 p3.1k_AACCTGAAGATGGTCG-1 p3.1k_AACGGGATCGTTATCT-1 
##                        2                        3                        3 
## Levels: 1 3 2 4 5
```

```r
alldata <- SetIdent(alldata, value = "kmeans_5")
print(alldata@active.ident[1:10])
```

```
## p3.1k_AAACCCAAGTGGTCAG-1 p3.1k_AAAGGTATCAACTACG-1 p3.1k_AAAGTCCAGCGTGTCC-1 p3.1k_AACACACTCAAGAGTA-1 p3.1k_AACACACTCGACGAGA-1 p3.1k_AACAGGGCAGGAGGTT-1 p3.1k_AACAGGGCAGTGTATC-1 
##                        1                        3                        3                        1                        2                        2                        2 
## p3.1k_AACAGGGTCAGAATAG-1 p3.1k_AACCTGAAGATGGTCG-1 p3.1k_AACGGGATCGTTATCT-1 
##                        2                        3                        3 
## Levels: 1 3 2 4 5
```

## Cell marker genes
***

Let us first compute a ranking for the highly differential genes in each cluster. There are many different tests and parameters to be chosen that can be used to refine your results. When looking for marker genes, we want genes that are positivelly expressed in a cell type and possibly not expressed in the others.


```r
# Compute differentiall expression
markers_genes <- FindAllMarkers(alldata, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
```

We can now select the top 25 up regulated genes for plotting.


```r
top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)
top25
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["p_val"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["avg_logFC"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["pct.1"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["pct.2"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["p_val_adj"],"name":[5],"type":["dbl"],"align":["right"]},{"label":["cluster"],"name":[6],"type":["fctr"],"align":["left"]},{"label":["gene"],"name":[7],"type":["chr"],"align":["left"]}],"data":[{"1":"4.732970e-20","2":"2.8043060","3":"0.982","4":"0.023","5":"7.647059e-16","6":"1","7":"CD79A"},{"1":"4.068382e-19","2":"2.2291954","3":"0.965","4":"0.099","5":"6.573285e-15","6":"1","7":"CD79B"},{"1":"4.671460e-19","2":"2.4351619","3":"0.932","4":"0.015","5":"7.547678e-15","6":"1","7":"MS4A1"},{"1":"2.556867e-17","2":"2.0741276","3":"1.000","4":"0.754","5":"4.131130e-13","6":"1","7":"CD74"},{"1":"5.041266e-17","2":"1.5299104","3":"0.997","4":"0.401","5":"8.145174e-13","6":"1","7":"HLA-DRB1"},{"1":"5.999527e-17","2":"3.6993687","3":"0.889","4":"0.026","5":"9.693435e-13","6":"1","7":"IGHM"},{"1":"1.109482e-16","2":"1.4546238","3":"0.985","4":"0.418","5":"1.792590e-12","6":"1","7":"HLA-DPA1"},{"1":"2.348711e-16","2":"1.7318111","3":"0.947","4":"0.226","5":"3.794812e-12","6":"1","7":"HLA-DQB1"},{"1":"3.009715e-16","2":"1.8284964","3":"0.914","4":"0.155","5":"4.862796e-12","6":"1","7":"HLA-DQA1"},{"1":"4.657565e-16","2":"1.4919689","3":"0.997","4":"0.444","5":"7.525227e-12","6":"1","7":"HLA-DRA"},{"1":"5.993230e-16","2":"1.5078652","3":"0.990","4":"0.410","5":"9.683262e-12","6":"1","7":"HLA-DPB1"},{"1":"2.624133e-15","2":"1.7187452","3":"0.816","4":"0.020","5":"4.239811e-11","6":"1","7":"BANK1"},{"1":"2.159618e-14","2":"1.3907589","3":"0.697","4":"0.010","5":"3.489294e-10","6":"1","7":"TNFRSF13C"},{"1":"2.753447e-14","2":"3.4520026","3":"0.765","4":"0.116","5":"4.448744e-10","6":"1","7":"IGKC"},{"1":"6.381388e-14","2":"2.3762014","3":"0.806","4":"0.006","5":"1.031041e-09","6":"1","7":"IGHD"},{"1":"1.848228e-13","2":"1.7769364","3":"0.811","4":"0.007","5":"2.986182e-09","6":"1","7":"LINC00926"},{"1":"1.278192e-12","2":"1.4510512","3":"0.707","4":"0.010","5":"2.065175e-08","6":"1","7":"CD22"},{"1":"1.459636e-12","2":"2.3363581","3":"0.677","4":"0.005","5":"2.358333e-08","6":"1","7":"TCL1A"},{"1":"4.330466e-12","2":"1.2640090","3":"0.826","4":"0.183","5":"6.996734e-08","6":"1","7":"HLA-DQA2"},{"1":"9.637204e-12","2":"1.2623166","3":"0.732","4":"0.162","5":"1.557083e-07","6":"1","7":"HVCN1"},{"1":"2.805533e-11","2":"1.0229111","3":"0.530","4":"0.030","5":"4.532900e-07","6":"1","7":"P2RX5"},{"1":"3.111000e-11","2":"0.9698941","3":"0.518","4":"0.011","5":"5.026442e-07","6":"1","7":"FAM129C"},{"1":"3.246174e-11","2":"1.1322508","3":"0.621","4":"0.103","5":"5.244843e-07","6":"1","7":"BCL11A"},{"1":"7.616269e-10","2":"1.1219209","3":"0.611","4":"0.014","5":"1.230561e-05","6":"1","7":"SPIB"},{"1":"1.259765e-09","2":"1.3608146","3":"0.619","4":"0.025","5":"2.035402e-05","6":"1","7":"FCER2"},{"1":"3.605454e-14","2":"1.3216389","3":"0.904","4":"0.185","5":"5.825332e-10","6":"3","7":"TRAC"},{"1":"1.183248e-13","2":"1.2584245","3":"0.932","4":"0.562","5":"1.911774e-09","6":"3","7":"LDHB"},{"1":"1.899628e-12","2":"1.1052596","3":"0.895","4":"0.152","5":"3.069229e-08","6":"3","7":"CD3D"},{"1":"2.539620e-10","2":"0.9377455","3":"0.667","4":"0.128","5":"4.103264e-06","6":"3","7":"CD27"},{"1":"2.750335e-10","2":"1.0412135","3":"0.897","4":"0.197","5":"4.443717e-06","6":"3","7":"CD3E"},{"1":"3.659374e-10","2":"1.0597191","3":"0.566","4":"0.142","5":"5.912451e-06","6":"3","7":"TRBC1"},{"1":"8.859758e-10","2":"0.7614416","3":"0.606","4":"0.123","5":"1.431471e-05","6":"3","7":"CD2"},{"1":"1.171037e-09","2":"1.0814736","3":"0.780","4":"0.292","5":"1.892045e-05","6":"3","7":"TRBC2"},{"1":"1.376109e-09","2":"0.9026696","3":"0.860","4":"0.189","5":"2.223379e-05","6":"3","7":"IL32"},{"1":"2.419008e-09","2":"0.4615654","3":"0.493","4":"0.181","5":"3.908392e-05","6":"3","7":"OCIAD2"},{"1":"3.484612e-09","2":"1.3830196","3":"0.817","4":"0.122","5":"5.630087e-05","6":"3","7":"IL7R"},{"1":"7.480748e-09","2":"0.8416513","3":"0.657","4":"0.191","5":"1.208664e-04","6":"3","7":"PIK3IP1"},{"1":"8.154392e-09","2":"1.0528031","3":"0.746","4":"0.115","5":"1.317505e-04","6":"3","7":"CD3G"},{"1":"8.947691e-09","2":"1.1470220","3":"0.741","4":"0.150","5":"1.445678e-04","6":"3","7":"TCF7"},{"1":"1.122471e-08","2":"0.8575523","3":"0.547","4":"0.069","5":"1.813576e-04","6":"3","7":"BCL11B"},{"1":"1.496267e-08","2":"0.8991642","3":"0.468","4":"0.015","5":"2.417518e-04","6":"3","7":"LEF1"},{"1":"2.901385e-08","2":"0.7163203","3":"0.615","4":"0.162","5":"4.687767e-04","6":"3","7":"LAT"},{"1":"3.584147e-08","2":"0.8530620","3":"0.951","4":"0.484","5":"5.790906e-04","6":"3","7":"LTB"},{"1":"1.208246e-07","2":"0.7221613","3":"0.508","4":"0.075","5":"1.952162e-03","6":"3","7":"CD6"},{"1":"1.321135e-07","2":"1.1118266","3":"0.806","4":"0.340","5":"2.134557e-03","6":"3","7":"NOSIP"},{"1":"1.443448e-07","2":"0.5563072","3":"0.590","4":"0.234","5":"2.332179e-03","6":"3","7":"ETS1"},{"1":"1.724671e-07","2":"0.5128230","3":"0.329","4":"0.024","5":"2.786551e-03","6":"3","7":"CAMK4"},{"1":"3.509572e-07","2":"0.7059181","3":"0.773","4":"0.222","5":"5.670415e-03","6":"3","7":"LCK"},{"1":"3.811300e-07","2":"0.8100437","3":"0.440","4":"0.006","5":"6.157918e-03","6":"3","7":"MAL"},{"1":"4.184401e-07","2":"0.7628937","3":"0.531","4":"0.125","5":"6.760737e-03","6":"3","7":"RCAN3"},{"1":"1.274688e-19","2":"4.1716709","3":"0.992","4":"0.115","5":"2.059514e-15","6":"2","7":"LYZ"},{"1":"1.598028e-19","2":"2.7169773","3":"0.988","4":"0.047","5":"2.581934e-15","6":"2","7":"CST3"},{"1":"6.235747e-19","2":"2.4008037","3":"0.986","4":"0.149","5":"1.007510e-14","6":"2","7":"LST1"},{"1":"1.197580e-18","2":"4.3304545","3":"0.966","4":"0.101","5":"1.934929e-14","6":"2","7":"S100A8"},{"1":"1.351885e-18","2":"2.4598313","3":"0.983","4":"0.135","5":"2.184240e-14","6":"2","7":"AIF1"},{"1":"1.690742e-18","2":"1.9276117","3":"0.918","4":"0.006","5":"2.731732e-14","6":"2","7":"CSTA"},{"1":"2.940298e-18","2":"2.6243436","3":"1.000","4":"0.424","5":"4.750640e-14","6":"2","7":"CTSS"},{"1":"3.244740e-18","2":"2.1329251","3":"0.984","4":"0.356","5":"5.242526e-14","6":"2","7":"PSAP"},{"1":"4.957770e-18","2":"4.7801798","3":"0.987","4":"0.136","5":"8.010269e-14","6":"2","7":"S100A9"},{"1":"5.973151e-18","2":"2.7690435","3":"0.966","4":"0.012","5":"9.650820e-14","6":"2","7":"FCN1"},{"1":"5.973151e-18","2":"2.5078751","3":"0.943","4":"0.016","5":"9.650820e-14","6":"2","7":"MNDA"},{"1":"1.376351e-17","2":"1.9567092","3":"0.865","4":"0.028","5":"2.223771e-13","6":"2","7":"KLF4"},{"1":"1.445201e-17","2":"1.7398010","3":"0.892","4":"0.112","5":"2.335011e-13","6":"2","7":"GRN"},{"1":"4.030706e-17","2":"2.3020871","3":"0.990","4":"0.225","5":"6.512411e-13","6":"2","7":"LGALS1"},{"1":"4.315311e-17","2":"3.2006692","3":"0.891","4":"0.014","5":"6.972248e-13","6":"2","7":"S100A12"},{"1":"6.945318e-17","2":"1.9087559","3":"0.881","4":"0.006","5":"1.122155e-12","6":"2","7":"SERPINA1"},{"1":"7.836333e-17","2":"1.6306006","3":"0.976","4":"0.545","5":"1.266116e-12","6":"2","7":"COTL1"},{"1":"1.373695e-16","2":"2.0449917","3":"0.959","4":"0.178","5":"2.219479e-12","6":"2","7":"TYMP"},{"1":"1.410948e-16","2":"2.6111435","3":"0.891","4":"0.008","5":"2.279668e-12","6":"2","7":"VCAN"},{"1":"2.287799e-16","2":"2.1411204","3":"0.882","4":"0.026","5":"3.696397e-12","6":"2","7":"AC020656.1"},{"1":"2.459287e-16","2":"1.6334837","3":"0.964","4":"0.427","5":"3.973470e-12","6":"2","7":"S100A11"},{"1":"4.237595e-16","2":"1.9525141","3":"0.902","4":"0.021","5":"6.846682e-12","6":"2","7":"FGL2"},{"1":"4.932930e-16","2":"1.4636294","3":"0.843","4":"0.185","5":"7.970135e-12","6":"2","7":"APLP2"},{"1":"6.906426e-16","2":"1.4258509","3":"0.855","4":"0.085","5":"1.115871e-11","6":"2","7":"RNF130"},{"1":"7.003520e-16","2":"2.1612191","3":"0.996","4":"0.161","5":"1.131559e-11","6":"2","7":"TYROBP"},{"1":"2.729269e-19","2":"3.4702000","3":"0.995","4":"0.086","5":"4.409680e-15","6":"4","7":"NKG7"},{"1":"2.880475e-17","2":"2.2116900","3":"0.889","4":"0.039","5":"4.653984e-13","6":"4","7":"CST7"},{"1":"4.706080e-17","2":"2.1867445","3":"0.936","4":"0.152","5":"7.603614e-13","6":"4","7":"CTSW"},{"1":"1.875469e-15","2":"2.4400183","3":"0.943","4":"0.059","5":"3.030195e-11","6":"4","7":"GZMA"},{"1":"2.321299e-15","2":"2.1869288","3":"0.856","4":"0.035","5":"3.750523e-11","6":"4","7":"PRF1"},{"1":"2.577118e-13","2":"1.5102666","3":"0.866","4":"0.187","5":"4.163849e-09","6":"4","7":"GZMM"},{"1":"1.103949e-12","2":"1.3237203","3":"0.676","4":"0.027","5":"1.783650e-08","6":"4","7":"MATK"},{"1":"1.899628e-12","2":"2.2245814","3":"0.835","4":"0.111","5":"3.069229e-08","6":"4","7":"KLRB1"},{"1":"9.650835e-12","2":"2.3042919","3":"0.866","4":"0.097","5":"1.559285e-07","6":"4","7":"CCL5"},{"1":"1.093340e-11","2":"1.3846658","3":"0.799","4":"0.245","5":"1.766510e-07","6":"4","7":"CD247"},{"1":"1.515661e-11","2":"3.8578389","3":"0.632","4":"0.037","5":"2.448854e-07","6":"4","7":"GNLY"},{"1":"2.274207e-10","2":"1.0764488","3":"0.715","4":"0.222","5":"3.674436e-06","6":"4","7":"APOBEC3G"},{"1":"4.762658e-10","2":"1.2165461","3":"0.969","4":"0.657","5":"7.695026e-06","6":"4","7":"HCST"},{"1":"5.708744e-10","2":"1.7518652","3":"0.792","4":"0.087","5":"9.223618e-06","6":"4","7":"HOPX"},{"1":"8.497877e-10","2":"1.1236887","3":"0.648","4":"0.094","5":"1.373002e-05","6":"4","7":"SAMD3"},{"1":"1.122157e-09","2":"1.6278557","3":"0.427","4":"0.027","5":"1.813069e-05","6":"4","7":"SPON2"},{"1":"1.122157e-09","2":"1.4865917","3":"0.468","4":"0.007","5":"1.813069e-05","6":"4","7":"KLRF1"},{"1":"2.704358e-09","2":"1.6511816","3":"0.419","4":"0.011","5":"4.369431e-05","6":"4","7":"GZMB"},{"1":"3.203437e-09","2":"0.9457311","3":"0.604","4":"0.108","5":"5.175793e-05","6":"4","7":"C12orf75"},{"1":"3.254137e-09","2":"1.4887908","3":"0.578","4":"0.127","5":"5.257710e-05","6":"4","7":"CMC1"},{"1":"4.000067e-09","2":"1.3713422","3":"0.707","4":"0.075","5":"6.462908e-05","6":"4","7":"KLRG1"},{"1":"4.888171e-09","2":"1.9770909","3":"0.589","4":"0.017","5":"7.897818e-05","6":"4","7":"KLRD1"},{"1":"1.015544e-08","2":"1.3530094","3":"0.532","4":"0.047","5":"1.640814e-04","6":"4","7":"IL2RB"},{"1":"2.016638e-08","2":"0.8274767","3":"0.632","4":"0.251","5":"3.258282e-04","6":"4","7":"ARPC5L"},{"1":"2.050895e-08","2":"1.4032752","3":"0.825","4":"0.296","5":"3.313631e-04","6":"4","7":"CD7"},{"1":"8.465147e-15","2":"2.3689367","3":"0.944","4":"0.006","5":"1.367714e-10","6":"5","7":"LILRA4"},{"1":"2.268992e-14","2":"2.9051890","3":"1.000","4":"0.063","5":"3.666011e-10","6":"5","7":"JCHAIN"},{"1":"3.122601e-14","2":"2.4405752","3":"0.944","4":"0.063","5":"5.045187e-10","6":"5","7":"PLD4"},{"1":"8.101092e-14","2":"1.9967267","3":"0.889","4":"0.017","5":"1.308893e-09","6":"5","7":"SERPINF1"},{"1":"8.101092e-14","2":"1.7806372","3":"0.889","4":"0.011","5":"1.308893e-09","6":"5","7":"TPM2"},{"1":"8.101092e-14","2":"1.5326684","3":"0.889","4":"0.007","5":"1.308893e-09","6":"5","7":"SMPD3"},{"1":"8.101092e-14","2":"1.4651406","3":"0.889","4":"0.000","5":"1.308893e-09","6":"5","7":"SCT"},{"1":"1.139709e-13","2":"2.0208807","3":"0.944","4":"0.052","5":"1.841427e-09","6":"5","7":"MZB1"},{"1":"1.223955e-13","2":"2.6036431","3":"1.000","4":"0.108","5":"1.977543e-09","6":"5","7":"ITM2C"},{"1":"3.646227e-13","2":"1.9692586","3":"0.944","4":"0.024","5":"5.891209e-09","6":"5","7":"DERL3"},{"1":"4.212910e-13","2":"1.6953386","3":"0.889","4":"0.008","5":"6.806799e-09","6":"5","7":"DNASE1L3"},{"1":"6.184215e-13","2":"2.6130189","3":"1.000","4":"0.126","5":"9.991836e-09","6":"5","7":"TCF4"},{"1":"7.283585e-13","2":"1.1678438","3":"0.833","4":"0.015","5":"1.176809e-08","6":"5","7":"GAS6"},{"1":"1.090582e-12","2":"1.9092540","3":"0.889","4":"0.025","5":"1.762054e-08","6":"5","7":"IL3RA"},{"1":"1.872633e-12","2":"1.5809879","3":"0.944","4":"0.099","5":"3.025614e-08","6":"5","7":"SPIB"},{"1":"3.041575e-12","2":"2.8134644","3":"0.944","4":"0.066","5":"4.914273e-08","6":"5","7":"GZMB"},{"1":"3.041575e-12","2":"2.0412797","3":"0.944","4":"0.091","5":"4.914273e-08","6":"5","7":"UGCG"},{"1":"3.242956e-12","2":"2.2862997","3":"0.889","4":"0.124","5":"5.239645e-08","6":"5","7":"CCDC50"},{"1":"6.164994e-12","2":"1.3351839","3":"0.778","4":"0.003","5":"9.960781e-08","6":"5","7":"CLEC4C"},{"1":"6.164994e-12","2":"1.3051923","3":"0.778","4":"0.000","5":"9.960781e-08","6":"5","7":"LRRC26"},{"1":"6.164994e-12","2":"1.0559316","3":"0.778","4":"0.006","5":"9.960781e-08","6":"5","7":"AL096865.1"},{"1":"6.164994e-12","2":"0.9020454","3":"0.778","4":"0.003","5":"9.960781e-08","6":"5","7":"TNFRSF21"},{"1":"1.215944e-11","2":"1.2614923","3":"0.889","4":"0.074","5":"1.964601e-07","6":"5","7":"TSPAN13"},{"1":"1.259429e-11","2":"1.0817647","3":"0.833","4":"0.083","5":"2.034859e-07","6":"5","7":"LILRB4"},{"1":"1.420112e-11","2":"1.6617392","3":"0.944","4":"0.178","5":"2.294475e-07","6":"5","7":"MAPKAPK2"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We can now select the top 25 up regulated genes for plotting.


```r
mypar(1, 5, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
    barplot(sort(setNames(top25$avg_logFC, top25$gene)[top25$cluster == i], F), horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
    abline(v = c(0, 0.25), lty = c(1, 2))
}
```

![](seurat_05_dge_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

We can visualize them as a heatmap. Here we are selecting the top 5.


```r
top5 <- markers_genes %>% group_by(cluster) %>% top_n(-5, p_val_adj)

alldata <- ScaleData(alldata, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(alldata, features = as.character(unique(top5$gene)), group.by = "kmeans_5", assay = "RNA")
```

![](seurat_05_dge_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

Another way is by representing the overal group expression and detection rates in a dot-plot.


```r
DotPlot(alldata, features = as.character(unique(top5$gene)), group.by = "kmeans_5", assay = "RNA") + coord_flip()
```

![](seurat_05_dge_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

We can also plot a violin plot for each gene.


```r
VlnPlot(alldata, features = as.character(unique(top5$gene)), ncol = 5, group.by = "kmeans_5", assay = "RNA")
```

![](seurat_05_dge_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

Take a screen shot of those results and re-run the same code above with another test: "wilcox" (Wilcoxon Rank Sum test), "bimod" (Likelihood-ratio test), "roc" (Identifies 'markers' of gene expression using ROC analysis),"t" (Student's t-test),"negbinom" (negative binomial generalized linear model),"poisson" (poisson generalized linear model), "LR" (logistic regression), "MAST" (hurdle model), "DESeq2" (negative binomial distribution).
</div>

## Differential expression across conditions
***

The second way of computing differential expression is to answer which genes are differentially expressed within a cluster. For example, in our case we have libraries comming from 2 different library preparation methods (batches) and we would like to know which genes are influenced the most in a particular cell type. The same concenpt applies if you have instead two or more biological groups (control vs treated, time#0 vs time#1 vs time#2, etc).

For this end, we will first subset our data for the desired cell cluster, then change the cell identities to the variable of comparison (which now in our case is the "Chemistry").


```r
cell_selection <- subset(alldata, cells = colnames(alldata)[alldata$kmeans_5 == 4])
cell_selection <- SetIdent(cell_selection, value = "Chemistry")
# Compute differentiall expression
DGE_cell_selection <- FindAllMarkers(cell_selection, logfc.threshold = 0.2, test.use = "wilcox", min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, assay = "RNA")
```

We can now plot the expression across the "Chemistry".


```r
top5_cell_selection <- DGE_cell_selection %>% group_by(cluster) %>% top_n(-5, p_val_adj)

VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)), ncol = 5, group.by = "Chemistry", assay = "RNA")
```

![](seurat_05_dge_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

We can clearly see some patterns across them. Those are the genes that impact the most on your batches (see the dimensionality reduction and integration exercises for more details). We can plot those genes using the integrated and non-integrated UMAP for ilustration.


```r
FeaturePlot(alldata, reduction = "UMAP_on_CCA", dims = 1:2, features = c("JUND", "RPS17", "CD81"), order = T, ncol = 3)
```

![](seurat_05_dge_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

Finally, lets save the integrated data for further analysis.


```r
saveRDS(alldata, "data/3pbmc_qc_dr_int_cl_dge.rds")
```


### Session Info
***


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.6
## 
## Matrix products: default
## BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] venn_1.7                    Seurat_3.0.2                dplyr_0.8.0.1               igraph_1.2.4.1              pheatmap_1.0.12             rafalib_1.0.0              
##  [7] cowplot_0.9.4               scran_1.10.2                scater_1.10.1               ggplot2_3.1.1               SingleCellExperiment_1.4.1  SummarizedExperiment_1.12.0
## [13] DelayedArray_0.8.0          BiocParallel_1.16.6         matrixStats_0.54.0          Biobase_2.42.0              GenomicRanges_1.34.0        GenomeInfoDb_1.18.2        
## [19] IRanges_2.16.0              S4Vectors_0.20.1            BiocGenerics_0.28.0        
## 
## loaded via a namespace (and not attached):
##   [1] Rtsne_0.15               ggbeeswarm_0.6.0         colorspace_1.4-1         ggridges_0.5.1           dynamicTreeCut_1.63-1    XVector_0.22.0           BiocNeighbors_1.0.0     
##   [8] rstudioapi_0.10          listenv_0.7.0            npsurv_0.4-0             ggrepel_0.8.0            codetools_0.2-16         splines_3.5.1            R.methodsS3_1.7.1       
##  [15] lsei_1.2-0               knitr_1.26               jsonlite_1.6             ica_1.0-2                cluster_2.0.9            png_0.1-7                R.oo_1.22.0             
##  [22] sctransform_0.2.0        HDF5Array_1.10.1         httr_1.4.0               compiler_3.5.1           assertthat_0.2.1         Matrix_1.2-17            lazyeval_0.2.2          
##  [29] limma_3.38.3             formatR_1.6              htmltools_0.4.0          tools_3.5.1              rsvd_1.0.0               gtable_0.3.0             glue_1.3.1              
##  [36] GenomeInfoDbData_1.2.0   RANN_2.6.1               reshape2_1.4.3           Rcpp_1.0.3               gdata_2.18.0             ape_5.3                  nlme_3.1-139            
##  [43] DelayedMatrixStats_1.4.0 gbRd_0.4-11              lmtest_0.9-37            xfun_0.11                stringr_1.4.0            globals_0.12.4           irlba_2.3.3             
##  [50] gtools_3.8.1             statmod_1.4.30           future_1.12.0            edgeR_3.24.3             zlibbioc_1.28.0          MASS_7.3-51.4            zoo_1.8-5               
##  [57] scales_1.0.0             rhdf5_2.26.2             RColorBrewer_1.1-2       yaml_2.2.0               reticulate_1.12          pbapply_1.4-0            gridExtra_2.3           
##  [64] stringi_1.4.3            caTools_1.17.1.2         bibtex_0.4.2             Rdpack_0.11-0            SDMTools_1.1-221.1       rlang_0.4.1              pkgconfig_2.0.2         
##  [71] bitops_1.0-6             evaluate_0.14            lattice_0.20-38          ROCR_1.0-7               purrr_0.3.2              Rhdf5lib_1.4.3           htmlwidgets_1.3         
##  [78] labeling_0.3             tidyselect_0.2.5         plyr_1.8.4               magrittr_1.5             R6_2.4.1                 gplots_3.0.1.1           pillar_1.3.1            
##  [85] withr_2.1.2              fitdistrplus_1.0-14      survival_2.44-1.1        RCurl_1.95-4.12          tsne_0.1-3               tibble_2.1.1             future.apply_1.2.0      
##  [92] crayon_1.3.4             KernSmooth_2.23-15       plotly_4.9.0             rmarkdown_1.17           viridis_0.5.1            locfit_1.5-9.1           grid_3.5.1              
##  [99] data.table_1.12.2        metap_1.1                digest_0.6.22            tidyr_0.8.3              R.utils_2.8.0            munsell_0.5.0            beeswarm_0.2.3          
## [106] viridisLite_0.3.0        vipor_0.4.5
```



















