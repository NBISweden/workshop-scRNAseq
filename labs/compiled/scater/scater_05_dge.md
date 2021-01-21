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

# Differential gene expression

In this tutorial we will cover about Differetial gene expression, which comprises an extensive range of topics and methods. In single cell, differential expresison can have multiple functionalities such as of identifying marker genes for cell populations, as well as differentially regulated genes across conditions (healthy vs control). We will also exercise on how to account the batch information in your test.

We can first load the data from the clustering session. Moreover, we can already decide which clustering resolution to use. First let's define using the `louvain` clustering to identifying differentially expressed genes.  


```r
suppressPackageStartupMessages({
    library(scater)
    library(scran)
    # library(venn)
    library(cowplot)
    library(ggplot2)
    # library(rafalib)
    library(pheatmap)
    library(igraph)
    library(dplyr)
})

sce <- readRDS("data/results/covid_qc_dr_int_cl.rds")
```

## Cell marker genes
***

Let us first compute a ranking for the highly differential genes in each cluster. There are many different tests and parameters to be chosen that can be used to refine your results. When looking for marker genes, we want genes that are positivelly expressed in a cell type and possibly not expressed in the others.


```r
# Compute differentiall expression
markers_genes <- scran::findMarkers(x = sce, groups = as.character(sce$louvain_SNNk15), 
    lfc = 0.5, pval.type = "all", direction = "up")

# List of dataFrames with the results for each cluster
markers_genes
```

```
## List of length 10
## names(10): 1 10 2 3 4 5 6 7 8 9
```

```r
# Visualizing the expression of one
markers_genes[["1"]]
```

```
## DataFrame with 18121 rows and 12 columns
##                p.value       FDR summary.logFC    logFC.10    logFC.2
##              <numeric> <numeric>     <numeric>   <numeric>  <numeric>
## CTSB       0.000029978  0.543231      0.698954     1.88904    2.00949
## FCN1       0.000403820  1.000000      0.683787     3.84707    3.74481
## CD14       0.001647679  1.000000      0.677626     2.69834    2.65354
## APLP2      0.004844297  1.000000      0.640151     2.50965    2.74728
## VCAN       0.007709632  1.000000      0.655263     3.24420    3.07794
## ...                ...       ...           ...         ...        ...
## AC011043.1           1         1    0.00238856  0.00238856 0.00238856
## AC007325.4           1         1    0.00617052  0.00617052 0.00715851
## AL354822.1           1         1    0.00272674 -0.00355950 0.00272674
## AC233755.1           1         1    0.00000000  0.00000000 0.00000000
## AC240274.1           1         1    0.00908758 -0.01389394 0.00908758
##                 logFC.3     logFC.4     logFC.5     logFC.6      logFC.7
##               <numeric>   <numeric>   <numeric>   <numeric>    <numeric>
## CTSB            1.90277     1.93142     1.91562    0.768169      1.79640
## FCN1            3.93218     3.93134     3.93130    0.975410      3.92790
## CD14            2.72778     2.72902     2.71753    2.418468      2.72469
## APLP2           2.50773     2.39105     2.20743    1.173234      2.59608
## VCAN            3.29866     3.30607     3.29558    3.026491      3.28994
## ...                 ...         ...         ...         ...          ...
## AC011043.1  0.000277072  0.00238856 -0.00254915  0.00238856  0.002388561
## AC007325.4  0.001027771  0.00316887  0.00313537 -0.01494578 -0.000270924
## AL354822.1 -0.002431376 -0.00984411 -0.01406130 -0.00136647 -0.003994630
## AC233755.1  0.000000000  0.00000000 -0.00562829  0.00000000  0.000000000
## AC240274.1  0.005880580  0.00366280 -0.00688946 -0.00440337 -0.000657017
##                logFC.8      logFC.9
##              <numeric>    <numeric>
## CTSB          0.698954      0.97215
## FCN1          0.683787      2.56441
## CD14          0.677626      2.45367
## APLP2         0.640151      1.53218
## VCAN          0.655263      2.47112
## ...                ...          ...
## AC011043.1  0.00238856 -0.005919522
## AC007325.4 -0.00116545  0.000426295
## AL354822.1  0.00177405  0.002726741
## AC233755.1  0.00000000  0.000000000
## AC240274.1 -0.00178147  0.000089034
```

We can now select the top 25 up regulated genes for plotting.


```r
# Colect the top 25 genes for each cluster and put the into a single table
top25 <- lapply(names(markers_genes), function(x) {
    temp <- markers_genes[[x]][1:25, 1:2]
    temp$gene <- rownames(markers_genes[[x]])[1:25]
    temp$cluster <- x
    return(temp)
})
top25 <- as_tibble(do.call(rbind, top25))
top25$p.value[top25$p.value == 0] <- 1e-300
top25
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["p.value"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["gene"],"name":[3],"type":["chr"],"align":["left"]},{"label":["cluster"],"name":[4],"type":["chr"],"align":["left"]}],"data":[{"1":"2.997796e-05","2":"5.432306e-01","3":"CTSB","4":"1"},{"1":"4.038203e-04","2":"1.000000e+00","3":"FCN1","4":"1"},{"1":"1.647679e-03","2":"1.000000e+00","3":"CD14","4":"1"},{"1":"4.844297e-03","2":"1.000000e+00","3":"APLP2","4":"1"},{"1":"7.709632e-03","2":"1.000000e+00","3":"VCAN","4":"1"},{"1":"8.239583e-02","2":"1.000000e+00","3":"LYZ","4":"1"},{"1":"1.664413e-01","2":"1.000000e+00","3":"CD36","4":"1"},{"1":"2.111126e-01","2":"1.000000e+00","3":"CEBPD","4":"1"},{"1":"3.220023e-01","2":"1.000000e+00","3":"ID1","4":"1"},{"1":"3.952893e-01","2":"1.000000e+00","3":"NCF2","4":"1"},{"1":"4.601536e-01","2":"1.000000e+00","3":"TALDO1","4":"1"},{"1":"5.715863e-01","2":"1.000000e+00","3":"CSTA","4":"1"},{"1":"6.242592e-01","2":"1.000000e+00","3":"MS4A6A","4":"1"},{"1":"6.293161e-01","2":"1.000000e+00","3":"GRN","4":"1"},{"1":"6.799593e-01","2":"1.000000e+00","3":"LGALS1","4":"1"},{"1":"6.962459e-01","2":"1.000000e+00","3":"TNFAIP2","4":"1"},{"1":"7.165509e-01","2":"1.000000e+00","3":"CAST","4":"1"},{"1":"8.777894e-01","2":"1.000000e+00","3":"METTL9","4":"1"},{"1":"8.981159e-01","2":"1.000000e+00","3":"S100A9","4":"1"},{"1":"9.115374e-01","2":"1.000000e+00","3":"LRP1","4":"1"},{"1":"9.127933e-01","2":"1.000000e+00","3":"ATP5MPL","4":"1"},{"1":"9.313998e-01","2":"1.000000e+00","3":"CYBB","4":"1"},{"1":"9.530949e-01","2":"1.000000e+00","3":"BLVRB","4":"1"},{"1":"9.610846e-01","2":"1.000000e+00","3":"MPEG1","4":"1"},{"1":"9.677807e-01","2":"1.000000e+00","3":"PKM","4":"1"},{"1":"1.000000e+00","2":"1.000000e+00","3":"IGHA1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"MTRNR2L12","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"HNRNPH1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"STMN1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"DDX17","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"CCND2","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"HIST1H4C","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"RNF213","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"OGT","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"MACF1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"HIST1H1C","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"SRSF11","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"PCNA","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"HMGB1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"ATRX","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"EPB41","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"PRDM1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"ATM","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"DDX46","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"XPO1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"SORL1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"HSPD1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"TYMS","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"MBNL1","4":"10"},{"1":"1.000000e+00","2":"1.000000e+00","3":"C1orf56","4":"10"},{"1":"8.718208e-46","2":"1.579827e-41","3":"PF4","4":"2"},{"1":"2.362160e-28","2":"2.140235e-24","3":"OST4","4":"2"},{"1":"1.949388e-19","2":"1.177495e-15","3":"PPBP","4":"2"},{"1":"2.698002e-16","2":"1.222262e-12","3":"H3F3A","4":"2"},{"1":"4.257611e-15","2":"1.543043e-11","3":"NRGN","4":"2"},{"1":"5.610571e-15","2":"1.694486e-11","3":"GNG11","4":"2"},{"1":"3.011734e-13","2":"7.796519e-10","3":"PDLIM1","4":"2"},{"1":"8.953589e-12","2":"2.028100e-08","3":"MYL9","4":"2"},{"1":"3.585651e-11","2":"7.219509e-08","3":"HIST1H2AC","4":"2"},{"1":"1.303039e-10","2":"2.361237e-07","3":"CAVIN2","4":"2"},{"1":"5.050671e-10","2":"8.320293e-07","3":"HIST1H3H","4":"2"},{"1":"7.856060e-10","2":"1.186331e-06","3":"TREML1","4":"2"},{"1":"1.591340e-09","2":"2.218206e-06","3":"GP9","4":"2"},{"1":"1.905318e-09","2":"2.466162e-06","3":"CA2","4":"2"},{"1":"3.566674e-09","2":"4.308780e-06","3":"CLEC1B","4":"2"},{"1":"4.445888e-09","2":"5.035246e-06","3":"SNCA","4":"2"},{"1":"6.915976e-09","2":"7.372023e-06","3":"TSC22D1","4":"2"},{"1":"8.612388e-09","2":"8.670283e-06","3":"RGS18","4":"2"},{"1":"1.376494e-08","2":"1.312813e-05","3":"CMTM5","4":"2"},{"1":"1.560002e-08","2":"1.413440e-05","3":"TUBB1","4":"2"},{"1":"1.909591e-08","2":"1.647795e-05","3":"ACRBP","4":"2"},{"1":"4.029158e-08","2":"3.318744e-05","3":"RSU1","4":"2"},{"1":"1.339924e-07","2":"1.055686e-04","3":"TAGLN2","4":"2"},{"1":"1.437030e-07","2":"1.085018e-04","3":"RGS10","4":"2"},{"1":"2.514815e-07","2":"1.822839e-04","3":"MTURN","4":"2"},{"1":"6.227074e-16","2":"1.128408e-11","3":"CD8A","4":"3"},{"1":"8.665600e-14","2":"7.851467e-10","3":"TRGC2","4":"3"},{"1":"2.286893e-13","2":"1.381360e-09","3":"DUSP2","4":"3"},{"1":"1.121009e-10","2":"5.078450e-07","3":"GZMK","4":"3"},{"1":"9.462334e-10","2":"3.429339e-06","3":"LYAR","4":"3"},{"1":"1.184235e-07","2":"3.576587e-04","3":"KLRG1","4":"3"},{"1":"1.861769e-03","2":"1.000000e+00","3":"CD8B","4":"3"},{"1":"6.504178e-01","2":"1.000000e+00","3":"CD3D","4":"3"},{"1":"7.942264e-01","2":"1.000000e+00","3":"PIK3R1","4":"3"},{"1":"8.728181e-01","2":"1.000000e+00","3":"TNFAIP3","4":"3"},{"1":"8.729507e-01","2":"1.000000e+00","3":"IL32","4":"3"},{"1":"8.923725e-01","2":"1.000000e+00","3":"SRSF7","4":"3"},{"1":"9.914247e-01","2":"1.000000e+00","3":"TUBA4A","4":"3"},{"1":"9.956414e-01","2":"1.000000e+00","3":"CD3G","4":"3"},{"1":"9.981542e-01","2":"1.000000e+00","3":"LINC01871","4":"3"},{"1":"9.983029e-01","2":"1.000000e+00","3":"CCL5","4":"3"},{"1":"9.999998e-01","2":"1.000000e+00","3":"PPP2R5C","4":"3"},{"1":"9.999999e-01","2":"1.000000e+00","3":"CD3E","4":"3"},{"1":"9.999999e-01","2":"1.000000e+00","3":"CD2","4":"3"},{"1":"1.000000e+00","2":"1.000000e+00","3":"ATG2A","4":"3"},{"1":"1.000000e+00","2":"1.000000e+00","3":"RPS26","4":"3"},{"1":"1.000000e+00","2":"1.000000e+00","3":"RNF19A","4":"3"},{"1":"1.000000e+00","2":"1.000000e+00","3":"RORA","4":"3"},{"1":"1.000000e+00","2":"1.000000e+00","3":"A2M-AS1","4":"3"},{"1":"1.000000e+00","2":"1.000000e+00","3":"HSP90AA1","4":"3"},{"1":"8.018089e-62","2":"1.452958e-57","3":"GNLY","4":"4"},{"1":"1.079846e-48","2":"9.783943e-45","3":"KLRF1","4":"4"},{"1":"3.288739e-41","2":"1.986508e-37","3":"FGFBP2","4":"4"},{"1":"3.027512e-34","2":"1.371539e-30","3":"CD247","4":"4"},{"1":"4.045455e-34","2":"1.466154e-30","3":"TRDC","4":"4"},{"1":"9.627926e-27","2":"2.907794e-23","3":"GZMB","4":"4"},{"1":"5.813728e-25","2":"1.505008e-21","3":"CTSW","4":"4"},{"1":"1.715792e-21","2":"3.886483e-18","3":"PRF1","4":"4"},{"1":"1.683407e-19","2":"3.389447e-16","3":"KLRD1","4":"4"},{"1":"3.891479e-17","2":"7.051749e-14","3":"SPON2","4":"4"},{"1":"7.093334e-16","2":"1.168530e-12","3":"NKG7","4":"4"},{"1":"3.628230e-12","2":"5.478929e-09","3":"MYOM2","4":"4"},{"1":"7.681339e-11","2":"1.070720e-07","3":"CD7","4":"4"},{"1":"1.776101e-08","2":"2.298908e-05","3":"KLRB1","4":"4"},{"1":"8.844454e-08","2":"1.068469e-04","3":"CLIC3","4":"4"},{"1":"1.207212e-05","2":"1.367243e-02","3":"HOPX","4":"4"},{"1":"1.608913e-03","2":"1.000000e+00","3":"APMAP","4":"4"},{"1":"1.748492e-03","2":"1.000000e+00","3":"IL2RB","4":"4"},{"1":"1.550639e-02","2":"1.000000e+00","3":"ABHD17A","4":"4"},{"1":"2.278209e-02","2":"1.000000e+00","3":"PRSS23","4":"4"},{"1":"3.995857e-02","2":"1.000000e+00","3":"SYNGR1","4":"4"},{"1":"9.754273e-02","2":"1.000000e+00","3":"GZMA","4":"4"},{"1":"2.029523e-01","2":"1.000000e+00","3":"CST7","4":"4"},{"1":"5.883317e-01","2":"1.000000e+00","3":"S1PR5","4":"4"},{"1":"6.451925e-01","2":"1.000000e+00","3":"KLRC2","4":"4"},{"1":"1.032045e-185","2":"1.870170e-181","3":"MS4A1","4":"5"},{"1":"4.237377e-185","2":"3.839276e-181","3":"CD79A","4":"5"},{"1":"6.963845e-124","2":"4.206394e-120","3":"LINC00926","4":"5"},{"1":"3.060675e-102","2":"1.386562e-98","3":"IGHD","4":"5"},{"1":"3.344243e-101","2":"1.212021e-97","3":"TNFRSF13C","4":"5"},{"1":"4.576046e-59","2":"1.382042e-55","3":"IGHM","4":"5"},{"1":"1.118940e-36","2":"2.896617e-33","3":"VPREB3","4":"5"},{"1":"4.164175e-34","2":"9.432378e-31","3":"CD79B","4":"5"},{"1":"7.400590e-34","2":"1.490068e-30","3":"RALGPS2","4":"5"},{"1":"2.071529e-33","2":"3.753818e-30","3":"CD22","4":"5"},{"1":"7.244453e-30","2":"1.193425e-26","3":"P2RX5","4":"5"},{"1":"1.846104e-27","2":"2.787772e-24","3":"JUND","4":"5"},{"1":"1.174335e-22","2":"1.636933e-19","3":"CD37","4":"5"},{"1":"2.327643e-20","2":"3.012801e-17","3":"FCER2","4":"5"},{"1":"3.317792e-18","2":"4.008114e-15","3":"NFKBID","4":"5"},{"1":"3.169576e-16","2":"3.589742e-13","3":"BLK","4":"5"},{"1":"2.003756e-15","2":"2.135886e-12","3":"TAGAP","4":"5"},{"1":"7.561358e-14","2":"7.612187e-11","3":"TCL1A","4":"5"},{"1":"2.753494e-13","2":"2.626108e-10","3":"BANK1","4":"5"},{"1":"1.173080e-11","2":"1.062869e-08","3":"FAM129C","4":"5"},{"1":"5.584700e-10","2":"4.819064e-07","3":"CD83","4":"5"},{"1":"4.843826e-07","2":"3.989772e-04","3":"FCRL1","4":"5"},{"1":"2.447968e-06","2":"1.928680e-03","3":"PAX5","4":"5"},{"1":"9.541769e-06","2":"6.941553e-03","3":"PRDM2","4":"5"},{"1":"9.576669e-06","2":"6.941553e-03","3":"EZR","4":"5"},{"1":"1.707907e-58","2":"3.094898e-54","3":"CDKN1C","4":"6"},{"1":"4.048937e-34","2":"3.668540e-30","3":"FCGR3A","4":"6"},{"1":"1.076491e-28","2":"6.502362e-25","3":"SMIM25","4":"6"},{"1":"3.282204e-23","2":"1.486920e-19","3":"MS4A7","4":"6"},{"1":"2.909746e-22","2":"1.054550e-18","3":"LRRC25","4":"6"},{"1":"4.172282e-21","2":"1.260099e-17","3":"LST1","4":"6"},{"1":"5.059392e-15","2":"1.309732e-11","3":"AIF1","4":"6"},{"1":"7.180112e-15","2":"1.626385e-11","3":"TCF7L2","4":"6"},{"1":"8.404017e-12","2":"1.692102e-08","3":"CTSC","4":"6"},{"1":"4.999599e-11","2":"9.059773e-08","3":"HES4","4":"6"},{"1":"1.209046e-09","2":"1.991738e-06","3":"PECAM1","4":"6"},{"1":"1.461928e-07","2":"2.207633e-04","3":"SLC2A6","4":"6"},{"1":"4.046769e-07","2":"5.640884e-04","3":"WARS","4":"6"},{"1":"2.930376e-06","2":"3.792953e-03","3":"SERPINA1","4":"6"},{"1":"4.155294e-06","2":"5.019872e-03","3":"IFITM3","4":"6"},{"1":"1.208420e-05","2":"1.368611e-02","3":"CSF1R","4":"6"},{"1":"6.112600e-05","2":"6.515672e-02","3":"MBD2","4":"6"},{"1":"1.723486e-04","2":"1.735072e-01","3":"TUBA1A","4":"6"},{"1":"2.072000e-04","2":"1.976142e-01","3":"DRAP1","4":"6"},{"1":"2.591573e-04","2":"2.348095e-01","3":"TNFRSF1B","4":"6"},{"1":"3.982349e-04","2":"3.436388e-01","3":"COTL1","4":"6"},{"1":"6.567524e-04","2":"5.310777e-01","3":"IER5","4":"6"},{"1":"6.740681e-04","2":"5.310777e-01","3":"SIGLEC10","4":"6"},{"1":"8.496742e-04","2":"6.239966e-01","3":"SPN","4":"6"},{"1":"8.608750e-04","2":"6.239966e-01","3":"LYPD2","4":"6"},{"1":"1.069543e-53","2":"1.938119e-49","3":"IL7R","4":"7"},{"1":"5.367164e-18","2":"4.862919e-14","3":"NOSIP","4":"7"},{"1":"2.884376e-16","2":"1.742259e-12","3":"MAL","4":"7"},{"1":"1.184593e-14","2":"5.366505e-11","3":"PIK3IP1","4":"7"},{"1":"6.261336e-14","2":"2.269233e-10","3":"SARAF","4":"7"},{"1":"1.031927e-12","2":"3.116591e-09","3":"LEPROTL1","4":"7"},{"1":"3.564069e-12","2":"9.226356e-09","3":"LDHB","4":"7"},{"1":"4.470369e-09","2":"1.012594e-05","3":"TCF7","4":"7"},{"1":"3.250354e-08","2":"6.544406e-05","3":"RCAN3","4":"7"},{"1":"9.230497e-06","2":"1.672658e-02","3":"TPT1","4":"7"},{"1":"4.774385e-05","2":"7.865148e-02","3":"RPS14","4":"7"},{"1":"6.662959e-04","2":"1.000000e+00","3":"RPS4X","4":"7"},{"1":"2.538595e-03","2":"1.000000e+00","3":"RPS12","4":"7"},{"1":"4.444296e-03","2":"1.000000e+00","3":"LTB","4":"7"},{"1":"1.605563e-02","2":"1.000000e+00","3":"EEF1A1","4":"7"},{"1":"2.671978e-02","2":"1.000000e+00","3":"AP3M2","4":"7"},{"1":"3.795444e-02","2":"1.000000e+00","3":"KLF2","4":"7"},{"1":"3.970998e-02","2":"1.000000e+00","3":"RPS15A","4":"7"},{"1":"7.314385e-02","2":"1.000000e+00","3":"RPL14","4":"7"},{"1":"9.012479e-02","2":"1.000000e+00","3":"RPL10","4":"7"},{"1":"9.197999e-02","2":"1.000000e+00","3":"ARHGAP15","4":"7"},{"1":"1.716454e-01","2":"1.000000e+00","3":"RPL3","4":"7"},{"1":"3.479385e-01","2":"1.000000e+00","3":"AQP3","4":"7"},{"1":"3.755028e-01","2":"1.000000e+00","3":"RPS3","4":"7"},{"1":"4.205351e-01","2":"1.000000e+00","3":"RPS25","4":"7"},{"1":"2.530389e-02","2":"1.000000e+00","3":"CXCL8","4":"8"},{"1":"5.624830e-01","2":"1.000000e+00","3":"G0S2","4":"8"},{"1":"6.743522e-01","2":"1.000000e+00","3":"S100A8","4":"8"},{"1":"8.404227e-01","2":"1.000000e+00","3":"NAMPT","4":"8"},{"1":"9.999566e-01","2":"1.000000e+00","3":"IFI27","4":"8"},{"1":"9.999773e-01","2":"1.000000e+00","3":"S100A12","4":"8"},{"1":"9.999975e-01","2":"1.000000e+00","3":"IL1B","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"RBP7","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"ISG15","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"SOD2","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"CCL3","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"MARCKS","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"BCL2A1","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"EGR1","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"RGS2","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"PTGS2","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"ALOX5AP","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"ATP2B1-AS1","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"NCF1","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"IFI6","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"GCA","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"NFKBIA","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"PLAC8","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"FOS","4":"8"},{"1":"1.000000e+00","2":"1.000000e+00","3":"DUSP1","4":"8"},{"1":"3.850809e-14","2":"6.978051e-10","3":"HLA-DQA1","4":"9"},{"1":"1.651203e-10","2":"1.089378e-06","3":"FCER1A","4":"9"},{"1":"1.803506e-10","2":"1.089378e-06","3":"HLA-DPA1","4":"9"},{"1":"5.602608e-08","2":"2.538121e-04","3":"HLA-DPB1","4":"9"},{"1":"3.310505e-07","2":"1.199793e-03","3":"HLA-DRA","4":"9"},{"1":"9.979741e-07","2":"3.014048e-03","3":"HLA-DRB1","4":"9"},{"1":"2.305512e-06","2":"5.968313e-03","3":"CLEC10A","4":"9"},{"1":"4.266912e-06","2":"9.665090e-03","3":"HLA-DMA","4":"9"},{"1":"8.416000e-05","2":"1.694515e-01","3":"ENHO","4":"9"},{"1":"6.219458e-04","2":"1.000000e+00","3":"CD1C","4":"9"},{"1":"8.173754e-04","2":"1.000000e+00","3":"KCNK6","4":"9"},{"1":"9.827081e-04","2":"1.000000e+00","3":"HLA-DQB1","4":"9"},{"1":"1.044303e-01","2":"1.000000e+00","3":"CST3","4":"9"},{"1":"1.099854e-01","2":"1.000000e+00","3":"NDRG2","4":"9"},{"1":"1.529525e-01","2":"1.000000e+00","3":"PPA1","4":"9"},{"1":"1.753818e-01","2":"1.000000e+00","3":"CPVL","4":"9"},{"1":"2.592021e-01","2":"1.000000e+00","3":"ACTG1","4":"9"},{"1":"2.728597e-01","2":"1.000000e+00","3":"HLA-DRB5","4":"9"},{"1":"2.729636e-01","2":"1.000000e+00","3":"MAT2A","4":"9"},{"1":"2.872657e-01","2":"1.000000e+00","3":"EIF3L","4":"9"},{"1":"3.025709e-01","2":"1.000000e+00","3":"PLD4","4":"9"},{"1":"3.344362e-01","2":"1.000000e+00","3":"LGALS2","4":"9"},{"1":"4.146244e-01","2":"1.000000e+00","3":"TMEM14C","4":"9"},{"1":"5.043426e-01","2":"1.000000e+00","3":"SLC25A3","4":"9"},{"1":"5.172306e-01","2":"1.000000e+00","3":"CD74","4":"9"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We can now select the top 25 up regulated genes for plotting.


```r
par(mfrow = c(1, 5), mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
    barplot(sort(setNames(-log10(top25$p.value), top25$gene)[top25$cluster == i], 
        F), horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", 
        yaxs = "i", xlab = "-log10FC")
    abline(v = c(0, -log10(0.05)), lty = c(1, 2))
}
```

![](scater_05_dge_files/figure-html/unnamed-chunk-4-1.png)<!-- -->![](scater_05_dge_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

We can visualize them as a heatmap. Here we are selecting the top 5.


```r
top5 <- as_tibble(top25) %>% group_by(cluster) %>% top_n(-5, p.value)

scater::plotHeatmap(sce[, order(sce$louvain_SNNk15)], features = unique(top5$gene), 
    center = T, zlim = c(-3, 3), colour_columns_by = "louvain_SNNk15", show_colnames = F, 
    cluster_cols = F, fontsize_row = 6, color = colorRampPalette(c("purple", "black", 
        "yellow"))(90))
```

![](scater_05_dge_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

We can also plot a violin plot for each gene.


```r
scater::plotExpression(sce, features = unique(top5$gene), x = "louvain_SNNk15", ncol = 5, 
    colour_by = "louvain_SNNk15", scales = "free")
```

![](scater_05_dge_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


## Differential expression across conditions
***

The second way of computing differential expression is to answer which genes are differentially expressed within a cluster. For example, in our case we have libraries comming from patients and controls and we would like to know which genes are influenced the most in a particular cell type.

For this end, we will first subset our data for the desired cell cluster, then change the cell identities to the variable of comparison (which now in our case is the "type", e.g. Covid/Ctrl).


```r
# Filter cells from that cluster
cell_selection <- sce[, sce$louvain_SNNk15 == 6]

# Compute differentiall expression
DGE_cell_selection <- findMarkers(x = cell_selection, groups = cell_selection@colData$type, 
    lfc = 0.25, pval.type = "all", direction = "any")
top5_cell_selection <- lapply(names(DGE_cell_selection), function(x) {
    temp <- DGE_cell_selection[[x]][1:5, 1:2]
    temp$gene <- rownames(DGE_cell_selection[[x]])[1:5]
    temp$cluster <- x
    return(temp)
})
top5_cell_selection <- as_tibble(do.call(rbind, top5_cell_selection))
top5_cell_selection
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["p.value"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["gene"],"name":[3],"type":["chr"],"align":["left"]},{"label":["cluster"],"name":[4],"type":["chr"],"align":["left"]}],"data":[{"1":"8.918461e-27","2":"1.616114e-22","3":"XIST","4":"Control"},{"1":"2.396754e-11","2":"2.171579e-07","3":"NFKBIA","4":"Control"},{"1":"8.170255e-11","2":"4.935107e-07","3":"SOD2","4":"Control"},{"1":"6.085548e-10","2":"2.756906e-06","3":"RPS4Y1","4":"Control"},{"1":"2.625642e-09","2":"9.515850e-06","3":"IFITM3","4":"Control"},{"1":"8.918461e-27","2":"1.616114e-22","3":"XIST","4":"Covid"},{"1":"2.396754e-11","2":"2.171579e-07","3":"NFKBIA","4":"Covid"},{"1":"8.170255e-11","2":"4.935107e-07","3":"SOD2","4":"Covid"},{"1":"6.085548e-10","2":"2.756906e-06","3":"RPS4Y1","4":"Covid"},{"1":"2.625642e-09","2":"9.515850e-06","3":"IFITM3","4":"Covid"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

We can now plot the expression across the "type".


```r
scater::plotExpression(cell_selection, features = unique(top5_cell_selection$gene), 
    x = "type", ncol = 5, colour_by = "type")
```

![](scater_05_dge_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

#DGE_ALL6.2:


```r
plotlist <- list()
for (i in unique(top5_cell_selection$gene)) {
    plotlist[[i]] <- plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = i, by_exprs_values = "logcounts") + 
        ggtitle(label = i) + theme(plot.title = element_text(size = 20))
}
plot_grid(ncol = 3, plotlist = plotlist)
```

![](scater_05_dge_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


## Gene Set Analysis
***

Hypergeometric enrichment test

Having a defined list of differentially expressed genes, you can now look for their combined function using hypergeometric test:


```r
# Load additional packages
library(enrichR)

# Check available databases to perform enrichment (then choose one)
enrichR::listEnrichrDbs()
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["geneCoverage"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["genesPerTerm"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["libraryName"],"name":[3],"type":["chr"],"align":["left"]},{"label":["link"],"name":[4],"type":["chr"],"align":["left"]},{"label":["numTerms"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"13362","2":"275","3":"Genome_Browser_PWMs","4":"http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/","5":"615"},{"1":"27884","2":"1284","3":"TRANSFAC_and_JASPAR_PWMs","4":"http://jaspar.genereg.net/html/DOWNLOAD/","5":"326"},{"1":"6002","2":"77","3":"Transcription_Factor_PPIs","4":"","5":"290"},{"1":"47172","2":"1370","3":"ChEA_2013","4":"http://amp.pharm.mssm.edu/lib/cheadownload.jsp","5":"353"},{"1":"47107","2":"509","3":"Drug_Perturbations_from_GEO_2014","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"701"},{"1":"21493","2":"3713","3":"ENCODE_TF_ChIP-seq_2014","4":"http://genome.ucsc.edu/ENCODE/downloads.html","5":"498"},{"1":"1295","2":"18","3":"BioCarta_2013","4":"https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways","5":"249"},{"1":"3185","2":"73","3":"Reactome_2013","4":"http://www.reactome.org/download/index.html","5":"78"},{"1":"2854","2":"34","3":"WikiPathways_2013","4":"http://www.wikipathways.org/index.php/Download_Pathways","5":"199"},{"1":"15057","2":"300","3":"Disease_Signatures_from_GEO_up_2014","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"142"},{"1":"4128","2":"48","3":"KEGG_2013","4":"http://www.kegg.jp/kegg/download/","5":"200"},{"1":"34061","2":"641","3":"TF-LOF_Expression_from_GEO","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"269"},{"1":"7504","2":"155","3":"TargetScan_microRNA","4":"http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_61","5":"222"},{"1":"16399","2":"247","3":"PPI_Hub_Proteins","4":"http://amp.pharm.mssm.edu/X2K","5":"385"},{"1":"12753","2":"57","3":"GO_Molecular_Function_2015","4":"http://www.geneontology.org/GO.downloads.annotations.shtml","5":"1136"},{"1":"23726","2":"127","3":"GeneSigDB","4":"http://genesigdb.org/genesigdb/downloadall.jsp","5":"2139"},{"1":"32740","2":"85","3":"Chromosome_Location","4":"http://software.broadinstitute.org/gsea/msigdb/index.jsp","5":"386"},{"1":"13373","2":"258","3":"Human_Gene_Atlas","4":"http://biogps.org/downloads/","5":"84"},{"1":"19270","2":"388","3":"Mouse_Gene_Atlas","4":"http://biogps.org/downloads/","5":"96"},{"1":"13236","2":"82","3":"GO_Cellular_Component_2015","4":"http://www.geneontology.org/GO.downloads.annotations.shtml","5":"641"},{"1":"14264","2":"58","3":"GO_Biological_Process_2015","4":"http://www.geneontology.org/GO.downloads.annotations.shtml","5":"5192"},{"1":"3096","2":"31","3":"Human_Phenotype_Ontology","4":"http://www.human-phenotype-ontology.org/","5":"1779"},{"1":"22288","2":"4368","3":"Epigenomics_Roadmap_HM_ChIP-seq","4":"http://www.roadmapepigenomics.org/","5":"383"},{"1":"4533","2":"37","3":"KEA_2013","4":"http://amp.pharm.mssm.edu/lib/keacommandline.jsp","5":"474"},{"1":"10231","2":"158","3":"NURSA_Human_Endogenous_Complexome","4":"https://www.nursa.org/nursa/index.jsf","5":"1796"},{"1":"2741","2":"5","3":"CORUM","4":"http://mips.helmholtz-muenchen.de/genre/proj/corum/","5":"1658"},{"1":"5655","2":"342","3":"SILAC_Phosphoproteomics","4":"http://amp.pharm.mssm.edu/lib/keacommandline.jsp","5":"84"},{"1":"10406","2":"715","3":"MGI_Mammalian_Phenotype_Level_3","4":"http://www.informatics.jax.org/","5":"71"},{"1":"10493","2":"200","3":"MGI_Mammalian_Phenotype_Level_4","4":"http://www.informatics.jax.org/","5":"476"},{"1":"11251","2":"100","3":"Old_CMAP_up","4":"http://www.broadinstitute.org/cmap/","5":"6100"},{"1":"8695","2":"100","3":"Old_CMAP_down","4":"http://www.broadinstitute.org/cmap/","5":"6100"},{"1":"1759","2":"25","3":"OMIM_Disease","4":"http://www.omim.org/downloads","5":"90"},{"1":"2178","2":"89","3":"OMIM_Expanded","4":"http://www.omim.org/downloads","5":"187"},{"1":"851","2":"15","3":"VirusMINT","4":"http://mint.bio.uniroma2.it/download.html","5":"85"},{"1":"10061","2":"106","3":"MSigDB_Computational","4":"http://www.broadinstitute.org/gsea/msigdb/collections.jsp","5":"858"},{"1":"11250","2":"166","3":"MSigDB_Oncogenic_Signatures","4":"http://www.broadinstitute.org/gsea/msigdb/collections.jsp","5":"189"},{"1":"15406","2":"300","3":"Disease_Signatures_from_GEO_down_2014","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"142"},{"1":"17711","2":"300","3":"Virus_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"323"},{"1":"17576","2":"300","3":"Virus_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"323"},{"1":"15797","2":"176","3":"Cancer_Cell_Line_Encyclopedia","4":"https://portals.broadinstitute.org/ccle/home\\n","5":"967"},{"1":"12232","2":"343","3":"NCI-60_Cancer_Cell_Lines","4":"http://biogps.org/downloads/","5":"93"},{"1":"13572","2":"301","3":"Tissue_Protein_Expression_from_ProteomicsDB","4":"https://www.proteomicsdb.org/","5":"207"},{"1":"6454","2":"301","3":"Tissue_Protein_Expression_from_Human_Proteome_Map","4":"http://www.humanproteomemap.org/index.php","5":"30"},{"1":"3723","2":"47","3":"HMDB_Metabolites","4":"http://www.hmdb.ca/downloads","5":"3906"},{"1":"7588","2":"35","3":"Pfam_InterPro_Domains","4":"ftp://ftp.ebi.ac.uk/pub/databases/interpro/","5":"311"},{"1":"7682","2":"78","3":"GO_Biological_Process_2013","4":"http://www.geneontology.org/GO.downloads.annotations.shtml","5":"941"},{"1":"7324","2":"172","3":"GO_Cellular_Component_2013","4":"http://www.geneontology.org/GO.downloads.annotations.shtml","5":"205"},{"1":"8469","2":"122","3":"GO_Molecular_Function_2013","4":"http://www.geneontology.org/GO.downloads.annotations.shtml","5":"402"},{"1":"13121","2":"305","3":"Allen_Brain_Atlas_up","4":"http://www.brain-map.org/","5":"2192"},{"1":"26382","2":"1811","3":"ENCODE_TF_ChIP-seq_2015","4":"http://genome.ucsc.edu/ENCODE/downloads.html","5":"816"},{"1":"29065","2":"2123","3":"ENCODE_Histone_Modifications_2015","4":"http://genome.ucsc.edu/ENCODE/downloads.html","5":"412"},{"1":"280","2":"9","3":"Phosphatase_Substrates_from_DEPOD","4":"http://www.koehn.embl.de/depod/","5":"59"},{"1":"13877","2":"304","3":"Allen_Brain_Atlas_down","4":"http://www.brain-map.org/","5":"2192"},{"1":"15852","2":"912","3":"ENCODE_Histone_Modifications_2013","4":"http://genome.ucsc.edu/ENCODE/downloads.html","5":"109"},{"1":"4320","2":"129","3":"Achilles_fitness_increase","4":"http://www.broadinstitute.org/achilles","5":"216"},{"1":"4271","2":"128","3":"Achilles_fitness_decrease","4":"http://www.broadinstitute.org/achilles","5":"216"},{"1":"10496","2":"201","3":"MGI_Mammalian_Phenotype_2013","4":"http://www.informatics.jax.org/","5":"476"},{"1":"1678","2":"21","3":"BioCarta_2015","4":"https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways","5":"239"},{"1":"756","2":"12","3":"HumanCyc_2015","4":"http://humancyc.org/","5":"125"},{"1":"3800","2":"48","3":"KEGG_2015","4":"http://www.kegg.jp/kegg/download/","5":"179"},{"1":"2541","2":"39","3":"NCI-Nature_2015","4":"http://pid.nci.nih.gov/","5":"209"},{"1":"1918","2":"39","3":"Panther_2015","4":"http://www.pantherdb.org/","5":"104"},{"1":"5863","2":"51","3":"WikiPathways_2015","4":"http://www.wikipathways.org/index.php/Download_Pathways","5":"404"},{"1":"6768","2":"47","3":"Reactome_2015","4":"http://www.reactome.org/download/index.html","5":"1389"},{"1":"25651","2":"807","3":"ESCAPE","4":"http://www.maayanlab.net/ESCAPE/","5":"315"},{"1":"19129","2":"1594","3":"HomoloGene","4":"http://www.ncbi.nlm.nih.gov/homologene","5":"12"},{"1":"23939","2":"293","3":"Disease_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"839"},{"1":"23561","2":"307","3":"Disease_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"839"},{"1":"23877","2":"302","3":"Drug_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"906"},{"1":"15886","2":"9","3":"Genes_Associated_with_NIH_Grants","4":"https://grants.nih.gov/grants/oer.htm\\n","5":"32876"},{"1":"24350","2":"299","3":"Drug_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"906"},{"1":"3102","2":"25","3":"KEA_2015","4":"http://amp.pharm.mssm.edu/Enrichr","5":"428"},{"1":"31132","2":"298","3":"Gene_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"2460"},{"1":"30832","2":"302","3":"Gene_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"2460"},{"1":"48230","2":"1429","3":"ChEA_2015","4":"http://amp.pharm.mssm.edu/Enrichr","5":"395"},{"1":"5613","2":"36","3":"dbGaP","4":"http://www.ncbi.nlm.nih.gov/gap","5":"345"},{"1":"9559","2":"73","3":"LINCS_L1000_Chem_Pert_up","4":"https://clue.io/","5":"33132"},{"1":"9448","2":"63","3":"LINCS_L1000_Chem_Pert_down","4":"https://clue.io/","5":"33132"},{"1":"16725","2":"1443","3":"GTEx_Tissue_Sample_Gene_Expression_Profiles_down","4":"http://www.gtexportal.org/","5":"2918"},{"1":"19249","2":"1443","3":"GTEx_Tissue_Sample_Gene_Expression_Profiles_up","4":"http://www.gtexportal.org/","5":"2918"},{"1":"15090","2":"282","3":"Ligand_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"261"},{"1":"16129","2":"292","3":"Aging_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"286"},{"1":"15309","2":"308","3":"Aging_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"286"},{"1":"15103","2":"318","3":"Ligand_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"261"},{"1":"15022","2":"290","3":"MCF7_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"401"},{"1":"15676","2":"310","3":"MCF7_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"401"},{"1":"15854","2":"279","3":"Microbe_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"312"},{"1":"15015","2":"321","3":"Microbe_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"312"},{"1":"3788","2":"159","3":"LINCS_L1000_Ligand_Perturbations_down","4":"https://clue.io/","5":"96"},{"1":"3357","2":"153","3":"LINCS_L1000_Ligand_Perturbations_up","4":"https://clue.io/","5":"96"},{"1":"12668","2":"300","3":"L1000_Kinase_and_GPCR_Perturbations_down","4":"https://clue.io/","5":"3644"},{"1":"12638","2":"300","3":"L1000_Kinase_and_GPCR_Perturbations_up","4":"https://clue.io/","5":"3644"},{"1":"8973","2":"64","3":"Reactome_2016","4":"http://www.reactome.org/download/index.html","5":"1530"},{"1":"7010","2":"87","3":"KEGG_2016","4":"http://www.kegg.jp/kegg/download/","5":"293"},{"1":"5966","2":"51","3":"WikiPathways_2016","4":"http://www.wikipathways.org/index.php/Download_Pathways","5":"437"},{"1":"15562","2":"887","3":"ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X","4":"","5":"104"},{"1":"17850","2":"300","3":"Kinase_Perturbations_from_GEO_down","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"285"},{"1":"17660","2":"300","3":"Kinase_Perturbations_from_GEO_up","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"285"},{"1":"1348","2":"19","3":"BioCarta_2016","4":"http://cgap.nci.nih.gov/Pathways/BioCarta_Pathways","5":"237"},{"1":"934","2":"13","3":"HumanCyc_2016","4":"http://humancyc.org/","5":"152"},{"1":"2541","2":"39","3":"NCI-Nature_2016","4":"http://pid.nci.nih.gov/","5":"209"},{"1":"2041","2":"42","3":"Panther_2016","4":"http://www.pantherdb.org/pathway/","5":"112"},{"1":"5209","2":"300","3":"DrugMatrix","4":"https://ntp.niehs.nih.gov/drugmatrix/","5":"7876"},{"1":"49238","2":"1550","3":"ChEA_2016","4":"http://amp.pharm.mssm.edu/Enrichr","5":"645"},{"1":"2243","2":"19","3":"huMAP","4":"http://proteincomplexes.org/","5":"995"},{"1":"19586","2":"545","3":"Jensen_TISSUES","4":"http://tissues.jensenlab.org/","5":"1842"},{"1":"22440","2":"505","3":"RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"1302"},{"1":"8184","2":"24","3":"MGI_Mammalian_Phenotype_2017","4":"http://www.informatics.jax.org/","5":"5231"},{"1":"18329","2":"161","3":"Jensen_COMPARTMENTS","4":"http://compartments.jensenlab.org/","5":"2283"},{"1":"15755","2":"28","3":"Jensen_DISEASES","4":"http://diseases.jensenlab.org/","5":"1811"},{"1":"10271","2":"22","3":"BioPlex_2017","4":"http://bioplex.hms.harvard.edu/","5":"3915"},{"1":"10427","2":"38","3":"GO_Cellular_Component_2017","4":"http://www.geneontology.org/","5":"636"},{"1":"10601","2":"25","3":"GO_Molecular_Function_2017","4":"http://www.geneontology.org/","5":"972"},{"1":"13822","2":"21","3":"GO_Biological_Process_2017","4":"http://www.geneontology.org/","5":"3166"},{"1":"8002","2":"143","3":"GO_Cellular_Component_2017b","4":"http://www.geneontology.org/","5":"816"},{"1":"10089","2":"45","3":"GO_Molecular_Function_2017b","4":"http://www.geneontology.org/","5":"3271"},{"1":"13247","2":"49","3":"GO_Biological_Process_2017b","4":"http://www.geneontology.org/","5":"10125"},{"1":"21809","2":"2316","3":"ARCHS4_Tissues","4":"http://amp.pharm.mssm.edu/archs4","5":"108"},{"1":"23601","2":"2395","3":"ARCHS4_Cell-lines","4":"http://amp.pharm.mssm.edu/archs4","5":"125"},{"1":"20883","2":"299","3":"ARCHS4_IDG_Coexp","4":"http://amp.pharm.mssm.edu/archs4","5":"352"},{"1":"19612","2":"299","3":"ARCHS4_Kinases_Coexp","4":"http://amp.pharm.mssm.edu/archs4","5":"498"},{"1":"25983","2":"299","3":"ARCHS4_TFs_Coexp","4":"http://amp.pharm.mssm.edu/archs4","5":"1724"},{"1":"19500","2":"137","3":"SysMyo_Muscle_Gene_Sets","4":"http://sys-myo.rhcloud.com/","5":"1135"},{"1":"14893","2":"128","3":"miRTarBase_2017","4":"http://mirtarbase.mbc.nctu.edu.tw/","5":"3240"},{"1":"17598","2":"1208","3":"TargetScan_microRNA_2017","4":"http://www.targetscan.org/","5":"683"},{"1":"5902","2":"109","3":"Enrichr_Libraries_Most_Popular_Genes","4":"http://amp.pharm.mssm.edu/Enrichr","5":"121"},{"1":"12486","2":"299","3":"Enrichr_Submissions_TF-Gene_Coocurrence","4":"http://amp.pharm.mssm.edu/Enrichr","5":"1722"},{"1":"1073","2":"100","3":"Data_Acquisition_Method_Most_Popular_Genes","4":"http://amp.pharm.mssm.edu/Enrichr","5":"12"},{"1":"19513","2":"117","3":"DSigDB","4":"http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/","5":"4026"},{"1":"14433","2":"36","3":"GO_Biological_Process_2018","4":"http://www.geneontology.org/","5":"5103"},{"1":"8655","2":"61","3":"GO_Cellular_Component_2018","4":"http://www.geneontology.org/","5":"446"},{"1":"11459","2":"39","3":"GO_Molecular_Function_2018","4":"http://www.geneontology.org/","5":"1151"},{"1":"19741","2":"270","3":"TF_Perturbations_Followed_by_Expression","4":"http://www.ncbi.nlm.nih.gov/geo/","5":"1958"},{"1":"27360","2":"802","3":"Chromosome_Location_hg19","4":"http://hgdownload.cse.ucsc.edu/downloads.html","5":"36"},{"1":"13072","2":"26","3":"NIH_Funded_PIs_2017_Human_GeneRIF","4":"https://www.ncbi.nlm.nih.gov/pubmed/","5":"5687"},{"1":"13464","2":"45","3":"NIH_Funded_PIs_2017_Human_AutoRIF","4":"https://www.ncbi.nlm.nih.gov/pubmed/","5":"12558"},{"1":"13787","2":"200","3":"Rare_Diseases_AutoRIF_ARCHS4_Predictions","4":"https://amp.pharm.mssm.edu/geneshot/","5":"3725"},{"1":"13929","2":"200","3":"Rare_Diseases_GeneRIF_ARCHS4_Predictions","4":"https://www.ncbi.nlm.nih.gov/gene/about-generif","5":"2244"},{"1":"16964","2":"200","3":"NIH_Funded_PIs_2017_AutoRIF_ARCHS4_Predictions","4":"https://www.ncbi.nlm.nih.gov/pubmed/","5":"12558"},{"1":"17258","2":"200","3":"NIH_Funded_PIs_2017_GeneRIF_ARCHS4_Predictions","4":"https://www.ncbi.nlm.nih.gov/pubmed/","5":"5684"},{"1":"10352","2":"58","3":"Rare_Diseases_GeneRIF_Gene_Lists","4":"https://www.ncbi.nlm.nih.gov/gene/about-generif","5":"2244"},{"1":"10471","2":"76","3":"Rare_Diseases_AutoRIF_Gene_Lists","4":"https://amp.pharm.mssm.edu/geneshot/","5":"3725"},{"1":"12419","2":"491","3":"SubCell_BarCode","4":"http://www.subcellbarcode.org/","5":"104"},{"1":"19378","2":"37","3":"GWAS_Catalog_2019","4":"https://www.ebi.ac.uk/gwas","5":"1737"},{"1":"6201","2":"45","3":"WikiPathways_2019_Human","4":"https://www.wikipathways.org/","5":"472"},{"1":"4558","2":"54","3":"WikiPathways_2019_Mouse","4":"https://www.wikipathways.org/","5":"176"},{"1":"3264","2":"22","3":"TRRUST_Transcription_Factors_2019","4":"https://www.grnpedia.org/trrust/","5":"571"},{"1":"7802","2":"92","3":"KEGG_2019_Human","4":"https://www.kegg.jp/","5":"308"},{"1":"8551","2":"98","3":"KEGG_2019_Mouse","4":"https://www.kegg.jp/","5":"303"},{"1":"12444","2":"23","3":"InterPro_Domains_2019","4":"https://www.ebi.ac.uk/interpro/","5":"1071"},{"1":"9000","2":"20","3":"Pfam_Domains_2019","4":"https://pfam.xfam.org/","5":"608"},{"1":"7744","2":"363","3":"DepMap_WG_CRISPR_Screens_Broad_CellLines_2019","4":"https://depmap.org/","5":"558"},{"1":"6204","2":"387","3":"DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019","4":"https://depmap.org/","5":"325"},{"1":"13420","2":"32","3":"MGI_Mammalian_Phenotype_Level_4_2019","4":"http://www.informatics.jax.org/","5":"5261"},{"1":"14148","2":"122","3":"UK_Biobank_GWAS_v1","4":"https://www.ukbiobank.ac.uk/tag/gwas/","5":"857"},{"1":"9813","2":"49","3":"BioPlanet_2019","4":"https://tripod.nih.gov/bioplanet/","5":"1510"},{"1":"1397","2":"13","3":"ClinVar_2019","4":"https://www.ncbi.nlm.nih.gov/clinvar/","5":"182"},{"1":"9116","2":"22","3":"PheWeb_2019","4":"http://pheweb.sph.umich.edu/","5":"1161"},{"1":"17464","2":"63","3":"DisGeNET","4":"https://www.disgenet.org","5":"9828"},{"1":"394","2":"73","3":"HMS_LINCS_KinomeScan","4":"http://lincs.hms.harvard.edu/kinomescan/","5":"148"},{"1":"11851","2":"586","3":"CCLE_Proteomics_2020","4":"https://portals.broadinstitute.org/ccle","5":"378"},{"1":"8189","2":"421","3":"ProteomicsDB_2020","4":"https://www.proteomicsdb.org/","5":"913"},{"1":"18704","2":"100","3":"lncHUB_lncRNA_Co-Expression","4":"https://amp.pharm.mssm.edu/lnchub/","5":"3729"},{"1":"5605","2":"39","3":"Virus-Host_PPI_P-HIPSTer_2020","4":"http://phipster.org/","5":"6715"},{"1":"5718","2":"31","3":"Elsevier_Pathway_Collection","4":"http://www.transgene.ru/disease-pathways/","5":"1721"},{"1":"14156","2":"40","3":"Table_Mining_of_CRISPR_Studies","4":"","5":"802"},{"1":"16979","2":"295","3":"COVID-19_Related_Gene_Sets","4":"https://amp.pharm.mssm.edu/covid19","5":"205"},{"1":"4383","2":"146","3":"MSigDB_Hallmark_2020","4":"https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp","5":"50"},{"1":"54974","2":"483","3":"Enrichr_Users_Contributed_Lists_2020","4":"https://maayanlab.cloud/Enrichr","5":"1482"},{"1":"12118","2":"448","3":"TG_GATES_2020","4":"https://toxico.nibiohn.go.jp/english/","5":"1190"},{"1":"12361","2":"124","3":"Allen_Brain_Atlas_10x_scRNA_2021","4":"https://portal.brain-map.org/","5":"766"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# Perform enrichment
top_DGE <- DGE_cell_selection$Covid[(DGE_cell_selection$Covid$p.value < 0.01) & (abs(DGE_cell_selection$Covid$summary.logFC) > 
    0.25), ]

enrich_results <- enrichr(genes = rownames(top_DGE), databases = "GO_Biological_Process_2017b")[[1]]
```

```
## Uploading data to Enrichr... Done.
##   Querying GO_Biological_Process_2017b... Done.
## Parsing results... Done.
```


Some databases of interest:

* `GO_Biological_Process_2017b`
* `KEGG_2019_Human`
* `KEGG_2019_Mouse`
* `WikiPathways_2019_Human`
* `WikiPathways_2019_Mouse`

You visualize your results using a simple barplot, for example:


```r
par(mar = c(3, 25, 2, 1))
barplot(height = -log10(enrich_results$P.value)[10:1], names.arg = enrich_results$Term[10:1], 
    horiz = TRUE, las = 1, border = FALSE, cex.names = 0.6)
abline(v = c(-log10(0.05)), lty = 2)
abline(v = 0, lty = 1)
```

![](scater_05_dge_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

Gene Set Enrichment Analysis (GSEA)

Besides the enrichment using hypergeometric test, we can also perform gene set enrichment analysis (GSEA), which scores ranked genes list (usually based on fold changes) and computes permutation test to check if a particular gene set is more present in the Up-regulated genes, amongthe DOWN_regulated genes or not differentially regulated.


```r
# Create a gene rank based on the gene expression fold change
gene_rank <- setNames(DGE_cell_selection$Covid$summary.logFC, casefold(rownames(DGE_cell_selection$Covid), 
    upper = T))
```

 Once our list of genes are sorted, we can proceed with the enrichment itself. We can use the package to get gene set from the Molecular Signature Database (MSigDB) and select KEGG pathways as an example.


```r
# install.packages('msigdbr')
library(msigdbr)

# Download gene sets
msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

# List available gene sets
unique(msigdbgmt$gs_subcat)
```

```
##  [1] "MIR:MIR_Legacy"  "TFT:TFT_Legacy"  "CGP"             "TFT:GTRD"       
##  [5] ""                "CP:BIOCARTA"     "CGN"             "GO:MF"          
##  [9] "GO:BP"           "GO:CC"           "HPO"             "CP:KEGG"        
## [13] "MIR:MIRDB"       "CM"              "CP"              "CP:PID"         
## [17] "CP:REACTOME"     "CP:WIKIPATHWAYS"
```

```r
# Subset which gene set you want to use.
msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat == "CP:WIKIPATHWAYS", ]
gmt <- lapply(unique(msigdbgmt_subset$gs_name), function(x) {
    msigdbgmt_subset[msigdbgmt_subset$gs_name == x, "gene_symbol"]
})
names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name, "_", msigdbgmt_subset$gs_exact_source))
```

 Next, we will be using the GSEA. This will result in a table containing information for several pathways. We can then sort and filter those pathways to visualize only the top ones. You can select/filter them by either `p-value` or normalized enrichemnet score (`NES`).


```r
library(fgsea)

# Perform enrichemnt analysis
fgseaRes <- fgsea(pathways = gmt, stats = gene_rank, minSize = 15, maxSize = 500)
fgseaRes <- fgseaRes[order(fgseaRes$NES, decreasing = T), ]

# Filter the results table to show only the top 10 UP or DOWN regulated processes
# (optional)
top10_UP <- fgseaRes$pathway[1:10]

# Nice summary table (shown as a plot)
plotGseaTable(gmt[top10_UP], gene_rank, fgseaRes, gseaParam = 0.5)
```

![](scater_05_dge_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

Which KEGG pathways are upregulated in this cluster?
Which KEGG pathways are dowregulated in this cluster?
Change the pathway source to another gene set (e.g. "CP:WIKIPATHWAYS" or "CP:REACTOME" or "CP:BIOCARTA" or "GO:BP") and check the if you get simmilar results?
</div>

Finally, lets save the integrated data for further analysis.



```r
saveRDS(sce, "data/results/covid_qc_dr_int_cl_dge.rds")
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
##  [1] fgsea_1.16.0                msigdbr_7.2.1              
##  [3] enrichR_2.1                 dplyr_1.0.3                
##  [5] igraph_1.2.6                pheatmap_1.0.12            
##  [7] reticulate_1.18             harmony_1.0                
##  [9] Rcpp_1.0.5                  venn_1.9                   
## [11] umap_0.2.7.0                rafalib_1.0.0              
## [13] scDblFinder_1.4.0           org.Hs.eg.db_3.12.0        
## [15] AnnotationDbi_1.52.0        cowplot_1.1.1              
## [17] scran_1.18.0                scater_1.18.0              
## [19] ggplot2_3.3.3               SingleCellExperiment_1.12.0
## [21] SummarizedExperiment_1.20.0 Biobase_2.50.0             
## [23] GenomicRanges_1.42.0        GenomeInfoDb_1.26.0        
## [25] IRanges_2.24.0              S4Vectors_0.28.0           
## [27] BiocGenerics_0.36.0         MatrixGenerics_1.2.0       
## [29] matrixStats_0.57.0          RJSONIO_1.3-1.4            
## [31] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] tidyselect_1.1.0          RSQLite_2.2.1            
##   [3] htmlwidgets_1.5.3         grid_4.0.3               
##   [5] BiocParallel_1.24.0       Rtsne_0.15               
##   [7] munsell_0.5.0             codetools_0.2-18         
##   [9] ica_1.0-2                 statmod_1.4.35           
##  [11] xgboost_1.3.0.1           future_1.21.0            
##  [13] miniUI_0.1.1.1            withr_2.3.0              
##  [15] batchelor_1.6.0           colorspace_2.0-0         
##  [17] knitr_1.30                Seurat_3.2.3             
##  [19] ROCR_1.0-11               tensor_1.5               
##  [21] listenv_0.8.0             labeling_0.4.2           
##  [23] GenomeInfoDbData_1.2.4    polyclip_1.10-0          
##  [25] bit64_4.0.5               farver_2.0.3             
##  [27] parallelly_1.23.0         vctrs_0.3.6              
##  [29] generics_0.1.0            xfun_0.19                
##  [31] R6_2.5.0                  ggbeeswarm_0.6.0         
##  [33] rsvd_1.0.3                locfit_1.5-9.4           
##  [35] hdf5r_1.3.3               bitops_1.0-6             
##  [37] spatstat.utils_1.17-0     DelayedArray_0.16.0      
##  [39] assertthat_0.2.1          promises_1.1.1           
##  [41] scales_1.1.1              beeswarm_0.2.3           
##  [43] gtable_0.3.0              beachmat_2.6.0           
##  [45] globals_0.14.0            goftest_1.2-2            
##  [47] rlang_0.4.10              splines_4.0.3            
##  [49] lazyeval_0.2.2            yaml_2.2.1               
##  [51] reshape2_1.4.4            abind_1.4-5              
##  [53] httpuv_1.5.4              tools_4.0.3              
##  [55] ellipsis_0.3.1            RColorBrewer_1.1-2       
##  [57] ggridges_0.5.2            plyr_1.8.6               
##  [59] sparseMatrixStats_1.2.0   zlibbioc_1.36.0          
##  [61] purrr_0.3.4               RCurl_1.98-1.2           
##  [63] rpart_4.1-15              openssl_1.4.3            
##  [65] deldir_0.2-3              pbapply_1.4-3            
##  [67] viridis_0.5.1             zoo_1.8-8                
##  [69] ggrepel_0.9.0             cluster_2.1.0            
##  [71] magrittr_2.0.1            data.table_1.13.6        
##  [73] RSpectra_0.16-0           scattermore_0.7          
##  [75] ResidualMatrix_1.0.0      lmtest_0.9-38            
##  [77] RANN_2.6.1                fitdistrplus_1.1-3       
##  [79] patchwork_1.1.1           mime_0.9                 
##  [81] evaluate_0.14             xtable_1.8-4             
##  [83] gridExtra_2.3             compiler_4.0.3           
##  [85] tibble_3.0.4              KernSmooth_2.23-18       
##  [87] crayon_1.3.4              htmltools_0.5.0          
##  [89] mgcv_1.8-33               later_1.1.0.1            
##  [91] tidyr_1.1.2               DBI_1.1.0                
##  [93] formatR_1.7               MASS_7.3-53              
##  [95] Matrix_1.3-0              getopt_1.20.3            
##  [97] pkgconfig_2.0.3           plotly_4.9.2.2           
##  [99] scuttle_1.0.0             vipor_0.4.5              
## [101] admisc_0.11               dqrng_0.2.1              
## [103] XVector_0.30.0            stringr_1.4.0            
## [105] digest_0.6.27             sctransform_0.3.2        
## [107] RcppAnnoy_0.0.18          spatstat.data_1.7-0      
## [109] fastmatch_1.1-0           rmarkdown_2.6            
## [111] leiden_0.3.6              uwot_0.1.10              
## [113] edgeR_3.32.0              DelayedMatrixStats_1.12.0
## [115] curl_4.3                  shiny_1.5.0              
## [117] rjson_0.2.20              lifecycle_0.2.0          
## [119] nlme_3.1-151              jsonlite_1.7.2           
## [121] BiocNeighbors_1.8.0       viridisLite_0.3.0        
## [123] askpass_1.1               limma_3.46.0             
## [125] pillar_1.4.7              lattice_0.20-41          
## [127] fastmap_1.0.1             httr_1.4.2               
## [129] survival_3.2-7            glue_1.4.2               
## [131] spatstat_1.64-1           png_0.1-7                
## [133] bluster_1.0.0             bit_4.0.4                
## [135] stringi_1.5.3             blob_1.2.1               
## [137] BiocSingular_1.6.0        memoise_1.1.0            
## [139] irlba_2.3.3               future.apply_1.7.0
```



















