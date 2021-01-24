---
title: "Scater/Scran:: Spatial Transcriptomics"
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

# Spatial transcriptomics
***


Spatial transcriptomic data with the Visium platform is in many ways similar to scRNAseq data. It contains UMI counts for 5-20 cells instead of single cells, but is still quite sparse in the same way as scRNAseq data is, but with the additional information about spatial location in the tissue. 

Here we will first run quality control in a similar manner to scRNAseq data, then QC filtering, dimensionality reduction, integration and clustering. Then we will use scRNAseq data from mouse cortex to run LabelTransfer to predict celltypes in the Visium spots. 

We will use two **Visium** spatial transcriptomics dataset of the mouse brain (Sagittal), which are publicly available from the [10x genomics website](https://support.10xgenomics.com/spatial-gene-expression/datasets/). Note, that these dataset have already been filtered for spots that does not overlap with the tissue.

### Load packages


```bash
mkdir -p data/visium/Posterior
mkdir -p data/visium/Anterior

cd data/visium/Posterior

curl -o V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.tar.gz https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior/V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.tar.gz
tar xvzf V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.tar.gz

curl -o V1_Mouse_Brain_Sagittal_Posterior_spatial.tar.gz https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Posterior/V1_Mouse_Brain_Sagittal_Posterior_spatial.tar.gz
tar xvzf V1_Mouse_Brain_Sagittal_Posterior_spatial.tar.gz
rm *.tar.gz

cd ../Anterior

curl -o V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.tar.gz https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.tar.gz
tar xvzf V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.tar.gz

curl -o V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz
tar xvzf V1_Mouse_Brain_Sagittal_Anterior_spatial.tar.gz
rm *.tar.gz
cd ..

```

```
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0 54.8M    0  270k    0     0  1557k      0  0:00:36 --:--:--  0:00:36 1548k  9 54.8M    9 5054k    0     0  4298k      0  0:00:13  0:00:01  0:00:12 4294k 17 54.8M   17 9726k    0     0  4474k      0  0:00:12  0:00:02  0:00:10 4472k 26 54.8M   26 14.4M    0     0  4655k      0  0:00:12  0:00:03  0:00:09 4655k 33 54.8M   33 18.5M    0     0  4541k      0  0:00:12  0:00:04  0:00:08 4540k 41 54.8M   41 22.7M    0     0  4491k      0  0:00:12  0:00:05  0:00:07 4593k 50 54.8M   50 27.5M    0     0  4567k      0  0:00:12  0:00:06  0:00:06 4630k 58 54.8M   58 31.9M    0     0  4563k      0  0:00:12  0:00:07  0:00:05 4602k 66 54.8M   66 36.6M    0     0  4588k      0  0:00:12  0:00:08  0:00:04 4546k 75 54.8M   75 41.4M    0     0  4619k      0  0:00:12  0:00:09  0:00:03 4685k 81 54.8M   81 44.8M    0     0  4515k      0  0:00:12  0:00:10  0:00:02 4539k 89 54.8M   89 49.1M    0     0  4503k      0  0:00:12  0:00:11  0:00:01 4423k 98 54.8M   98 53.7M    0     0  4521k      0  0:00:12  0:00:12 --:--:-- 4460k100 54.8M  100 54.8M    0     0  4535k      0  0:00:12  0:00:12 --:--:-- 4432k
## filtered_feature_bc_matrix/
## filtered_feature_bc_matrix/barcodes.tsv.gz
## filtered_feature_bc_matrix/matrix.mtx.gz
## filtered_feature_bc_matrix/features.tsv.gz
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  7 9480k    7  754k    0     0  2328k      0  0:00:04 --:--:--  0:00:04 2328k 59 9480k   59 5679k    0     0  4296k      0  0:00:02  0:00:01  0:00:01 4293k100 9480k  100 9480k    0     0  4436k      0  0:00:02  0:00:02 --:--:-- 4438k
## spatial/
## spatial/tissue_positions_list.csv
## spatial/tissue_hires_image.png
## spatial/scalefactors_json.json
## spatial/aligned_fiducials.jpg
## spatial/detected_tissue_image.jpg
## spatial/tissue_lowres_image.png
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  5 54.8M    5 3006k    0     0  2965k      0  0:00:18  0:00:01  0:00:17 2968k 13 54.8M   13 7438k    0     0  3749k      0  0:00:14  0:00:01  0:00:13 3747k 21 54.8M   21 11.7M    0     0  4041k      0  0:00:13  0:00:02  0:00:11 4039k 30 54.8M   30 16.8M    0     0  4321k      0  0:00:12  0:00:03  0:00:09 4321k 39 54.8M   39 21.8M    0     0  4482k      0  0:00:12  0:00:04  0:00:08 4481k 48 54.8M   48 26.7M    0     0  4570k      0  0:00:12  0:00:05  0:00:07 4897k 57 54.8M   57 31.6M    0     0  4631k      0  0:00:12  0:00:06  0:00:06 4981k 66 54.8M   66 36.2M    0     0  4651k      0  0:00:12  0:00:07  0:00:05 5014k 75 54.8M   75 41.1M    0     0  4685k      0  0:00:11  0:00:08  0:00:03 4975k 84 54.8M   84 46.2M    0     0  4739k      0  0:00:11  0:00:09  0:00:02 4995k 92 54.8M   92 50.5M    0     0  4710k      0  0:00:11  0:00:10  0:00:01 4877k100 54.8M  100 54.8M    0     0  4709k      0  0:00:11  0:00:11 --:--:-- 4820k
## filtered_feature_bc_matrix/
## filtered_feature_bc_matrix/barcodes.tsv.gz
## filtered_feature_bc_matrix/features.tsv.gz
## filtered_feature_bc_matrix/matrix.mtx.gz
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0 24 9017k   24 2174k    0     0  3446k      0  0:00:02 --:--:--  0:00:02 3446k 81 9017k   81 7342k    0     0  4480k      0  0:00:02  0:00:01  0:00:01 4480k100 9017k  100 9017k    0     0  4612k      0  0:00:01  0:00:01 --:--:-- 4612k
## spatial/
## spatial/tissue_positions_list.csv
## spatial/tissue_hires_image.png
## spatial/scalefactors_json.json
## spatial/detected_tissue_image.jpg
## spatial/tissue_lowres_image.png
## spatial/aligned_fiducials.jpg
```



```r
BiocManager::install("DropletUtils", update = F)
devtools::install_github("RachelQueen1/Spaniel", ref = "Development", upgrade = F, 
    dependencies = F)

library(Spaniel)
library(biomaRt)

suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(scran))
suppressPackageStartupMessages(require(SingleR))
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(cowplot))
```


### Load ST data

We can first load and merge the objects into one SCE object.


```r
sce.a <- Spaniel::createVisiumSCE(tenXDir = "data/visium/Anterior", resolution = "Low")
sce.p <- Spaniel::createVisiumSCE(tenXDir = "data/visium/Posterior", resolution = "Low")

sce <- cbind(sce.a, sce.p)
sce$Sample <- sub(".*[/]", "", sub("/filtered_feature_bc_matrix", "", sce$Sample))

lll <- list(sce.a, sce.p)
lll <- lapply(lll, function(x) x@metadata)
names(lll) <- c("Anterior", "Posterior")
sce@metadata <- lll
```

We can further convert the gene ensembl IDs to gene names.


```r
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
annot <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
    mart = mart, useCache = F)

gene_names <- as.character(annot[match(rownames(sce), annot[, "ensembl_gene_id"]), 
    "external_gene_name"])
gene_names[is.na(gene_names)] <- ""

sce <- sce[gene_names != "", ]
rownames(sce) <- gene_names[gene_names != ""]
dim(sce)
```

```
## [1] 32275  6050
```


##Quality control
***

Similar to scRNAseq we use statistics on number of counts, number of features and percent mitochondria for quality control. 

Now the counts and feature counts are calculated on the Spatial assay, so they are named  "nCount_Spatial" and "nFeature_Spatial".


```r
# Mitochondrial genes
mito_genes <- rownames(sce)[grep("^mt-", rownames(sce))]

# Ribosomal genes
ribo_genes <- rownames(sce)[grep("^Rp[sl]", rownames(sce))]

# Hemoglobin genes - includes all genes starting with HB except HBP.
hb_genes <- rownames(sce)[grep("^Hb[^(p)]", rownames(sce))]

sce <- addPerCellQC(sce, flatten = T, subsets = list(mt = mito_genes, hb = hb_genes, 
    ribo = ribo_genes))

head(colData(sce))
```

```
## DataFrame with 6 rows and 32 columns
##        Sample            Barcode   Section    Spot_Y    Spot_X   Image_Y
##   <character>        <character> <integer> <integer> <integer> <integer>
## 1    Anterior AAACAAGTATCTCCCA-1         1        50       102      7474
## 2    Anterior AAACACCAATAACTGC-1         1        59        19      8552
## 3    Anterior AAACAGAGCGACTCCT-1         1        14        94      3163
## 4    Anterior AAACAGCTTTCAGAAG-1         1        43         9      6636
## 5    Anterior AAACAGGGTCTATATT-1         1        47        13      7115
## 6    Anterior AAACATGGTGAGAGGA-1         1        62         0      8912
##     Image_X       pixel_x       pixel_y       sum  detected   percent_top_50
##   <integer>     <numeric>     <numeric> <numeric> <integer>        <numeric>
## 1      8500   438.8984605 214.079165438     13991      4462 25.5235508541205
## 2      2788 143.958695044 158.416513624     39797      8126 20.7603588210167
## 3      7950  410.49914835 436.678137581     29951      6526 28.4230910487129
## 4      2100   108.4337373 257.349390132     42333      8190 15.6591784187277
## 5      2375 122.633393375 232.616171005     35700      8090 13.1960784313725
## 6      1480   76.41996724 139.827872944     22148      6518 20.3946180242008
##    percent_top_100  percent_top_200  percent_top_500     total       sum
##          <numeric>        <numeric>        <numeric> <numeric> <numeric>
## 1 31.7489814880995  40.132942605961 54.1133585876635     13991     13963
## 2 26.1326230620399  33.768877051034 47.1442571048069     39797     39764
## 3 33.6115655570766 40.7498914894327 53.1935494641247     29951     29919
## 4 21.7608012661517 30.0309451255522  44.157040606619     42333     42293
## 5 18.9859943977591 27.1652661064426 40.9635854341737     35700     35686
## 6  26.363554271266 33.9804948528084 47.0200469568358     22148     22115
##    detected   percent_top_50  percent_top_100  percent_top_200  percent_top_500
##   <integer>        <numeric>        <numeric>        <numeric>        <numeric>
## 1      4461 25.5317625152188 31.7123827257753 40.0773472749409 54.0571510420397
## 2      8124 20.7775877678302 26.1543104315461 33.7792978573584 47.1356000402374
## 3      6525 28.4534910926167 33.6308031685551 40.7433403522845 53.1735686353153
## 4      8187  15.673988603315   21.78138226184  30.035703307876  44.142056605112
## 5      8089 13.2012553942723 18.9934428067029 27.1759233312784 40.9740514487474
## 6      6517 20.4205290526792 26.3486321501244 33.9407641872033  46.972643002487
##   subsets_mt_sum subsets_mt_detected subsets_mt_percent subsets_hb_sum
##        <numeric>           <integer>          <numeric>      <numeric>
## 1           1521                  12   10.8930745541789             60
## 2           3977                  12   10.0015089025249            831
## 3           4265                  12   14.2551555867509            111
## 4           2870                  12   6.78599295391672            117
## 5           1831                  13   5.13086364400605             73
## 6           2390                  12   10.8071444720778            134
##   subsets_hb_detected subsets_hb_percent subsets_ribo_sum subsets_ribo_detected
##             <integer>          <numeric>        <numeric>             <integer>
## 1                   4  0.429707083005085              826                    85
## 2                   6    2.0898299969822             2199                    89
## 3                   5  0.371001704602427             1663                    88
## 4                   5   0.27664152460218             3129                    88
## 5                   5  0.204562013114387             2653                    90
## 6                   5  0.605923581279674             1478                    84
##   subsets_ribo_percent     total
##              <numeric> <numeric>
## 1     5.91563417603667     13963
## 2     5.53012775374711     39764
## 3     5.55834085363816     29919
## 4     7.39838744000189     42293
## 5      7.4342879560612     35686
## 6     6.68324666515939     22115
```

```r
plot_grid(plotColData(sce, y = "detected", x = "Sample", colour_by = "Sample"), plotColData(sce, 
    y = "total", x = "Sample", colour_by = "Sample"), plotColData(sce, y = "subsets_mt_percent", 
    x = "Sample", colour_by = "Sample"), plotColData(sce, y = "subsets_ribo_percent", 
    x = "Sample", colour_by = "Sample"), plotColData(sce, y = "subsets_hb_percent", 
    x = "Sample", colour_by = "Sample"), ncol = 3)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

We can also plot the same data onto the tissue section.


```r
samples <- c("Anterior", "Posterior")
to_plot <- c("detected", "total", "subsets_mt_percent", "subsets_ribo_percent", "subsets_hb_percent")

plist <- list()
n = 1
for (j in to_plot) {
    for (i in samples) {
        temp <- sce[, sce$Sample == i]
        temp@metadata <- temp@metadata[[i]]
        plist[[n]] <- spanielPlot(object = temp, plotType = "Cluster", clusterRes = j, 
            customTitle = j, techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
        n <- n + 1
    }
}

plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


As you can see, the spots with low number of counts/features and high mitochondrial content is mainly towards the edges of the tissue. It is quite likely that these regions are damaged tissue. You may also see regions within a tissue with low quality if you have tears or folds in your section. 

But remember, for some tissue types, the amount of genes expressed and proportion mitochondria may also be a biological features, so bear in mind what tissue you are working on and what these features mean.

### Filter

Select all spots with less than 25% mitocondrial reads, less than 20% hb-reads and 1000 detected genes. You must judge for yourself based on your knowledge of the tissue what are appropriate filtering criteria for your dataset.



```r
sce <- sce[, sce$detected > 500 & sce$subsets_mt_percent < 25 & sce$subsets_hb_percent < 
    20]
dim(sce)
```

```
## [1] 32275  5805
```

And replot onto tissue section:


```r
samples <- c("Anterior", "Posterior")
to_plot <- c("detected", "total", "subsets_mt_percent", "subsets_mt_percent", "subsets_hb_percent")

plist <- list()
n = 1
for (j in to_plot) {
    for (i in samples) {
        temp <- sce[, sce$Sample == i]
        temp@metadata <- temp@metadata[[i]]
        plist[[n]] <- spanielPlot(object = temp, plotType = "Cluster", clusterRes = j, 
            customTitle = j, techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
        n <- n + 1
    }
}

plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

### Top expressed genes
As for scRNAseq data, we will look at what the top expressed genes are.


```r
C = counts(sce)
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", 
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
rm(C)
```

As you can see, the mitochondrial genes are among the top expressed. Also the lncRNA gene Bc1 (brain cytoplasmic RNA 1). Also one hemoglobin gene.

### Filter genes
We will remove the Bc1 gene, hemoglobin genes (blood contamination) and the mitochondrial genes.


```r
dim(sce)
```

```
## [1] 32275  5805
```

```r
# Filter Bl1
sce <- sce[!grepl("Bc1", rownames(sce)), ]

# Filter Mitocondrial
sce <- sce[!grepl("^mt-", rownames(sce)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
sce <- sce[!grepl("^Hb.*-", rownames(sce)), ]

dim(sce)
```

```
## [1] 32253  5805
```

## Analysis
***


```r
sce <- computeSumFactors(sce, sizes = c(20, 40, 60, 80))
sce <- logNormCounts(sce)
```


Now we can plot gene expression of individual genes, the gene Hpca is a strong hippocampal marker and Ttr is a marker of the choroid plexus.


```r
samples <- c("Anterior", "Posterior")
to_plot <- c("Hpca", "Ttr")

plist <- list()
n = 1
for (j in to_plot) {
    for (i in samples) {
        temp <- sce[, sce$Sample == i]
        temp@metadata <- temp@metadata[[i]]
        plist[[n]] <- spanielPlot(object = temp, plotType = "Gene", gene = j, customTitle = j, 
            techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
        n <- n + 1
    }
}

plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-11-1.png)<!-- -->


### Dimensionality reduction and clustering
We can then now run dimensionality reduction and clustering using the same workflow as we use for scRNA-seq analysis. 

But make sure you run it on the `SCT` assay.


```r
var.out <- modelGeneVar(sce, method = "loess")
hvgs = getTopHVGs(var.out, n = 2000)
sce <- runPCA(sce, exprs_values = "logcounts", subset_row = hvgs, ncomponents = 50, 
    ntop = 100, scale = T)
g <- buildSNNGraph(sce, k = 5, use.dimred = "PCA")
sce$louvain_SNNk5 <- factor(igraph::cluster_louvain(g)$membership)
sce <- runUMAP(sce, dimred = "PCA", n_dimred = 50, ncomponents = 2, min_dist = 0.1, 
    spread = 0.3, metric = "correlation", name = "UMAP_on_PCA")
```

We can then plot clusters onto umap or onto the tissue section.


```r
samples <- c("Anterior", "Posterior")
to_plot <- c("louvain_SNNk5")

plist <- list()
n = 1
for (j in to_plot) {
    for (i in samples) {
        temp <- sce[, sce$Sample == i]
        temp@metadata <- temp@metadata[[i]]
        plist[[n]] <- spanielPlot(object = temp, plotType = "Cluster", clusterRes = j, 
            customTitle = j, techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
        n <- n + 1
    }
}

plist[[3]] <- plotReducedDim(sce, dimred = "UMAP_on_PCA", colour_by = "louvain_SNNk5")
plist[[4]] <- plotReducedDim(sce, dimred = "UMAP_on_PCA", colour_by = "Sample")

plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


## Integration

Quite often there are strong batch effects between different ST sections, so it may be a good idea to integrate the data across sections.

We will do a similar integration as in the Data Integration lab.


```r
mnn_out <- batchelor::fastMNN(sce, subset.row = hvgs, batch = factor(sce$Sample), 
    k = 20, d = 50)

reducedDim(sce, "MNN") <- reducedDim(mnn_out, "corrected")
rm(mnn_out)
gc()
```

```
##             used   (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells   7709695  411.8   11801654  630.3  11801654  630.3
## Vcells 187630638 1431.6  343344438 2619.6 343049786 2617.3
```


Then we run dimensionality reduction and clustering as before.


```r
g <- buildSNNGraph(sce, k = 5, use.dimred = "MNN")
sce$louvain_SNNk5 <- factor(igraph::cluster_louvain(g)$membership)
sce <- runUMAP(sce, dimred = "MNN", n_dimred = 50, ncomponents = 2, min_dist = 0.1, 
    spread = 0.3, metric = "correlation", name = "UMAP_on_MNN")
```


```r
samples <- c("Anterior", "Posterior")
to_plot <- c("louvain_SNNk5")

plist <- list()
n = 1
for (j in to_plot) {
    for (i in samples) {
        temp <- sce[, sce$Sample == i]
        temp@metadata <- temp@metadata[[i]]
        plist[[n]] <- spanielPlot(object = temp, plotType = "Cluster", clusterRes = j, 
            customTitle = j, techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
        n <- n + 1
    }
}

plist[[3]] <- plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = "louvain_SNNk5")
plist[[4]] <- plotReducedDim(sce, dimred = "UMAP_on_MNN", colour_by = "Sample")

plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

Do you see any differences between the integrated and non-integrated clusering? Judge for yourself, which of the clusterings do you think looks best? 
As a reference, you can compare to brain regions in the [Allen brain atlas](https://mouse.brain-map.org/experiment/thumbnails/100042147?image_type=atlas). 

## Identification of Spatially Variable Features

 There are two main workflows to identify molecular features that correlate with spatial location within a tissue. The first is to perform differential expression based on spatially distinct clusters, the other is to find features that are have spatial patterning without taking clusters or spatial annotation into account. 

First, we will do differential expression between clusters just as we did for the scRNAseq data before.


```r
# differential expression between cluster 4 and cluster 6
cell_selection <- sce[, sce$louvain_SNNk5 %in% c(4, 6)]
cell_selection$louvain_SNNk5 <- factor(cell_selection$louvain_SNNk5)

markers_genes <- scran::findMarkers(x = cell_selection, groups = cell_selection$louvain_SNNk5, 
    lfc = 0.25, pval.type = "all", direction = "up")

# List of dataFrames with the results for each cluster
top5_cell_selection <- lapply(names(markers_genes), function(x) {
    temp <- markers_genes[[x]][1:5, 1:2]
    temp$gene <- rownames(markers_genes[[x]])[1:5]
    temp$cluster <- x
    return(temp)
})
top5_cell_selection <- as_tibble(do.call(rbind, top5_cell_selection))
top5_cell_selection
```

<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":["p.value"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["gene"],"name":[3],"type":["chr"],"align":["left"]},{"label":["cluster"],"name":[4],"type":["chr"],"align":["left"]}],"data":[{"1":"2.103220e-22","2":"6.783515e-18","3":"Mbp","4":"4"},{"1":"2.600718e-20","2":"4.194048e-16","3":"Mobp","4":"4"},{"1":"8.586278e-19","2":"9.231107e-15","3":"Tshz2","4":"4"},{"1":"5.096770e-17","2":"4.109653e-13","3":"Tmsb10","4":"4"},{"1":"1.492235e-15","2":"9.625808e-12","3":"Plp1","4":"4"},{"1":"1.962206e-21","2":"6.328704e-17","3":"Wfs1","4":"6"},{"1":"8.839835e-17","2":"1.425556e-12","3":"Nptxr","4":"6"},{"1":"3.154485e-16","2":"3.391387e-12","3":"Calb1","4":"6"},{"1":"6.342810e-16","2":"5.114366e-12","3":"Lamp5","4":"6"},{"1":"1.340552e-14","2":"8.647365e-11","3":"Tmem215","4":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>

```r
# plot top markers
samples <- c("Anterior", "Posterior")
to_plot <- top5_cell_selection$gene[1:5]

plist <- list()
n = 1
for (j in to_plot) {
    for (i in samples) {
        temp <- sce[, sce$Sample == i]
        temp@metadata <- temp@metadata[[i]]
        plist[[n]] <- spanielPlot(object = temp, plotType = "Gene", gene = j, customTitle = j, 
            techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
        n <- n + 1
    }
}
plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


## Single cell data

We can also perform data integration between one scRNA-seq dataset and one spatial transcriptomics dataset. Such task is particularly useful because it allows us to transfer cell type labels to the Visium dataset, which were dentified from the scRNA-seq dataset. 

We will use a reference scRNA-seq dataset of ~14,000 adult mouse cortical cell taxonomy from the Allen Institute, generated with the SMART-Seq2 protocol.


First dowload the seurat data from: https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1 to folder `data/spatial/` with command:



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

For speed, and for a more fair comparison of the celltypes, we will subsample all celltypes to a maximum of 200 cells per class (`subclass`).


```r
allen_reference <- readRDS("data/spatial/allen_cortex.rds")
allen_reference_sce <- Seurat::as.SingleCellExperiment(allen_reference)

# check number of cells per subclass
allen_reference_sce$subclass <- sub("/", "_", sub(" ", "_", allen_reference_sce$subclass))
table(allen_reference_sce$subclass)
```

```
## 
##      Astro         CR       Endo    L2_3_IT         L4      L5_IT      L5_PT 
##        368          7         94        982       1401        880        544 
##      L6_CT      L6_IT        L6b      Lamp5 Macrophage      Meis2         NP 
##        960       1872        358       1122         51         45        362 
##      Oligo       Peri      Pvalb   Serpinf1        SMC       Sncg        Sst 
##         91         32       1337         27         55        125       1741 
##        Vip       VLMC 
##       1728         67
```

```r
# select 20 cells per subclass, fist set subclass ass active.ident
subset_cells <- lapply(unique(allen_reference_sce$subclass), function(x) {
    if (sum(allen_reference_sce$subclass == x) > 20) {
        temp <- sample(colnames(allen_reference_sce)[allen_reference_sce$subclass == 
            x], size = 20)
    } else {
        temp <- colnames(allen_reference_sce)[allen_reference_sce$subclass == x]
    }
})
allen_reference_sce <- allen_reference_sce[, unlist(subset_cells)]

# check again number of cells per subclass
table(allen_reference_sce$subclass)
```

```
## 
##      Astro         CR       Endo    L2_3_IT         L4      L5_IT      L5_PT 
##         20          7         20         20         20         20         20 
##      L6_CT      L6_IT        L6b      Lamp5 Macrophage      Meis2         NP 
##         20         20         20         20         20         20         20 
##      Oligo       Peri      Pvalb   Serpinf1        SMC       Sncg        Sst 
##         20         20         20         20         20         20         20 
##        Vip       VLMC 
##         20         20
```

Then run normalization and dimensionality reduction.


```r
allen_reference_sce <- computeSumFactors(allen_reference_sce, sizes = c(20, 40, 60, 
    80))
allen_reference_sce <- logNormCounts(allen_reference_sce)
allen.var.out <- modelGeneVar(allen_reference_sce, method = "loess")
allen.hvgs = getTopHVGs(allen.var.out, n = 2000)
```


### Subset ST for cortex
Since the scRNAseq dataset was generated from the mouse cortex, we will subset the visium dataset in order to select mainly the spots part of the cortex. Note that the integration can also be performed on the whole brain slice, but it would give rise to false positive cell type assignments and and therefore it should be interpreted with more care.


### Integrate with scRNAseq

Here, will use SingleR for prediciting which cell types are present in the dataset.

We can first select the anterior part as an example (to speed up predictions).


```r
sce.anterior <- sce[, sce$Sample == "Anterior"]
sce.anterior@metadata <- sce.anterior@metadata[["Anterior"]]
```

Next, we select the highly variable genes that are present in both datasets.


```r
# Find common highly variable genes
common_hvgs <- allen.hvgs[allen.hvgs %in% hvgs]

# Predict cell classes
pred.grun <- SingleR(test = sce.anterior[common_hvgs, ], ref = allen_reference_sce[common_hvgs, 
    ], labels = allen_reference_sce$subclass)

# Transfer the classes to the SCE object
sce.anterior$cell_prediction <- pred.grun$labels
sce.anterior@colData <- cbind(sce.anterior@colData, as.data.frame.matrix(table(list(1:ncol(sce.anterior), 
    sce.anterior$cell_prediction))))
```

Then we can plot the predicted cell populations back to tissue.


```r
# Plot cell predictions
spanielPlot(object = sce.anterior, plotType = "Cluster", clusterRes = "cell_prediction", 
    customTitle = "cell_prediction", techType = "Visium", ptSizeMax = 1, ptSizeMin = 0.1)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
plist <- list()
n = 1
for (i in c("L2_3_IT", "L4", "L5_IT", "L6_IT")) {
    plist[[n]] <- spanielPlot(object = sce.anterior, plotType = "Cluster", clusterRes = i, 
        customTitle = i, techType = "Visium", ptSize = 0.3, ptSizeMax = 1, ptSizeMin = 0.1)
    n <- n + 1
}
plot_grid(ncol = 2, plotlist = plist)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-22-2.png)<!-- -->


Keep in mind, that the scores are "just" prediction scores, and do not correspond to proportion of cells that are of a certain celltype or similar. It mainly tells you that gene expression in a certain spot is hihgly similar/dissimilar to gene expression of a celltype.

If we look at the scores, we see that some spots got really clear predictions by celltype, while others did not have high scores for any of the celltypes.


We can also plot the gene expression and add filters together, too:


```r
spanielPlot(object = sce.anterior, plotType = "Gene", gene = "Wfs1", showFilter = sce.anterior$L4, 
    customTitle = "", techType = "Visium", ptSize = 0, ptSizeMin = -0.3, ptSizeMax = 1)
```

![](scater_07_spatial_files/figure-html/unnamed-chunk-23-1.png)<!-- -->


### Session info


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
##  [1] cowplot_1.1.1               patchwork_1.1.1            
##  [3] scater_1.14.0               ggplot2_3.3.3              
##  [5] SingleR_1.0.0               scran_1.14.1               
##  [7] SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.0
##  [9] DelayedArray_0.12.0         BiocParallel_1.20.1        
## [11] matrixStats_0.57.0          Biobase_2.46.0             
## [13] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
## [15] IRanges_2.20.2              S4Vectors_0.24.4           
## [17] BiocGenerics_0.32.0         dplyr_1.0.3                
## [19] Matrix_1.3-2                biomaRt_2.42.1             
## [21] Spaniel_1.2.0               RJSONIO_1.3-1.4            
## [23] optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] R.utils_2.10.1                reticulate_1.18              
##   [3] tidyselect_1.1.0              RSQLite_2.2.2                
##   [5] AnnotationDbi_1.48.0          htmlwidgets_1.5.3            
##   [7] grid_3.6.1                    Rtsne_0.15                   
##   [9] DropletUtils_1.6.1            devtools_2.3.2               
##  [11] munsell_0.5.0                 codetools_0.2-18             
##  [13] ica_1.0-2                     statmod_1.4.35               
##  [15] future_1.21.0                 miniUI_0.1.1.1               
##  [17] batchelor_1.2.1               withr_2.4.0                  
##  [19] colorspace_2.0-0              knitr_1.30                   
##  [21] Seurat_3.2.3                  ROCR_1.0-11                  
##  [23] tensor_1.5                    listenv_0.8.0                
##  [25] labeling_0.4.2                GenomeInfoDbData_1.2.2       
##  [27] polyclip_1.10-0               farver_2.0.3                 
##  [29] bit64_4.0.5                   rhdf5_2.30.1                 
##  [31] rprojroot_2.0.2               parallelly_1.23.0            
##  [33] vctrs_0.3.6                   generics_0.1.0               
##  [35] xfun_0.20                     BiocFileCache_1.10.0         
##  [37] R6_2.5.0                      ggbeeswarm_0.6.0             
##  [39] rsvd_1.0.3                    locfit_1.5-9.4               
##  [41] bitops_1.0-6                  spatstat.utils_1.20-2        
##  [43] assertthat_0.2.1              promises_1.1.1               
##  [45] scales_1.1.1                  beeswarm_0.2.3               
##  [47] gtable_0.3.0                  globals_0.14.0               
##  [49] processx_3.4.5                goftest_1.2-2                
##  [51] rlang_0.4.10                  splines_3.6.1                
##  [53] lazyeval_0.2.2                BiocManager_1.30.10          
##  [55] yaml_2.2.1                    reshape2_1.4.4               
##  [57] abind_1.4-5                   httpuv_1.5.5                 
##  [59] tools_3.6.1                   usethis_2.0.0                
##  [61] ellipsis_0.3.1                RColorBrewer_1.1-2           
##  [63] sessioninfo_1.1.1             ggridges_0.5.3               
##  [65] Rcpp_1.0.6                    plyr_1.8.6                   
##  [67] progress_1.2.2                zlibbioc_1.32.0              
##  [69] purrr_0.3.4                   RCurl_1.98-1.2               
##  [71] ps_1.5.0                      prettyunits_1.1.1            
##  [73] rpart_4.1-15                  openssl_1.4.3                
##  [75] deldir_0.2-9                  viridis_0.5.1                
##  [77] pbapply_1.4-3                 zoo_1.8-8                    
##  [79] ggrepel_0.9.1                 cluster_2.1.0                
##  [81] fs_1.5.0                      magrittr_2.0.1               
##  [83] RSpectra_0.16-0               data.table_1.13.6            
##  [85] scattermore_0.7               lmtest_0.9-38                
##  [87] RANN_2.6.1                    fitdistrplus_1.1-3           
##  [89] pkgload_1.1.0                 hms_1.0.0                    
##  [91] mime_0.9                      evaluate_0.14                
##  [93] xtable_1.8-4                  XML_3.99-0.3                 
##  [95] gridExtra_2.3                 testthat_3.0.1               
##  [97] compiler_3.6.1                tibble_3.0.5                 
##  [99] KernSmooth_2.23-18            crayon_1.3.4                 
## [101] R.oo_1.24.0                   htmltools_0.5.1              
## [103] mgcv_1.8-33                   later_1.1.0.1                
## [105] tidyr_1.1.2                   DBI_1.1.1                    
## [107] ExperimentHub_1.12.0          formatR_1.7                  
## [109] dbplyr_2.0.0                  MASS_7.3-53                  
## [111] rappdirs_0.3.1                getopt_1.20.3                
## [113] cli_2.2.0                     R.methodsS3_1.8.1            
## [115] igraph_1.2.6                  pkgconfig_2.0.3              
## [117] plotly_4.9.3                  vipor_0.4.5                  
## [119] dqrng_0.2.1                   XVector_0.26.0               
## [121] stringr_1.4.0                 callr_3.5.1                  
## [123] digest_0.6.27                 sctransform_0.3.2            
## [125] RcppAnnoy_0.0.18              spatstat.data_1.7-0          
## [127] rmarkdown_2.6                 leiden_0.3.6                 
## [129] edgeR_3.28.1                  uwot_0.1.10                  
## [131] DelayedMatrixStats_1.8.0      curl_4.3                     
## [133] shiny_1.5.0                   lifecycle_0.2.0              
## [135] nlme_3.1-150                  jsonlite_1.7.2               
## [137] Rhdf5lib_1.8.0                BiocNeighbors_1.4.0          
## [139] limma_3.42.2                  desc_1.2.0                   
## [141] viridisLite_0.3.0             askpass_1.1                  
## [143] fansi_0.4.2                   pillar_1.4.7                 
## [145] lattice_0.20-41               fastmap_1.0.1                
## [147] httr_1.4.2                    pkgbuild_1.2.0               
## [149] survival_3.2-7                interactiveDisplayBase_1.24.0
## [151] glue_1.4.2                    remotes_2.2.0                
## [153] spatstat_1.64-1               png_0.1-7                    
## [155] BiocVersion_3.10.1            bit_4.0.4                    
## [157] HDF5Array_1.14.0              stringi_1.5.3                
## [159] blob_1.2.1                    AnnotationHub_2.18.0         
## [161] BiocSingular_1.2.0            memoise_1.1.0                
## [163] irlba_2.3.3                   future.apply_1.7.0
```
