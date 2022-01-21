---
title: "Scater/Scran:: Spatial Transcriptomics"
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
## 
  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
  3 54.8M    3 2141k    0     0  2288k      0  0:00:24 --:--:--  0:00:24 2288k
 16 54.8M   16 9405k    0     0  4864k      0  0:00:11  0:00:01  0:00:10 4863k
 32 54.8M   32 18.0M    0     0  6315k      0  0:00:08  0:00:02  0:00:06 6314k
 46 54.8M   46 25.5M    0     0  6650k      0  0:00:08  0:00:03  0:00:05 6649k
 57 54.8M   57 31.4M    0     0  6524k      0  0:00:08  0:00:04  0:00:04 6524k
 76 54.8M   76 41.9M    0     0  7205k      0  0:00:07  0:00:05  0:00:02 8122k
 92 54.8M   92 50.9M    0     0  7527k      0  0:00:07  0:00:06  0:00:01 8556k
100 54.8M  100 54.8M    0     0  7661k      0  0:00:07  0:00:07 --:--:-- 8560k
## x filtered_feature_bc_matrix/
## x filtered_feature_bc_matrix/barcodes.tsv.gz
## x filtered_feature_bc_matrix/matrix.mtx.gz
## x filtered_feature_bc_matrix/features.tsv.gz
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
## 
  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
  0 9480k    0  6363    0     0   7255      0  0:22:18 --:--:--  0:22:18  7255
 64 9480k   64 6093k    0     0  4034k      0  0:00:02  0:00:01  0:00:01 4032k
100 9480k  100 9480k    0     0  5003k      0  0:00:01  0:00:01 --:--:-- 5003k
## x spatial/
## x spatial/tissue_positions_list.csv
## x spatial/tissue_hires_image.png
## x spatial/scalefactors_json.json
## x spatial/aligned_fiducials.jpg
## x spatial/detected_tissue_image.jpg
## x spatial/tissue_lowres_image.png
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
## 
  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
  4 54.8M    4 2669k    0     0  4767k      0  0:00:11 --:--:--  0:00:11 4767k
 26 54.8M   26 14.3M    0     0  9398k      0  0:00:05  0:00:01  0:00:04 9393k
 51 54.8M   51 28.4M    0     0  10.8M      0  0:00:05  0:00:02  0:00:03 10.8M
 80 54.8M   80 44.2M    0     0  12.4M      0  0:00:04  0:00:03  0:00:01 12.4M
100 54.8M  100 54.8M    0     0  13.1M      0  0:00:04  0:00:04 --:--:-- 13.1M
## x filtered_feature_bc_matrix/
## x filtered_feature_bc_matrix/barcodes.tsv.gz
## x filtered_feature_bc_matrix/features.tsv.gz
## x filtered_feature_bc_matrix/matrix.mtx.gz
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
## 
  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
  2 9017k    2  223k    0     0   711k      0  0:00:12 --:--:--  0:00:12  710k
100 9017k  100 9017k    0     0  9811k      0 --:--:-- --:--:-- --:--:-- 9801k
## x spatial/
## x spatial/tissue_positions_list.csv
## x spatial/tissue_hires_image.png
## x spatial/scalefactors_json.json
## x spatial/detected_tissue_image.jpg
## x spatial/tissue_lowres_image.png
## x spatial/aligned_fiducials.jpg
```



```r
# BiocManager::install('DropletUtils',update = F)
devtools::install_github("RachelQueen1/Spaniel", ref = "Development", upgrade = F,
    dependencies = F)
```

```
##   
   checking for file ‘/private/var/folders/1s/j9ck5c_162s487xcprlxtmdh0000gp/T/RtmpVgCm2g/remotesf57f685ebf1/RachelQueen1-Spaniel-bb1bb99/DESCRIPTION’ ...
  
✔  checking for file ‘/private/var/folders/1s/j9ck5c_162s487xcprlxtmdh0000gp/T/RtmpVgCm2g/remotesf57f685ebf1/RachelQueen1-Spaniel-bb1bb99/DESCRIPTION’
## 
  
─  preparing ‘Spaniel’:
## 
  
   checking DESCRIPTION meta-information ...
  
✔  checking DESCRIPTION meta-information
## 
  
─  checking for LF line-endings in source and make files and shell scripts
## 
  
─  checking for empty or unneeded directories
## 
  
   Omitted ‘LazyData’ from DESCRIPTION
## 
  
─  building ‘Spaniel_1.2.0.tar.gz’
## 
  
   
## 
```

```r
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
## [1] 32146  6050
```


# Quality control
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
## DataFrame with 6 rows and 24 columns
##        Sample            Barcode   Section    Spot_Y    Spot_X   Image_Y
##   <character>        <character> <integer> <integer> <integer> <integer>
## 1    Anterior AAACAAGTATCTCCCA-1         1        50       102      7474
## 2    Anterior AAACACCAATAACTGC-1         1        59        19      8552
## 3    Anterior AAACAGAGCGACTCCT-1         1        14        94      3163
## 4    Anterior AAACAGCTTTCAGAAG-1         1        43         9      6636
## 5    Anterior AAACAGGGTCTATATT-1         1        47        13      7115
## 6    Anterior AAACATGGTGAGAGGA-1         1        62         0      8912
##     Image_X   pixel_x   pixel_y       sum  detected     total       sum
##   <integer> <numeric> <numeric> <numeric> <integer> <numeric> <numeric>
## 1      8500   438.898   214.079     13991      4462     13991     13959
## 2      2788   143.959   158.417     39797      8126     39797     39740
## 3      7950   410.499   436.678     29951      6526     29951     29905
## 4      2100   108.434   257.349     42333      8190     42333     42262
## 5      2375   122.633   232.616     35700      8090     35700     35659
## 6      1480    76.420   139.828     22148      6518     22148     22098
##    detected subsets_mt_sum subsets_mt_detected subsets_mt_percent
##   <integer>      <numeric>           <integer>          <numeric>
## 1      4457           1521                  12           10.89620
## 2      8116           3977                  12           10.00755
## 3      6520           4265                  12           14.26183
## 4      8181           2870                  12            6.79097
## 5      8082           1831                  13            5.13475
## 6      6511           2390                  12           10.81546
##   subsets_hb_sum subsets_hb_detected subsets_hb_percent subsets_ribo_sum
##        <numeric>           <integer>          <numeric>        <numeric>
## 1             60                   4           0.429830              826
## 2            831                   6           2.091092             2199
## 3            111                   5           0.371175             1663
## 4            117                   5           0.276844             3129
## 5             73                   5           0.204717             2653
## 6            134                   5           0.606390             1478
##   subsets_ribo_detected subsets_ribo_percent     total
##               <integer>            <numeric> <numeric>
## 1                    85              5.91733     13959
## 2                    89              5.53347     39740
## 3                    88              5.56094     29905
## 4                    88              7.40381     42262
## 5                    90              7.43992     35659
## 6                    84              6.68839     22098
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
## [1] 32146  5803
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
## [1] 32146  5803
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
## [1] 32124  5803
```

# Analysis
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


### Integration

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
##             used   (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells  11463709  612.3   20655705  1103.2   20655705  1103.2
## Vcells 900767995 6872.4 1537221656 11728.1 1322374276 10089.0
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

### Identification of Spatially Variable Features

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
{"columns":[{"label":["p.value"],"name":[1],"type":["dbl"],"align":["right"]},{"label":["FDR"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["gene"],"name":[3],"type":["chr"],"align":["left"]},{"label":["cluster"],"name":[4],"type":["chr"],"align":["left"]}],"data":[{"1":"1.578515e-145","2":"5.070823e-141","3":"Gng4","4":"4"},{"1":"5.907770e-127","2":"9.489061e-123","3":"Gpsm1","4":"4"},{"1":"1.740317e-108","2":"1.863531e-104","3":"Synpr","4":"4"},{"1":"3.401817e-107","2":"2.731999e-103","3":"Pcbp3","4":"4"},{"1":"1.033999e-98","2":"6.643237e-95","3":"Pbx3","4":"4"},{"1":"1.238583e-133","2":"3.978823e-129","3":"Ptgds","4":"6"},{"1":"3.001991e-103","2":"4.821797e-99","3":"Atp1a2","4":"6"},{"1":"3.361581e-76","2":"3.599581e-72","3":"Apoe","4":"6"},{"1":"1.126059e-70","2":"9.043379e-67","3":"Id3","4":"6"},{"1":"6.333538e-70","2":"4.069172e-66","3":"Slc1a2","4":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
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


# Single cell data

We can use a scRNA-seq dataset as a referenced to predict the proportion of different celltypes in the Visium spots. 

Keep in mind that it is important to have a reference that contains all the celltypes you expect to find in your spots. Ideally it should be a scRNAseq reference from the exact same tissue. 

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


# Subset ST for cortex
Since the scRNAseq dataset was generated from the mouse cortex, we will subset the visium dataset in order to select mainly the spots part of the cortex. Note that the integration can also be performed on the whole brain slice, but it would give rise to false positive cell type assignments and and therefore it should be interpreted with more care.


# Integrate with scRNAseq

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


Keep in mind, that the scores are "just" prediction scores, and do not correspond to proportion of cells that are of a certain celltype or similar. It mainly tell you that gene expression in a certain spot is hihgly similar/dissimilar to gene expression of a celltype.

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
##  [1] patchwork_1.1.1             SingleR_1.8.0              
##  [3] Matrix_1.4-0                biomaRt_2.50.0             
##  [5] Spaniel_1.2.0               caret_6.0-90               
##  [7] lattice_0.20-45             SeuratObject_4.0.4         
##  [9] Seurat_4.0.6                scmap_1.16.0               
## [11] scPred_1.9.2                fgsea_1.20.0               
## [13] msigdbr_7.4.1               enrichR_3.0                
## [15] dplyr_1.0.7                 igraph_1.2.11              
## [17] pheatmap_1.0.12             clustree_0.4.4             
## [19] ggraph_2.0.5                reticulate_1.22            
## [21] harmony_1.0                 Rcpp_1.0.8                 
## [23] umap_0.2.7.0                rafalib_1.0.0              
## [25] scDblFinder_1.8.0           DoubletFinder_2.0.3        
## [27] org.Hs.eg.db_3.14.0         AnnotationDbi_1.56.1       
## [29] cowplot_1.1.1               scran_1.22.0               
## [31] scater_1.22.0               ggplot2_3.3.5              
## [33] scuttle_1.4.0               SingleCellExperiment_1.16.0
## [35] SummarizedExperiment_1.24.0 Biobase_2.54.0             
## [37] GenomicRanges_1.46.0        GenomeInfoDb_1.30.0        
## [39] IRanges_2.28.0              S4Vectors_0.32.0           
## [41] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
## [43] matrixStats_0.61.0          RJSONIO_1.3-1.6            
## [45] optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] rsvd_1.0.5                ica_1.0-2                
##   [3] class_7.3-20              ps_1.6.0                 
##   [5] foreach_1.5.1             lmtest_0.9-39            
##   [7] rprojroot_2.0.2           crayon_1.4.2             
##   [9] rhdf5filters_1.6.0        spatstat.core_2.3-2      
##  [11] MASS_7.3-55               nlme_3.1-155             
##  [13] backports_1.4.1           rlang_0.4.12             
##  [15] XVector_0.34.0            ROCR_1.0-11              
##  [17] irlba_2.3.5               callr_3.7.0              
##  [19] limma_3.50.0              filelock_1.0.2           
##  [21] xgboost_1.5.0.1           BiocParallel_1.28.0      
##  [23] rjson_0.2.21              bit64_4.0.5              
##  [25] glue_1.6.0                sctransform_0.3.3        
##  [27] parallel_4.1.2            processx_3.5.2           
##  [29] vipor_0.4.5               spatstat.sparse_2.1-0    
##  [31] spatstat.geom_2.3-1       tidyselect_1.1.1         
##  [33] usethis_2.1.5             fitdistrplus_1.1-6       
##  [35] XML_3.99-0.8              tidyr_1.1.4              
##  [37] zoo_1.8-9                 xtable_1.8-4             
##  [39] magrittr_2.0.1            evaluate_0.14            
##  [41] cli_3.1.0                 zlibbioc_1.40.0          
##  [43] miniUI_0.1.1.1            bslib_0.3.1              
##  [45] rpart_4.1-15              fastmatch_1.1-3          
##  [47] shiny_1.7.1               BiocSingular_1.10.0      
##  [49] xfun_0.29                 askpass_1.1              
##  [51] pkgbuild_1.3.1            cluster_2.1.2            
##  [53] tidygraph_1.2.0           KEGGREST_1.34.0          
##  [55] tibble_3.1.6              ggrepel_0.9.1            
##  [57] listenv_0.8.0             Biostrings_2.62.0        
##  [59] png_0.1-7                 future_1.23.0            
##  [61] ipred_0.9-12              withr_2.4.3              
##  [63] bitops_1.0-7              ggforce_0.3.3            
##  [65] plyr_1.8.6                e1071_1.7-9              
##  [67] dqrng_0.3.0               pROC_1.18.0              
##  [69] pillar_1.6.4              cachem_1.0.6             
##  [71] fs_1.5.2                  kernlab_0.9-29           
##  [73] hdf5r_1.3.5               googleVis_0.6.11         
##  [75] DelayedMatrixStats_1.16.0 vctrs_0.3.8              
##  [77] ellipsis_0.3.2            generics_0.1.1           
##  [79] lava_1.6.10               devtools_2.4.3           
##  [81] tools_4.1.2               beeswarm_0.4.0           
##  [83] munsell_0.5.0             tweenr_1.0.2             
##  [85] proxy_0.4-26              DelayedArray_0.20.0      
##  [87] fastmap_1.1.0             compiler_4.1.2           
##  [89] pkgload_1.2.4             abind_1.4-5              
##  [91] httpuv_1.6.5              sessioninfo_1.2.2        
##  [93] plotly_4.10.0             GenomeInfoDbData_1.2.7   
##  [95] prodlim_2019.11.13        gridExtra_2.3            
##  [97] edgeR_3.36.0              deldir_1.0-6             
##  [99] utf8_1.2.2                later_1.2.0              
## [101] BiocFileCache_2.2.0       recipes_0.1.17           
## [103] jsonlite_1.7.2            scales_1.1.1             
## [105] ScaledMatrix_1.2.0        pbapply_1.5-0            
## [107] sparseMatrixStats_1.6.0   lazyeval_0.2.2           
## [109] promises_1.2.0.1          R.utils_2.11.0           
## [111] goftest_1.2-3             spatstat.utils_2.3-0     
## [113] checkmate_2.0.0           rmarkdown_2.11           
## [115] statmod_1.4.36            Rtsne_0.15               
## [117] uwot_0.1.11               HDF5Array_1.22.0         
## [119] survival_3.2-13           ResidualMatrix_1.4.0     
## [121] yaml_2.2.1                htmltools_0.5.2          
## [123] memoise_2.0.1             locfit_1.5-9.4           
## [125] graphlayouts_0.8.0        here_1.0.1               
## [127] viridisLite_0.4.0         digest_0.6.29            
## [129] assertthat_0.2.1          mime_0.12                
## [131] rappdirs_0.3.3            RSQLite_2.2.8            
## [133] future.apply_1.8.1        remotes_2.4.2            
## [135] data.table_1.14.2         blob_1.2.2               
## [137] R.oo_1.24.0               splines_4.1.2            
## [139] labeling_0.4.2            Rhdf5lib_1.16.0          
## [141] RCurl_1.98-1.5            hms_1.1.1                
## [143] rhdf5_2.38.0              DropletUtils_1.14.0      
## [145] colorspace_2.0-2          BiocManager_1.30.16      
## [147] ggbeeswarm_0.6.0          nnet_7.3-17              
## [149] sass_0.4.0                RANN_2.6.1               
## [151] fansi_1.0.0               parallelly_1.30.0        
## [153] ModelMetrics_1.2.2.2      R6_2.5.1                 
## [155] grid_4.1.2                ggridges_0.5.3           
## [157] lifecycle_1.0.1           formatR_1.11             
## [159] bluster_1.4.0             curl_4.3.2               
## [161] leiden_0.3.9              testthat_3.1.1           
## [163] getopt_1.20.3             jquerylib_0.1.4          
## [165] desc_1.4.0                RcppAnnoy_0.0.19         
## [167] RColorBrewer_1.1-2        iterators_1.0.13         
## [169] stringr_1.4.0             gower_0.2.2              
## [171] htmlwidgets_1.5.4         beachmat_2.10.0          
## [173] polyclip_1.10-0           purrr_0.3.4              
## [175] mgcv_1.8-38               globals_0.14.0           
## [177] openssl_1.4.6             batchelor_1.10.0         
## [179] codetools_0.2-18          lubridate_1.8.0          
## [181] FNN_1.1.3                 metapod_1.2.0            
## [183] randomForest_4.6-14       prettyunits_1.1.1        
## [185] dbplyr_2.1.1              R.methodsS3_1.8.1        
## [187] RSpectra_0.16-0           gtable_0.3.0             
## [189] DBI_1.1.2                 tensor_1.5               
## [191] httr_1.4.2                highr_0.9                
## [193] KernSmooth_2.23-20        stringi_1.7.6            
## [195] progress_1.2.2            reshape2_1.4.4           
## [197] farver_2.1.0              viridis_0.6.2            
## [199] timeDate_3043.102         xml2_1.3.3               
## [201] BiocNeighbors_1.12.0      scattermore_0.7          
## [203] bit_4.0.4                 spatstat.data_2.1-2      
## [205] pkgconfig_2.0.3           babelgene_21.4           
## [207] knitr_1.37
```
