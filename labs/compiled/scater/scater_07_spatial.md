---
title: "Scater/Scran:: Spatial Transcriptomics"
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
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0 36 54.8M   36 19.9M    0     0  34.3M      0  0:00:01 --:--:--  0:00:01 34.2M100 54.8M  100 54.8M    0     0  47.6M      0  0:00:01  0:00:01 --:--:-- 47.7M
## x filtered_feature_bc_matrix/
## x filtered_feature_bc_matrix/barcodes.tsv.gz
## x filtered_feature_bc_matrix/matrix.mtx.gz
## x filtered_feature_bc_matrix/features.tsv.gz
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0100 9480k  100 9480k    0     0  35.3M      0 --:--:-- --:--:-- --:--:-- 35.3M
## x spatial/
## x spatial/tissue_positions_list.csv
## x spatial/tissue_hires_image.png
## x spatial/scalefactors_json.json
## x spatial/aligned_fiducials.jpg
## x spatial/detected_tissue_image.jpg
## x spatial/tissue_lowres_image.png
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0 54.8M    0  2148    0     0  29424      0  0:32:33 --:--:--  0:32:33 29027 81 54.8M   81 44.8M    0     0  42.9M      0  0:00:01  0:00:01 --:--:-- 42.9M100 54.8M  100 54.8M    0     0  43.6M      0  0:00:01  0:00:01 --:--:-- 43.5M
## x filtered_feature_bc_matrix/
## x filtered_feature_bc_matrix/barcodes.tsv.gz
## x filtered_feature_bc_matrix/features.tsv.gz
## x filtered_feature_bc_matrix/matrix.mtx.gz
##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
##                                  Dload  Upload   Total   Spent    Left  Speed
##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0100 9017k  100 9017k    0     0  31.9M      0 --:--:-- --:--:-- --:--:-- 31.9M
## x spatial/
## x spatial/tissue_positions_list.csv
## x spatial/tissue_hires_image.png
## x spatial/scalefactors_json.json
## x spatial/detected_tissue_image.jpg
## x spatial/tissue_lowres_image.png
## x spatial/aligned_fiducials.jpg
```



```r
devtools::install_github("RachelQueen1/Spaniel", ref = "Development")

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
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", 
    host = "https://nov2020.archive.ensembl.org")
annot <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"), 
    mart = mart)

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
## 1      8500   438.898   214.079     13991      4462     13991     13963
## 2      2788   143.959   158.417     39797      8126     39797     39764
## 3      7950   410.499   436.678     29951      6526     29951     29919
## 4      2100   108.434   257.349     42333      8190     42333     42293
## 5      2375   122.633   232.616     35700      8090     35700     35686
## 6      1480    76.420   139.828     22148      6518     22148     22115
##    detected subsets_mt_sum subsets_mt_detected subsets_mt_percent
##   <integer>      <numeric>           <integer>          <numeric>
## 1      4461           1521                  12           10.89307
## 2      8124           3977                  12           10.00151
## 3      6525           4265                  12           14.25516
## 4      8187           2870                  12            6.78599
## 5      8089           1831                  13            5.13086
## 6      6517           2390                  12           10.80714
##   subsets_hb_sum subsets_hb_detected subsets_hb_percent subsets_ribo_sum
##        <numeric>           <integer>          <numeric>        <numeric>
## 1             60                   4           0.429707              826
## 2            831                   6           2.089830             2199
## 3            111                   5           0.371002             1663
## 4            117                   5           0.276642             3129
## 5             73                   5           0.204562             2653
## 6            134                   5           0.605924             1478
##   subsets_ribo_detected subsets_ribo_percent     total
##               <integer>            <numeric> <numeric>
## 1                    85              5.91563     13963
## 2                    89              5.53013     39764
## 3                    88              5.55834     29919
## 4                    88              7.39839     42293
## 5                    90              7.43429     35686
## 6                    84              6.68325     22115
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
## Ncells   8678320  463.5   15881785  848.2  12399795  662.3
## Vcells 189347732 1444.7  332601989 2537.6 332601697 2537.6
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
pred.grun <- SingleR(test = sce.anterior[common_hvgs, ], de.n = 20, ref = allen_reference_sce[common_hvgs, 
    ], labels = allen_reference_sce$subclass, de.method = "wilcox")

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
##  [1] cowplot_1.1.1               patchwork_1.1.1            
##  [3] scater_1.18.0               ggplot2_3.3.3              
##  [5] SingleR_1.4.0               scran_1.18.0               
##  [7] SingleCellExperiment_1.12.0 SummarizedExperiment_1.20.0
##  [9] Biobase_2.50.0              GenomicRanges_1.42.0       
## [11] GenomeInfoDb_1.26.0         IRanges_2.24.0             
## [13] S4Vectors_0.28.0            BiocGenerics_0.36.0        
## [15] MatrixGenerics_1.2.0        matrixStats_0.57.0         
## [17] dplyr_1.0.3                 Matrix_1.3-2               
## [19] biomaRt_2.46.0              Spaniel_1.2.0              
## [21] RJSONIO_1.3-1.4             optparse_1.6.6             
## 
## loaded via a namespace (and not attached):
##   [1] R.utils_2.10.1            reticulate_1.18          
##   [3] tidyselect_1.1.0          RSQLite_2.2.2            
##   [5] AnnotationDbi_1.52.0      htmlwidgets_1.5.3        
##   [7] BiocParallel_1.24.0       grid_4.0.3               
##   [9] Rtsne_0.15                DropletUtils_1.10.2      
##  [11] devtools_2.3.2            munsell_0.5.0            
##  [13] codetools_0.2-18          ica_1.0-2                
##  [15] statmod_1.4.35            future_1.21.0            
##  [17] miniUI_0.1.1.1            batchelor_1.6.0          
##  [19] withr_2.4.0               colorspace_2.0-0         
##  [21] knitr_1.30                Seurat_3.2.3             
##  [23] ROCR_1.0-11               tensor_1.5               
##  [25] listenv_0.8.0             labeling_0.4.2           
##  [27] GenomeInfoDbData_1.2.4    polyclip_1.10-0          
##  [29] farver_2.0.3              bit64_4.0.5              
##  [31] rhdf5_2.34.0              rprojroot_2.0.2          
##  [33] parallelly_1.23.0         vctrs_0.3.6              
##  [35] generics_0.1.0            xfun_0.20                
##  [37] BiocFileCache_1.14.0      R6_2.5.0                 
##  [39] ggbeeswarm_0.6.0          rsvd_1.0.3               
##  [41] locfit_1.5-9.4            rhdf5filters_1.2.0       
##  [43] bitops_1.0-6              spatstat.utils_1.20-2    
##  [45] DelayedArray_0.16.0       assertthat_0.2.1         
##  [47] promises_1.1.1            scales_1.1.1             
##  [49] beeswarm_0.2.3            gtable_0.3.0             
##  [51] beachmat_2.6.0            globals_0.14.0           
##  [53] processx_3.4.5            goftest_1.2-2            
##  [55] rlang_0.4.10              splines_4.0.3            
##  [57] lazyeval_0.2.2            yaml_2.2.1               
##  [59] reshape2_1.4.4            abind_1.4-5              
##  [61] httpuv_1.5.5              tools_4.0.3              
##  [63] usethis_1.6.3             ellipsis_0.3.1           
##  [65] RColorBrewer_1.1-2        sessioninfo_1.1.1        
##  [67] ggridges_0.5.3            Rcpp_1.0.6               
##  [69] plyr_1.8.6                sparseMatrixStats_1.2.0  
##  [71] progress_1.2.2            zlibbioc_1.36.0          
##  [73] purrr_0.3.4               RCurl_1.98-1.2           
##  [75] ps_1.5.0                  prettyunits_1.1.1        
##  [77] rpart_4.1-15              openssl_1.4.3            
##  [79] deldir_0.2-9              viridis_0.5.1            
##  [81] pbapply_1.4-3             zoo_1.8-8                
##  [83] ggrepel_0.9.1             cluster_2.1.0            
##  [85] fs_1.5.0                  magrittr_2.0.1           
##  [87] RSpectra_0.16-0           data.table_1.13.6        
##  [89] scattermore_0.7           ResidualMatrix_1.0.0     
##  [91] lmtest_0.9-38             RANN_2.6.1               
##  [93] fitdistrplus_1.1-3        pkgload_1.1.0            
##  [95] hms_1.0.0                 mime_0.9                 
##  [97] evaluate_0.14             xtable_1.8-4             
##  [99] XML_3.99-0.5              gridExtra_2.3            
## [101] testthat_3.0.1            compiler_4.0.3           
## [103] tibble_3.0.5              KernSmooth_2.23-18       
## [105] crayon_1.3.4              R.oo_1.24.0              
## [107] htmltools_0.5.1           mgcv_1.8-33              
## [109] later_1.1.0.1             tidyr_1.1.2              
## [111] DBI_1.1.1                 formatR_1.7              
## [113] dbplyr_2.0.0              MASS_7.3-53              
## [115] rappdirs_0.3.1            getopt_1.20.3            
## [117] cli_2.2.0                 R.methodsS3_1.8.1        
## [119] igraph_1.2.6              pkgconfig_2.0.3          
## [121] scuttle_1.0.0             plotly_4.9.3             
## [123] xml2_1.3.2                vipor_0.4.5              
## [125] dqrng_0.2.1               XVector_0.30.0           
## [127] stringr_1.4.0             callr_3.5.1              
## [129] digest_0.6.27             sctransform_0.3.2        
## [131] RcppAnnoy_0.0.18          spatstat.data_1.7-0      
## [133] rmarkdown_2.6             leiden_0.3.6             
## [135] edgeR_3.32.0              uwot_0.1.10              
## [137] DelayedMatrixStats_1.12.0 curl_4.3                 
## [139] shiny_1.5.0               lifecycle_0.2.0          
## [141] nlme_3.1-151              jsonlite_1.7.2           
## [143] Rhdf5lib_1.12.0           BiocNeighbors_1.8.0      
## [145] limma_3.46.0              desc_1.2.0               
## [147] viridisLite_0.3.0         askpass_1.1              
## [149] fansi_0.4.2               pillar_1.4.7             
## [151] lattice_0.20-41           fastmap_1.0.1            
## [153] httr_1.4.2                pkgbuild_1.2.0           
## [155] survival_3.2-7            glue_1.4.2               
## [157] remotes_2.2.0             spatstat_1.64-1          
## [159] png_0.1-7                 bluster_1.0.0            
## [161] bit_4.0.4                 HDF5Array_1.18.0         
## [163] stringi_1.5.3             blob_1.2.1               
## [165] BiocSingular_1.6.0        memoise_1.1.0            
## [167] irlba_2.3.3               future.apply_1.7.0
```
