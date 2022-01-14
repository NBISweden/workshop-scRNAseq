---
title: "Scater/Scran: Quality control"
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

***
# Get data

In this tutorial, we will run all tutorials with a set of 6 PBMC 10x datasets from 3 covid-19 patients and 3 healthy controls, the samples have been subsampled to 1500 cells per sample. They are part of the github repo and if you have cloned the repo they should be available in folder: `labs/data/covid_data_GSE149689`. Instructions on how to download them can also be found in the Precourse material.


```r
webpath <- "https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/new_dataset/labs/data/covid_data_GSE149689/sub/"
dir.create("./data/raw", recursive = T)
```

```
## Warning in dir.create("./data/raw", recursive = T): './data/raw' already exists
```

```r
file_list <- c("Normal_PBMC_13.h5", "Normal_PBMC_14.h5", "Normal_PBMC_5.h5", "nCoV_PBMC_15.h5",
    "nCoV_PBMC_17.h5", "nCoV_PBMC_1.h5")
for (i in file_list) {
    download.file(url = paste0(webpath, i), destfile = paste0("./data/raw/", i))
}
```

With data in place, now we can start loading libraries we will use in this tutorial.


```r
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(cowplot))
suppressMessages(require(org.Hs.eg.db))
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = F)
```

```
## Skipping install of 'DoubletFinder' from a github remote, the SHA1 (554097ba) has not changed since last install.
##   Use `force = TRUE` to force installation
```

```r
suppressMessages(require(DoubletFinder))
```

We can first load the data individually by reading directly from HDF5 file format (.h5).


```r
cov.15 <- Seurat::Read10X_h5(filename = "data/raw/nCoV_PBMC_15.h5", use.names = T)
```

```
## Warning in sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x =
## counts[]), : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

```r
cov.1 <- Seurat::Read10X_h5(filename = "data/raw/nCoV_PBMC_1.h5", use.names = T)
```

```
## Warning in sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x =
## counts[]), : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

```r
cov.17 <- Seurat::Read10X_h5(filename = "data/raw/nCoV_PBMC_17.h5", use.names = T)
```

```
## Warning in sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x =
## counts[]), : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

```r
ctrl.5 <- Seurat::Read10X_h5(filename = "data/raw/Normal_PBMC_5.h5", use.names = T)
```

```
## Warning in sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x =
## counts[]), : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

```r
ctrl.13 <- Seurat::Read10X_h5(filename = "data/raw/Normal_PBMC_13.h5", use.names = T)
```

```
## Warning in sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x =
## counts[]), : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

```r
ctrl.14 <- Seurat::Read10X_h5(filename = "data/raw/Normal_PBMC_14.h5", use.names = T)
```

```
## Warning in sparseMatrix(i = indices[] + 1, p = indptr[], x = as.numeric(x =
## counts[]), : 'giveCsparse' has been deprecated; setting 'repr = "T"' for you
```

***
# Create one merged object

We can now load the expression matricies into objects and then merge them into a single merged object. Each analysis workflow (Seurat, Scater, Scranpy, etc) has its own way of storing data. We will add dataset labels as cell.ids just in case you have overlapping barcodes between the datasets. After that we add a column `Chemistry` in the metadata for plotting later on.


```r
sce <- SingleCellExperiment(assays = list(counts = cbind(cov.1, cov.15, cov.17, ctrl.5,
    ctrl.13, ctrl.14)))
dim(sce)
```

```
## [1] 33538  9000
```

```r
# Adding metadata
sce@colData$sample <- unlist(sapply(c("cov.1", "cov.15", "cov.17", "ctrl.5", "ctrl.13",
    "ctrl.14"), function(x) rep(x, ncol(get(x)))))
sce@colData$type <- ifelse(grepl("cov", sce@colData$sample), "Covid", "Control")
```

Once you have created the merged Seurat object, the count matrices and individual count matrices and objects are not needed anymore. It is a good idea to remove them and run garbage collect to free up some memory.


```r
# remove all objects that will not be used.
rm(cov.15, cov.1, cov.17, ctrl.5, ctrl.13, ctrl.14)

# run garbage collect to free up memory
gc()
```

```
##            used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells  8515416 454.8   15528642 829.4 12165947 649.8
## Vcells 34181126 260.8   96012892 732.6 71370890 544.6
```


 Here it is how the count matrix and the metatada look like for every cell.


```r
head(counts(sce)[, 1:10])

head(sce@colData, 10)
```

```
## 6 x 10 sparse Matrix of class "dgCMatrix"
##                                
## MIR1302-2HG . . . . . . . . . .
## FAM138A     . . . . . . . . . .
## OR4F5       . . . . . . . . . .
## AL627309.1  . . . . . . . . . .
## AL627309.3  . . . . . . . . . .
## AL627309.2  . . . . . . . . . .
## DataFrame with 10 rows and 2 columns
##                         sample        type
##                    <character> <character>
## AGGGTCCCATGACCCG-1       cov.1       Covid
## TACCCACAGCGGGTTA-1       cov.1       Covid
## CCCAACTTCATATGGC-1       cov.1       Covid
## TCAAGTGTCCGAACGC-1       cov.1       Covid
## ATTCCTAGTGACTGTT-1       cov.1       Covid
## GTGTTCCGTGGGCTCT-1       cov.1       Covid
## CCTAAGACAGATTAAG-1       cov.1       Covid
## AATAGAGAGGGTTAGC-1       cov.1       Covid
## GGGTCACTCACCTACC-1       cov.1       Covid
## TCCTCTTGTACAGTCT-1       cov.1       Covid
```


***
# Calculate QC

Having the data in a suitable format, we can start calculating some quality metrics. We can for example calculate the percentage of mitocondrial and ribosomal genes per cell and add to the metadata. This will be helpfull to visualize them across different metadata parameteres (i.e. datasetID and chemistry version). There are several ways of doing this, and here manually calculate the proportion of mitochondrial reads and add to the metadata table.

Citing from "Simple Single Cell" workflows (Lun, McCarthy & Marioni, 2017): "High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane."

 First, let Scran calculate some general qc-stats for genes and cells with the function `perCellQCMetrics`. It can also calculate proportion of counts for specific gene subsets, so first we need to define which genes are mitochondrial, ribosomal and hemoglogin.


```r
# Mitochondrial genes
mito_genes <- rownames(sce)[grep("^MT-", rownames(sce))]

# Ribosomal genes
ribo_genes <- rownames(sce)[grep("^RP[SL]", rownames(sce))]

# Hemoglobin genes - includes all genes starting with HB except HBP.
hb_genes <- rownames(sce)[grep("^HB[^(P)]", rownames(sce))]
```


```r
sce <- addPerCellQC(sce, flatten = T, subsets = list(mt = mito_genes, hb = hb_genes,
    ribo = ribo_genes))

head(colData(sce))
```

```
## DataFrame with 6 rows and 14 columns
##                         sample        type       sum  detected subsets_mt_sum
##                    <character> <character> <numeric> <integer>      <numeric>
## AGGGTCCCATGACCCG-1       cov.1       Covid      7698      2140            525
## TACCCACAGCGGGTTA-1       cov.1       Covid     13416      3391            952
## CCCAACTTCATATGGC-1       cov.1       Covid     16498      3654           1253
## TCAAGTGTCCGAACGC-1       cov.1       Covid      1425       608            141
## ATTCCTAGTGACTGTT-1       cov.1       Covid      7535      1808            470
## GTGTTCCGTGGGCTCT-1       cov.1       Covid      4378      1345            352
##                    subsets_mt_detected subsets_mt_percent subsets_hb_sum
##                              <integer>          <numeric>      <numeric>
## AGGGTCCCATGACCCG-1                  11            6.81995              2
## TACCCACAGCGGGTTA-1                  11            7.09600              6
## CCCAACTTCATATGGC-1                  12            7.59486              1
## TCAAGTGTCCGAACGC-1                  10            9.89474              1
## ATTCCTAGTGACTGTT-1                  11            6.23756              4
## GTGTTCCGTGGGCTCT-1                  10            8.04020              1
##                    subsets_hb_detected subsets_hb_percent subsets_ribo_sum
##                              <integer>          <numeric>        <numeric>
## AGGGTCCCATGACCCG-1                   1         0.02598077             2564
## TACCCACAGCGGGTTA-1                   2         0.04472272             2264
## CCCAACTTCATATGGC-1                   1         0.00606134             2723
## TCAAGTGTCCGAACGC-1                   1         0.07017544              444
## ATTCCTAGTGACTGTT-1                   3         0.05308560             3397
## GTGTTCCGTGGGCTCT-1                   1         0.02284148             1588
##                    subsets_ribo_detected subsets_ribo_percent     total
##                                <integer>            <numeric> <numeric>
## AGGGTCCCATGACCCG-1                    82              33.3074      7698
## TACCCACAGCGGGTTA-1                    85              16.8754     13416
## CCCAACTTCATATGGC-1                    87              16.5050     16498
## TCAAGTGTCCGAACGC-1                    68              31.1579      1425
## ATTCCTAGTGACTGTT-1                    81              45.0829      7535
## GTGTTCCGTGGGCTCT-1                    79              36.2723      4378
```

 Here is an example on how to calculate proportion mitochondria in another way:


```r
# Way2: Doing it manually
sce@colData$percent_mito <- Matrix::colSums(counts(sce)[mito_genes, ])/sce@colData$total
```


***
# Plot QC

Now we can plot some of the QC-features as violin plots.


```r
# total is total UMIs per cell detected is number of detected genes.  the
# different gene subset percentages are listed as subsets_mt_percent etc.

plot_grid(plotColData(sce, y = "detected", x = "sample", colour_by = "sample"), plotColData(sce,
    y = "total", x = "sample", colour_by = "sample"), plotColData(sce, y = "subsets_mt_percent",
    x = "sample", colour_by = "sample"), plotColData(sce, y = "subsets_ribo_percent",
    x = "sample", colour_by = "sample"), plotColData(sce, y = "subsets_hb_percent",
    x = "sample", colour_by = "sample"), ncol = 3)
```

![](scater_01_qc_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

As you can see, there is quite some difference in quality for the 4 datasets, with for instance the covid_15 sample having fewer cells with many detected genes and more mitochondrial content. As the ribosomal proteins are highly expressed they will make up a larger proportion of the transcriptional landscape when fewer of the lowly expressed genes are detected. And we can plot the different QC-measures as scatter plots.


```r
plotColData(sce, x = "total", y = "detected", colour_by = "sample")
```

![](scater_01_qc_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

<style>
div.blue { background-color:#e6f0ff; border-radius: 5px; padding: 10px;}
</style>
<div class = "blue">
**Your turn**

Plot additional QC stats that we have calculated as scatter plots. How are the different measures correlated? Can you explain why?
</div>

***
# Filtering

## Detection-based filtering

In scran, we can use the function `quickPerCellQC` to filter out outliers from distributions of qc stats, such as dected genes, gene subsets etc. But in this case, we will take one setting at a time and run through the steps of filtering cells.

A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. Here we will only consider cells with at least 200 detected genes and genes need to be expressed in at least 3 cells. Please note that those values are highly dependent on the library preparation method used.



```r
dim(sce)

selected_c <- colnames(sce)[sce$detected > 200]
selected_f <- rownames(sce)[Matrix::rowSums(counts(sce)) > 3]

sce.filt <- sce[selected_f, selected_c]
dim(sce.filt)
```

```
## [1] 33538  9000
## [1] 18147  7973
```

 Extremely high number of detected genes could indicate doublets. However, depending on the cell type composition in your sample, you may have cells with higher number of genes (and also higher counts) from one cell type. <br>In this case, we will run doublet prediction further down, so we will skip this step now, but the code below is an example of how it can be run:


```r
# skip for now and run doublet detection instead...

# high.det.v3 <- sce.filt$nFeatures > 4100 high.det.v2 <- (sce.filt$nFeatures >
# 2000) & (sce.filt$sample_id == 'v2.1k')

# remove these cells sce.filt <- sce.filt[ , (!high.det.v3) & (!high.det.v2)]

# check number of cells ncol(sce.filt)
```

Additionally, we can also see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.

In scater, you can also use the function `plotHighestExprs()` to plot the gene contribution, but the function is quite slow.


```r
# Compute the relative expression of each gene per cell

# Use sparse matrix operations, if your dataset is large, doing matrix
# devisions the regular way will take a very long time.
C = counts(sce)
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
```

![](scater_01_qc_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
rm(C)

# also, there is the option of running the function 'plotHighestExprs' in the
# scater package, however, this function takes very long to execute.
```

As you can see, MALAT1 constitutes up to 30% of the UMIs from a single cell and the other top genes are mitochondrial and ribosomal genes. It is quite common that nuclear lincRNAs have correlation with quality and mitochondrial reads, so high detection of MALAT1 may be a technical issue. Let us assemble some information about such genes, which are important for quality control and downstream filtering.

## Mito/Ribo filtering

We also have quite a lot of cells with high proportion of mitochondrial and low proportion ofribosomal reads. It could be wise to remove those cells, if we have enough cells left after filtering. <br>Another option would be to either remove all mitochondrial reads from the dataset and hope that the remaining genes still have enough biological signal. <br>A third option would be to just regress out the `percent_mito` variable during scaling. In this case we had as much as 99.7% mitochondrial reads in some of the cells, so it is quite unlikely that there is much cell type signature left in those. <br>Looking at the plots, make reasonable decisions on where to draw the cutoff. In this case, the bulk of the cells are below 20% mitochondrial reads and that will be used as a cutoff. We will also remove cells with less than 5% ribosomal reads.


```r
selected_mito <- sce.filt$subsets_mt_percent < 30
selected_ribo <- sce.filt$subsets_ribo_percent > 5

# and subset the object to only keep those cells
sce.filt <- sce.filt[, selected_mito & selected_ribo]
dim(sce.filt)
```

```
## [1] 18147  5896
```

 As you can see, a large proportion of sample covid_15 is filtered out. Also, there is still quite a lot of variation in `percent_mito`, so it will have to be dealt with in the data analysis step. We can also notice that the `percent_ribo` are also highly variable, but that is expected since different cell types have different proportions of ribosomal content, according to their function.

## Plot filtered QC

Lets plot the same QC-stats another time.


```r
plot_grid(plotColData(sce, y = "detected", x = "sample", colour_by = "sample"), plotColData(sce,
    y = "total", x = "sample", colour_by = "sample"), plotColData(sce, y = "subsets_mt_percent",
    x = "sample", colour_by = "sample"), plotColData(sce, y = "subsets_ribo_percent",
    x = "sample", colour_by = "sample"), plotColData(sce, y = "subsets_hb_percent",
    x = "sample", colour_by = "sample"), ncol = 3)
```

![](scater_01_qc_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

## Filter genes

As the level of expression of mitochondrial and MALAT1 genes are judged as mainly technical, it can be wise to remove them from the dataset bofore any further analysis.


```r
dim(sce.filt)

# Filter MALAT1
sce.filt <- sce.filt[!grepl("MALAT1", rownames(sce.filt)), ]

# Filter Mitocondrial
sce.filt <- sce.filt[!grepl("^MT-", rownames(sce.filt)), ]

# Filter Ribossomal gene (optional if that is a problem on your data) sce.filt
# <- sce.filt[ ! grepl('^RP[SL]', rownames(sce.filt)), ]

# Filter Hemoglobin gene
sce.filt <- sce.filt[!grepl("^HB[^(P)]", rownames(sce.filt)), ]

dim(sce.filt)
```

```
## [1] 18147  5896
## [1] 18121  5896
```


# Sample sex

When working with human or animal samples, you should ideally constrain you experiments to a single sex to avoid including sex bias in the conclusions. However this may not always be possible. By looking at reads from chromosomeY (males) and XIST (X-inactive specific transcript) expression (mainly female) it is quite easy to determine per sample which sex it is. It can also bee a good way to detect if there has been any sample mixups, if the sample metadata sex does not agree with the computational predictions.

To get choromosome information for all genes, you should ideally parse the information from the gtf file that you used in the mapping pipeline as it has the exact same annotation version/gene naming. However, it may not always be available, as in this case where we have downloaded public data. Hence, we will use biomart to fetch chromosome information. 
As the biomart instances quite often are unresponsive, you can try the code below, but if it fails, we have the file with gene annotations on github [here](https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/labs/misc/genes.table.csv). Make sure you put it at the correct location for the path `genes.file` to work. 


```r
genes.file = "data/results/genes.table.csv"

if (!file.exists(genes.file)) {
    suppressMessages(require(biomaRt))

    # initialize connection to mart, may take some time if the sites are
    # unresponsive.
    mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

    # fetch chromosome info plus some other annotations
    genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
        "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
        useCache = F))

    if (!dir.exists("data/results")) {
        dir.create("data/results")
    }
    if (is.data.frame(genes.table)) {
        write.csv(genes.table, file = genes.file)
    }

    if (!file.exists(genes.file)) {
        download.file("https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/misc/genes.table.csv",
            destfile = "data/results/genes.table.csv")
        genes.table = read.csv(genes.file)
    }
} else {
    genes.table = read.csv(genes.file)
}

genes.table <- genes.table[genes.table$external_gene_name %in% rownames(sce.filt),
    ]
```

Now that we have the chromosome information, we can calculate per cell the proportion of reads that comes from chromosome Y.


```r
chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]

sce.filt@colData$pct_chrY = Matrix::colSums(counts(sce.filt)[chrY.gene, ])/colSums(counts(sce.filt))
```

Then plot XIST expression vs chrY proportion. As you can see, the samples are clearly on either side, even if some cells do not have detection of either.


```r
# as plotColData cannot take an expression vs metadata, we need to add in XIST
# expression to colData
sce.filt@colData$XIST = counts(sce.filt)["XIST", ]/colSums(counts(sce.filt)) * 10000

plotColData(sce.filt, "XIST", "pct_chrY")
```

![](scater_01_qc_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

Plot as violins.


```r
plot_grid(plotColData(sce.filt, y = "XIST", x = "sample", colour_by = "sample"),
    plotColData(sce.filt, y = "pct_chrY", x = "sample", colour_by = "sample"), ncol = 2)
```

![](scater_01_qc_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

Here, we can see clearly that we have two males and 4 females, can you see which samples they are? 
Do you think this will cause any problems for downstream analysis? Discuss with your group: what would be the best way to deal with this type of sex bias?

# Calculate cell-cycle scores

We here perform cell cycle scoring. To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list. Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.


```r
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package = "scran"))
anno <- select(org.Hs.eg.db, keys = rownames(sce.filt), keytype = "SYMBOL", column = "ENSEMBL")
```

```
## 'select()' returned 1:many mapping between keys and columns
```

```r
ensembl <- anno$ENSEMBL[match(rownames(sce.filt), anno$SYMBOL)]

# Use only genes related to biological process cell cycle to speed up
# https://www.ebi.ac.uk/QuickGO/term/GO:0007049 = cell cycle (BP,Biological
# Process)
GOs <- na.omit(select(org.Hs.eg.db, keys = na.omit(ensembl), keytype = "ENSEMBL",
    column = "GO"))
```

```
## 'select()' returned many:many mapping between keys and columns
```

```r
GOs <- GOs[GOs$GO == "GO:0007049", "ENSEMBL"]
hs.pairs <- lapply(hs.pairs, function(x) {
    x[rowSums(apply(x, 2, function(i) i %in% GOs)) >= 1, ]
})
str(hs.pairs)
cc.ensembl <- ensembl[ensembl %in% GOs]  #This is the fastest (less genes), but less accurate too
# cc.ensembl <- ensembl[ ensembl %in% unique(unlist(hs.pairs))]


assignments <- cyclone(sce.filt[ensembl %in% cc.ensembl, ], hs.pairs, gene.names = ensembl[ensembl %in%
    cc.ensembl])
sce.filt$G1.score <- assignments$scores$G1
sce.filt$G2M.score <- assignments$scores$G2M
sce.filt$S.score <- assignments$scores$S
```

```
## List of 3
##  $ G1 :'data.frame':	6916 obs. of  2 variables:
##   ..$ first : chr [1:6916] "ENSG00000100519" "ENSG00000100519" "ENSG00000100519" "ENSG00000100519" ...
##   ..$ second: chr [1:6916] "ENSG00000065135" "ENSG00000080345" "ENSG00000101266" "ENSG00000124486" ...
##  $ S  :'data.frame':	8889 obs. of  2 variables:
##   ..$ first : chr [1:8889] "ENSG00000255302" "ENSG00000119969" "ENSG00000179051" "ENSG00000127586" ...
##   ..$ second: chr [1:8889] "ENSG00000100519" "ENSG00000100519" "ENSG00000100519" "ENSG00000136856" ...
##  $ G2M:'data.frame':	6894 obs. of  2 variables:
##   ..$ first : chr [1:6894] "ENSG00000100519" "ENSG00000136856" "ENSG00000136856" "ENSG00000136856" ...
##   ..$ second: chr [1:6894] "ENSG00000146457" "ENSG00000007968" "ENSG00000227268" "ENSG00000101265" ...
```

We can now plot a violin plot for the cell cycle scores as well.


```r
plot_grid(plotColData(sce.filt, y = "G2M.score", x = "G1.score", colour_by = "sample"),
    plotColData(sce.filt, y = "G2M.score", x = "sample", colour_by = "sample"), plotColData(sce.filt,
        y = "G1.score", x = "sample", colour_by = "sample"), plotColData(sce.filt,
        y = "S.score", x = "sample", colour_by = "sample"), ncol = 4)
```

```
## Warning: Removed 491 rows containing missing values (geom_point).
```

```
## Warning: Removed 418 rows containing non-finite values (stat_ydensity).
```

```
## Warning: Removed 418 rows containing missing values (position_quasirandom).
```

```
## Warning: Removed 456 rows containing non-finite values (stat_ydensity).
```

```
## Warning: Removed 456 rows containing missing values (position_quasirandom).
```

```
## Warning: Removed 336 rows containing non-finite values (stat_ydensity).
```

```
## Warning: Removed 336 rows containing missing values (position_quasirandom).
```

![](scater_01_qc_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

Cyclone predicts most cells as G1, but also quite a lot of cells with high S-Phase scores. Compare to results with Seurat and Scanpy and see how different predictors will give clearly different results.

# Predict doublets

Doublets/Mulitples of cells in the same well/droplet is a common issue in scRNAseq protocols. Especially in droplet-based methods whith overloading of cells. In a typical 10x experiment the proportion of doublets is linearly dependent on the amount of loaded cells. As  indicated from the Chromium user guide, doublet rates are about as follows:
![](../../figs/10x_doublet_rate.png)
Most doublet detectors simulates doublets by merging cell counts and predicts doublets as cells that have similar embeddings as the simulated doublets. Most such packages need an assumption about the number/proportion of expected doublets in the dataset. The data you are using is subsampled, but the orignial datasets contained about 5 000 cells per sample, hence we can assume that they loaded about 9 000 cells and should have a doublet rate at about 4%.
OBS! Ideally doublet prediction should be run on each sample separately, especially if your different samples have different proportions of celltypes. In this case, the data is subsampled so we have very few cells per sample and all samples are sorted PBMCs so it is okay to run them together.

There is a method to predict if a cluster consists of mainly doublets `findDoubletClusters()`, but we can also predict individual cells based on simulations using the function `computeDoubletDensity()` which we will do here.
Doublet detection will be performed using PCA, so we need to first normalize the data and run variable gene detection, as well as UMAP for visualization. These steps will be explored in more detail in coming exercises.
There is a method to predict if a cluster consists of mainly doublets `findDoubletClusters()`, but we can also predict individual cells based on simulations using the function `computeDoubletDensity()` which we will do here.
Doublet detection will be performed using PCA, so we need to first normalize the data and run variable gene detection, as well as UMAP for visualization. These steps will be explored in more detail in coming exercises.


```r
sce.filt <- logNormCounts(sce.filt)
dec <- modelGeneVar(sce.filt, block = sce.filt$sample)
hvgs = getTopHVGs(dec, n = 2000)

sce.filt <- runPCA(sce.filt, subset_row = hvgs)

sce.filt <- runUMAP(sce.filt, pca = 10)
```


```r
BiocManager::install("scDblFinder", update = F)
suppressPackageStartupMessages(require(scDblFinder))

# run computeDoubletDensity with 10 principal components.
sce.filt <- scDblFinder(sce.filt, dims = 10)
```


```r
plot_grid(plotUMAP(sce.filt, colour_by = "scDblFinder.score"), plotUMAP(sce.filt,
    colour_by = "scDblFinder.class"), plotUMAP(sce.filt, colour_by = "sample"), ncol = 3)
```

![](scater_01_qc_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

Now, lets remove all predicted doublets from our data.


```r
sce.filt = sce.filt[, sce.filt$scDblFinder.score < 2]

dim(sce.filt)
```

```
## [1] 18121  5896
```

# Save data
Finally, lets save the QC-filtered data for further analysis. Create output directory `results` and save data to that folder.


```r
dir.create("data/results", showWarnings = F)

saveRDS(sce.filt, "data/results/covid_qc.rds")
```

We should expect that two cells have more detected genes than a single cell, lets check if our predicted doublets also have more detected genes in general.

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
##  [1] scDblFinder_1.8.0           DoubletFinder_2.0.3        
##  [3] org.Hs.eg.db_3.14.0         AnnotationDbi_1.56.1       
##  [5] cowplot_1.1.1               scran_1.22.0               
##  [7] scater_1.22.0               ggplot2_3.3.5              
##  [9] scuttle_1.4.0               SingleCellExperiment_1.16.0
## [11] SummarizedExperiment_1.24.0 Biobase_2.54.0             
## [13] GenomicRanges_1.46.0        GenomeInfoDb_1.30.0        
## [15] IRanges_2.28.0              S4Vectors_0.32.0           
## [17] BiocGenerics_0.40.0         MatrixGenerics_1.6.0       
## [19] matrixStats_0.61.0          RJSONIO_1.3-1.6            
## [21] optparse_1.7.1             
## 
## loaded via a namespace (and not attached):
##   [1] utf8_1.2.2                reticulate_1.22          
##   [3] tidyselect_1.1.1          RSQLite_2.2.8            
##   [5] htmlwidgets_1.5.4         grid_4.1.2               
##   [7] BiocParallel_1.28.0       Rtsne_0.15               
##   [9] munsell_0.5.0             ScaledMatrix_1.2.0       
##  [11] codetools_0.2-18          ica_1.0-2                
##  [13] xgboost_1.5.0.1           statmod_1.4.36           
##  [15] future_1.23.0             miniUI_0.1.1.1           
##  [17] withr_2.4.3               colorspace_2.0-2         
##  [19] highr_0.9                 knitr_1.37               
##  [21] Seurat_4.0.6              ROCR_1.0-11              
##  [23] tensor_1.5                listenv_0.8.0            
##  [25] labeling_0.4.2            GenomeInfoDbData_1.2.7   
##  [27] polyclip_1.10-0           bit64_4.0.5              
##  [29] farver_2.1.0              parallelly_1.30.0        
##  [31] vctrs_0.3.8               generics_0.1.1           
##  [33] xfun_0.29                 R6_2.5.1                 
##  [35] ggbeeswarm_0.6.0          rsvd_1.0.5               
##  [37] locfit_1.5-9.4            hdf5r_1.3.5              
##  [39] bitops_1.0-7              spatstat.utils_2.3-0     
##  [41] cachem_1.0.6              DelayedArray_0.20.0      
##  [43] assertthat_0.2.1          promises_1.2.0.1         
##  [45] scales_1.1.1              beeswarm_0.4.0           
##  [47] gtable_0.3.0              beachmat_2.10.0          
##  [49] globals_0.14.0            goftest_1.2-3            
##  [51] rlang_0.4.12              splines_4.1.2            
##  [53] lazyeval_0.2.2            spatstat.geom_2.3-1      
##  [55] BiocManager_1.30.16       yaml_2.2.1               
##  [57] reshape2_1.4.4            abind_1.4-5              
##  [59] httpuv_1.6.5              tools_4.1.2              
##  [61] ellipsis_0.3.2            spatstat.core_2.3-2      
##  [63] jquerylib_0.1.4           RColorBrewer_1.1-2       
##  [65] ggridges_0.5.3            Rcpp_1.0.8               
##  [67] plyr_1.8.6                sparseMatrixStats_1.6.0  
##  [69] zlibbioc_1.40.0           purrr_0.3.4              
##  [71] RCurl_1.98-1.5            rpart_4.1-15             
##  [73] deldir_1.0-6              pbapply_1.5-0            
##  [75] viridis_0.6.2             zoo_1.8-9                
##  [77] SeuratObject_4.0.4        ggrepel_0.9.1            
##  [79] cluster_2.1.2             magrittr_2.0.1           
##  [81] RSpectra_0.16-0           data.table_1.14.2        
##  [83] scattermore_0.7           lmtest_0.9-39            
##  [85] RANN_2.6.1                fitdistrplus_1.1-6       
##  [87] patchwork_1.1.1           mime_0.12                
##  [89] evaluate_0.14             xtable_1.8-4             
##  [91] gridExtra_2.3             compiler_4.1.2           
##  [93] tibble_3.1.6              KernSmooth_2.23-20       
##  [95] crayon_1.4.2              htmltools_0.5.2          
##  [97] mgcv_1.8-38               later_1.2.0              
##  [99] tidyr_1.1.4               DBI_1.1.2                
## [101] formatR_1.11              MASS_7.3-55              
## [103] Matrix_1.4-0              getopt_1.20.3            
## [105] parallel_4.1.2            metapod_1.2.0            
## [107] igraph_1.2.11             pkgconfig_2.0.3          
## [109] plotly_4.10.0             spatstat.sparse_2.1-0    
## [111] vipor_0.4.5               bslib_0.3.1              
## [113] dqrng_0.3.0               XVector_0.34.0           
## [115] stringr_1.4.0             digest_0.6.29            
## [117] sctransform_0.3.3         RcppAnnoy_0.0.19         
## [119] spatstat.data_2.1-2       Biostrings_2.62.0        
## [121] rmarkdown_2.11            leiden_0.3.9             
## [123] uwot_0.1.11               edgeR_3.36.0             
## [125] DelayedMatrixStats_1.16.0 curl_4.3.2               
## [127] shiny_1.7.1               lifecycle_1.0.1          
## [129] nlme_3.1-155              jsonlite_1.7.2           
## [131] BiocNeighbors_1.12.0      viridisLite_0.4.0        
## [133] limma_3.50.0              fansi_1.0.0              
## [135] pillar_1.6.4              lattice_0.20-45          
## [137] KEGGREST_1.34.0           fastmap_1.1.0            
## [139] httr_1.4.2                survival_3.2-13          
## [141] glue_1.6.0                remotes_2.4.2            
## [143] png_0.1-7                 bluster_1.4.0            
## [145] bit_4.0.4                 stringi_1.7.6            
## [147] sass_0.4.0                blob_1.2.2               
## [149] BiocSingular_1.10.0       memoise_2.0.1            
## [151] dplyr_1.0.7               irlba_2.3.5              
## [153] future.apply_1.8.1
```
