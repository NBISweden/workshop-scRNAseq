Quality Control
================

Created by: Åsa Björklund

# Overview

#### Quality control of data for filtering cells using Seurat and Scater packages.

In this tutorial we will look at different ways of doing filtering and
cell and exploring variablility in the data.

The first part is using Seurat (<https://satijalab.org/seurat/>) for
visualizing QC-measures and filtering cells. However, we will not go
into depth in how to use the Seurat package as this will be covered in
other tutorials.

The second part will explore the scater package
(<https://bioconductor.org/packages/release/bioc/html/scater.html>) in
some more detail. Looking at different ways of visualizing QC-stats and
exploring variation in the data.

### Dataset

For this tutorial we will use 3 different PBMC datasets from the 10x
Genomics website
(<https://support.10xgenomics.com/single-cell-gene-expression/datasets>).

  - 1k PBMCs using 10x v2 chemistry
  - 1k PBMCs using 10x v3 chemistry
  - 1k PBMCs using 10x v3 chemistry in combination with cell surface
    proteins, but disregarding the protein data and only looking at gene
    expression.

To use the exact code, put the files in folder `data/3pbmc/`.

    cd data
    curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5
    curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
    curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5

Load required packages:

``` r
suppressMessages(require(Seurat))
suppressMessages(require(scater))
suppressMessages(require(Matrix))
suppressMessages(require(reshape2))
```

#### Read data

Here, we use the function Read10X\_h5 to read in the expression
matrices.

``` r
v3.1k <- Read10X_h5("../data/3pbmc/pbmc_1k_v3_filtered_feature_bc_matrix.h5", use.names = T)
v2.1k <- Read10X_h5("../data/3pbmc/pbmc_1k_v2_filtered_feature_bc_matrix.h5", use.names = T)
p3.1k <- Read10X_h5("../data/3pbmc/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5", use.names = T)
```

    ## Genome matrix has multiple modalities, returning a list of matrices for this genome

``` r
# select only gene expression data from the CITE-seq data.
p3.1k <- p3.1k$`Gene Expression`
```

# Seurat

### Create Seurat object

First, create Seurat objects for each of the datasets, and then merge
into one large seurat object.

``` r
sdata.v2.1k <- CreateSeuratObject(v2.1k, project = "v2.1k")
sdata.v3.1k <- CreateSeuratObject(v3.1k, project = "v3.1k")
sdata.p3.1k <- CreateSeuratObject(p3.1k, project = "p3.1k")

# merge into one single seurat object. Add cell ids just in case you have overlapping barcodes between the datasets.
alldata <- merge(sdata.v2.1k, c(sdata.v3.1k,sdata.p3.1k), add.cell.ids=c("v2.1k","v3.1k","p3.1k"))

# also add in a metadata column that indicates v2 vs v3 chemistry
chemistry <- rep("v3",ncol(alldata))
chemistry[Idents(alldata) == "v2.1k"] <- "v2"
alldata <- AddMetaData(alldata, chemistry, col.name = "Chemistry")
alldata
```

    ## An object of class Seurat 
    ## 33538 features across 2931 samples within 1 assay 
    ## Active assay: RNA (33538 features)

``` r
# check number of cells from each sample, is stored in the orig.ident slot of metadata and is autmatically set as active ident.
table(Idents(alldata))
```

    ## 
    ## p3.1k v2.1k v3.1k 
    ##   713   996  1222

#### Calculate mitochondrial proportion

Seurat automatically calculates some QC-stats, like number of UMIs and
features per cell. Stored in columns nCount\_RNA & nFeature\_RNA of the
metadata.

``` r
head(alldata@meta.data)
```

    ##                          orig.ident nCount_RNA nFeature_RNA Chemistry
    ## v2.1k_AAACCTGAGCGCTCCA-1      v2.1k       6631         2029        v2
    ## v2.1k_AAACCTGGTGATAAAC-1      v2.1k       2196          881        v2
    ## v2.1k_AAACGGGGTTTGTGTG-1      v2.1k       2700          791        v2
    ## v2.1k_AAAGATGAGTACTTGC-1      v2.1k       3551         1183        v2
    ## v2.1k_AAAGCAAGTCTCTTAT-1      v2.1k       3080         1333        v2
    ## v2.1k_AAAGCAATCCACGAAT-1      v2.1k       5769         1556        v2

We will manually calculate the proportion of mitochondrial reads and add
to the metadata table.

``` r
mt.genes <- rownames(alldata)[grep("^MT-",rownames(alldata))]
C<-GetAssayData(object = alldata, slot = "counts")

percent.mito <- Matrix::colSums(C[mt.genes,])/Matrix::colSums(C)*100
alldata <- AddMetaData(alldata, percent.mito, col.name = "percent.mito")
```

#### Calculate ribosomal proportion

In the same manner we will calculate the proportion gene expression that
comes from ribosomal proteins. NOTE - add text on why\!

``` r
rb.genes <- rownames(alldata)[grep("^RP[SL]",rownames(alldata))]
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
alldata <- AddMetaData(alldata, percent.ribo, col.name = "percent.ribo")
```

### Plot QC

Now we can plot some of the QC-features as violin plots

``` r
VlnPlot(alldata, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot-1.png)<!-- -->

``` r
VlnPlot(alldata, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot-2.png)<!-- -->

``` r
VlnPlot(alldata, features = "percent.mito", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot-3.png)<!-- -->

``` r
VlnPlot(alldata, features = "percent.ribo", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot-4.png)<!-- -->

As you can see, the v2 chemistry gives lower gene detection, but higher
detection of ribosomal proteins. As the ribosomal proteins are highly
expressed they will make up a larger proportion of the transcriptional
landscape when fewer of the lowly expressed genes are detected.

And we can plot the different QC-measures as scatter plots

``` r
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
```

![](Quality_control_files/figure-gfm/scatter-1.png)<!-- -->

``` r
FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mito")
```

![](Quality_control_files/figure-gfm/scatter-2.png)<!-- -->

``` r
FeatureScatter(alldata, feature1="percent.ribo", feature2="nFeature_RNA")
```

![](Quality_control_files/figure-gfm/scatter-3.png)<!-- -->

We can also subset the data to only plot one sample.

``` r
FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cells = WhichCells(alldata, expression = orig.ident == "v3.1k") )
```

![](Quality_control_files/figure-gfm/scatter2-1.png)<!-- -->

### Plot top expressed genes

There is currently no function in Seurat to plot the top expressed genes
as both Scater and Scanpy does, so here is some code for a custom plot.

``` r
total_genes = Matrix::rowSums(C)
top = order(total_genes, decreasing = T)[1:20]
c.top = C[top,]
# normalize by total counts to get fraction of umis.
n = Matrix::colSums(C)
c.norm = sweep(c.top, 2, n, FUN="/")
c.norm <- data.frame(c.norm)
c.norm <- melt(t(c.norm))


p <- ggplot(c.norm, aes(x=Var2, y=value)) + 
  geom_boxplot() + coord_flip()
print(p)
```

![](Quality_control_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

As you can see, MALAT1 consitutes up to 30% of the umis from a single
cell and the other top genes are mitochondrial.

### Filtering

#### Mitochondrial filtering

We have quite a lot of cells with high proportion of mitochondrial
reads. It could be wise to remove those cells, if we have enough cells
left after filtering. Another option would be to either remove all
mitochondrial reads from the dataset and hope that the remaining genes
still have enough biological signal. A third option would be to just
regress out the percent.mito variable during scaling.

In this case we have as much as 99.7% mitochondrial reads in some of the
cells, so it is quite unlikely that there is much celltype signature
left in those.

Looking at the plots, make resonable decisions on where to draw the
cutoff. In this case, the bulk of the cells are below 25% mitochondrial
reads and that will be used as a cutoff.

``` r
#select cells with percent.mito < 25
selected <- WhichCells(alldata, expression = percent.mito < 25)
length(selected)
```

    ## [1] 2703

``` r
# and subset the object to only keep those cells
data.filt <- subset(alldata, cells = selected)

# plot violins for new data
VlnPlot(data.filt, features = "percent.mito")
```

![](Quality_control_files/figure-gfm/mito.filt-1.png)<!-- -->

As you can see, there is still quite a lot of variation in percent mito,
so it will have to be dealt with in the data analysis step.

#### Gene detection filtering

Extremely high number of detected genes could indicate doublets.
However, depending on the celltype composition in your sample, you may
have cells with higher number of genes (and also higher counts) from one
celltype.

In these datasets, there is also a clear difference between the v2 vs v3
10x chemistry with regards to gene detection, so it may not be fair to
apply the same cutoffs to all of them.

Also, in the protein assay data there is a lot of cells with few
detected genes giving a bimodal distribution. This type of distribution
is not seen in the other 2 datasets. Considering that they are all pbmc
datasets it makes sense to regard this distribution as low quality
libraries.

Filter the cells with high gene detection (putative doublets) with
cutoffs 4100 for v3 chemistry and 2000 for v2.

``` r
#start with cells with many genes detected.
high.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA > 4100)
high.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA > 2000 & orig.ident == "v2.1k")

# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(high.det.v2,high.det.v3)))

# check number of cells
ncol(data.filt)
```

    ## [1] 2631

Filter the cells with low gene detection (low quality libraries) with
less than 1000 genes for v2 and \< 500 for v2.

``` r
#start with cells with many genes detected.
low.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA < 1000 & orig.ident != "v2.1k")
low.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA < 500 & orig.ident == "v2.1k")

# remove these cells
data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(low.det.v2,low.det.v3)))

# check number of cells
ncol(data.filt)
```

    ## [1] 2531

#### Plot QC-stats again

Lets plot the same qc-stats another time.

``` r
VlnPlot(data.filt, features = "nFeature_RNA", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot2-1.png)<!-- -->

``` r
VlnPlot(data.filt, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot2-2.png)<!-- -->

``` r
VlnPlot(data.filt, features = "percent.mito", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot2-3.png)<!-- -->

``` r
VlnPlot(data.filt, features = "percent.ribo", pt.size = 0.1) + NoLegend()
```

![](Quality_control_files/figure-gfm/vln.plot2-4.png)<!-- -->

``` r
# and check the number of cells per sample before and after filtering
table(Idents(alldata))
```

    ## 
    ## p3.1k v2.1k v3.1k 
    ##   713   996  1222

``` r
table(Idents(data.filt))
```

    ## 
    ## p3.1k v2.1k v3.1k 
    ##   526   933  1072

### Calculate cell-cycle scores

Seurat has a function for calculating cell cycle scores based on a list
of know S-phase and G2/M-phase genes.

``` r
data.filt <- CellCycleScoring(
  object = data.filt,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

VlnPlot(data.filt, features = c("S.Score","G2M.Score"))
```

![](Quality_control_files/figure-gfm/cc-1.png)<!-- -->

In this case it looks like we only have a few cycling cells in the
datasets.

### Save the filtered data

For coming analyses, remove the mitochondiral genes and MALAT1 from the
matrix and save to a file.

``` r
keep = setdiff(rownames(data.filt),c(mt.genes, "MALAT1"))
data.filt <- data.filt[keep,]

savefile = "../write/seurat/filtered_3pbmc.Rdata"
save(data.filt, file=savefile)
```

### Session info

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS  10.14.6
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /Users/asbj/Programs/miniconda3_4.2.12/envs/elixir-course/lib/R/lib/libRblas.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] reshape2_1.4.3              Matrix_1.2-17              
    ##  [3] scater_1.10.1               ggplot2_3.1.1              
    ##  [5] SingleCellExperiment_1.4.0  SummarizedExperiment_1.12.0
    ##  [7] DelayedArray_0.8.0          BiocParallel_1.16.6        
    ##  [9] matrixStats_0.54.0          Biobase_2.42.0             
    ## [11] GenomicRanges_1.34.0        GenomeInfoDb_1.18.1        
    ## [13] IRanges_2.16.0              S4Vectors_0.20.1           
    ## [15] BiocGenerics_0.28.0         Seurat_3.0.0               
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rtsne_0.15               ggbeeswarm_0.6.0        
    ##   [3] colorspace_1.4-1         ggridges_0.5.1          
    ##   [5] XVector_0.22.0           listenv_0.7.0           
    ##   [7] npsurv_0.4-0             bit64_0.9-7             
    ##   [9] ggrepel_0.8.1            codetools_0.2-16        
    ##  [11] splines_3.5.1            R.methodsS3_1.7.1       
    ##  [13] lsei_1.2-0               knitr_1.20              
    ##  [15] jsonlite_1.6             ica_1.0-2               
    ##  [17] cluster_2.0.9            png_0.1-7               
    ##  [19] R.oo_1.22.0              HDF5Array_1.10.1        
    ##  [21] sctransform_0.2.0        compiler_3.5.1          
    ##  [23] httr_1.4.0               assertthat_0.2.1        
    ##  [25] lazyeval_0.2.2           htmltools_0.3.6         
    ##  [27] tools_3.5.1              rsvd_1.0.0              
    ##  [29] igraph_1.2.4.1           gtable_0.3.0            
    ##  [31] glue_1.3.1               GenomeInfoDbData_1.2.1  
    ##  [33] RANN_2.6                 dplyr_0.8.0.1           
    ##  [35] Rcpp_1.0.1               gdata_2.18.0            
    ##  [37] ape_5.3                  nlme_3.1-139            
    ##  [39] DelayedMatrixStats_1.4.0 gbRd_0.4-11             
    ##  [41] lmtest_0.9-36            stringr_1.4.0           
    ##  [43] globals_0.12.4           irlba_2.3.3             
    ##  [45] gtools_3.8.1             future_1.12.0           
    ##  [47] MASS_7.3-51.4            zlibbioc_1.28.0         
    ##  [49] zoo_1.8-5                scales_1.0.0            
    ##  [51] rhdf5_2.26.2             RColorBrewer_1.1-2      
    ##  [53] yaml_2.2.0               reticulate_1.12         
    ##  [55] pbapply_1.4-0            gridExtra_2.3           
    ##  [57] stringi_1.4.3            caTools_1.17.1.2        
    ##  [59] bibtex_0.4.2             Rdpack_0.10-1           
    ##  [61] SDMTools_1.1-221.1       rlang_0.3.4             
    ##  [63] pkgconfig_2.0.2          bitops_1.0-6            
    ##  [65] evaluate_0.13            lattice_0.20-38         
    ##  [67] ROCR_1.0-7               purrr_0.3.2             
    ##  [69] Rhdf5lib_1.4.3           labeling_0.3            
    ##  [71] htmlwidgets_1.3          bit_1.1-14              
    ##  [73] cowplot_0.9.4            tidyselect_0.2.5        
    ##  [75] plyr_1.8.4               magrittr_1.5            
    ##  [77] R6_2.4.0                 gplots_3.0.1.1          
    ##  [79] pillar_1.3.1             withr_2.1.2             
    ##  [81] fitdistrplus_1.0-14      survival_2.44-1.1       
    ##  [83] RCurl_1.95-4.12          tibble_2.1.1            
    ##  [85] future.apply_1.1.0       tsne_0.1-3              
    ##  [87] crayon_1.3.4             hdf5r_1.2.0             
    ##  [89] KernSmooth_2.23-15       plotly_4.8.0            
    ##  [91] rmarkdown_1.11           viridis_0.5.1           
    ##  [93] grid_3.5.1               data.table_1.12.2       
    ##  [95] metap_1.1                digest_0.6.18           
    ##  [97] tidyr_0.8.3              R.utils_2.8.0           
    ##  [99] munsell_0.5.0            beeswarm_0.2.3          
    ## [101] viridisLite_0.3.0        vipor_0.4.5
