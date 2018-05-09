Methods for differential gene expression in scRNAseq data
=========================================================

Here are some examples of methods commonly used for differential expression in scRNAseq data. Included in this tutorial are:

-   SCDE
-   MAST
-   SC3 package - Kruskall-Wallis test
-   Pagoda package - 4 different tests

The dataset used is single-cell RNA-seq data (SmartSeq) from mouse embryonic development from Deng. et al. Science 2014, Vol. 343 no. 6167 pp. 193-196, "Single-Cell RNA-Seq Reveals Dynamic, Random Monoallelic Gene Expression in Mammalian Cells".

We will check for differentially expressed gene between 8-cell and 16-cell stage embryos, it is quite a small dataset, for faster calculations.

All data you need is available in the course uppmax folder with subfolder:

`scrnaseq_course/data/mouse_embryo/`

Load data
---------

First read in the data.

``` r
# we will read both counts and rpkms as different method are more adopted for different type of data
R<-read.table("data/mouse_embryo/rpkms_Deng2014_preinplantation.txt")
C<-read.table("data/mouse_embryo/counts_Deng2014_preinplantation.txt")
M <- read.table("data/mouse_embryo/Metadata_mouse_embryo.csv",sep=",",header=T)

# select only 8 and 16 stage embryos.
selection <- c(grep("8cell",M$Stage),grep("16cell",M$Stage))

# check number of cells
table(M$Stage[selection])
```

    ## 
    ## early2cell earlyblast  late2cell  lateblast   mid2cell   midblast 
    ##          0          0          0          0          0          0 
    ##  MIIoocyte    X16cell     X4cell     X8cell         zy 
    ##          0         50          0         28          0

``` r
# select those cells only from the data frames
M<-M[selection,]
R <- R[,selection]
C <- C[,selection]
```

SCDE
====

We will first run DE with SCDE - see more at: <http://hms-dbmi.github.io/scde/diffexp.html>

``` r
require(scde)
```

    ## Warning: package 'scde' was built under R version 3.4.2

``` r
# make a count dataset and clean up
cd <- clean.counts(C, min.lib.size=1000, min.reads = 1, min.detected = 1)

# make factor for stage
stages <- factor(M$Stage,levels=c("X8cell","X16cell"))
names(stages) <- colnames(C)

# fit error models - takes a while to run (~20 mins).
# you can skip this step and load the precomputed file.
savefile<-"data/mouse_embryo/DE/scde_error_models.Rdata"
if (file.exists(savefile)){
  load(savefile)
}else{
  o.ifm <- scde.error.models(counts = cd, groups = stages, n.cores = 1, 
                             threshold.segmentation = TRUE, save.crossfit.plots = FALSE, 
                             save.model.plots = FALSE, verbose = 0)
  save(o.ifm,file=savefile)
}

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)

# run differential expression test
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  stages, 
                                    n.randomizations  =  100, n.cores  =  1, verbose  =  1)
```

    ## comparing groups:
    ## 
    ## X16cell  X8cell 
    ##      50      28 
    ## calculating difference posterior
    ## summarizing differences

``` r
# convert Z-score and corrected Z-score to 2-sided p-values 
ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
ediff$p.adjust <- 2*pnorm(-abs(ediff$cZ))

# sort and look at top results
ediff <- ediff[order(abs(ediff$Z), decreasing  =  TRUE), ]
head(ediff)
```

    ##                   lb        mle        ub         ce         Z        cZ
    ## Actg1     -1.9809972 -1.4707403 -1.020514 -1.0205137 -7.033208 -5.518869
    ## Elf3      -3.6618432 -2.4612389 -1.530771 -1.5307705 -6.855648 -5.418799
    ## Odc1       0.5702871  0.9004533  1.170589  0.5702871  6.246663  4.719634
    ## Dab2      -3.4217224 -2.3411785 -1.320665 -1.3206648 -6.077668 -4.558448
    ## Socs3      0.9904986  1.5908008  2.191103  0.9904986  5.979057  4.476501
    ## Serpinb9b -4.5022663 -3.7218735 -2.431224 -2.4312238 -5.371590 -3.697397
    ##                 pvalue     p.adjust
    ## Actg1     2.018385e-12 3.411879e-08
    ## Elf3      7.098988e-12 6.000065e-08
    ## Odc1      4.193144e-10 2.362697e-06
    ## Dab2      1.219431e-09 5.153314e-06
    ## Socs3     2.244328e-09 7.587623e-06
    ## Serpinb9b 7.804546e-08 2.178212e-04

``` r
# write output to a file with top DE genes on top.
write.table(ediff,file="data/mouse_embryo/DE/scde_8cell_vs_16_cell.tab",sep="\t",quote=F)
```

If you want to convert SCDE Z-scores into p-values use the pnorm function in R. More info on Z-scores and p-values can be found here: <http://pro.arcgis.com/en/pro-app/tool-reference/spatial-statistics/what-is-a-z-score-what-is-a-p-value.htm>

The `mle` column (maximum likelihood estimate) corresponds to the log fold-change. Positive/negative mle's and Z-scores will indicate if a gene is up/down-regulated.

### SCDE with Batch info

``` r
# include also batch information in the test
# in this case we will consider each embryo as a batch just for testing
# OBS! in this case there are no batch that belongs to both groups so the correction may be pointless.

batch <- factor(M$Embryo, levels <- unique(M$Embryo))
ediff.batch <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  stages, 
                                          n.randomizations  =  100, n.cores  =  1, batch=batch)
```

    ## WARNING: strong interaction between groups and batches! Correction may be ineffective:
    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  bgti
    ## p-value < 2.2e-16
    ## alternative hypothesis: two.sided

``` r
# now scde.expression.difference returns a list with the corrected values as well as the DE results.

de.batch <- ediff.batch$results
de.batch$pvalue <- 2*pnorm(-abs(de.batch$Z))
de.batch$p.adjust <- 2*pnorm(-abs(de.batch$cZ))

# look at top results
head(de.batch[order(abs(de.batch$Z), decreasing  =  TRUE), ])
```

    ##                   lb        mle        ub         ce         Z        cZ
    ## Actg1     -1.9809972 -1.4707403 -1.020514 -1.0205137 -7.033208 -5.518869
    ## Elf3      -3.6618432 -2.4612389 -1.530771 -1.5307705 -6.855648 -5.418799
    ## Odc1       0.5702871  0.9004533  1.170589  0.5702871  6.246663  4.719634
    ## Dab2      -3.4217224 -2.3411785 -1.320665 -1.3206648 -6.077668 -4.558448
    ## Socs3      0.9904986  1.5908008  2.191103  0.9904986  5.979057  4.476501
    ## Serpinb9b -4.5022663 -3.7218735 -2.431224 -2.4312238 -5.371590 -3.697397
    ##                 pvalue     p.adjust
    ## Actg1     2.018385e-12 3.411879e-08
    ## Elf3      7.098988e-12 6.000065e-08
    ## Odc1      4.193144e-10 2.362697e-06
    ## Dab2      1.219431e-09 5.153314e-06
    ## Socs3     2.244328e-09 7.587623e-06
    ## Serpinb9b 7.804546e-08 2.178212e-04

In this case we get the same results after batch correction, but when you have clear batch effects that may influence your data, you should include a batch information to the test.

MAST
====

Running the Hurdle model from MAST for DE. For more info, see: <http://bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html>

``` r
library(MAST)

fdata <- data.frame(rownames(R))
rownames(fdata)<-rownames(R)

# Make a single cell assay object
sca <- FromMatrix(log2(as.matrix(R)+1),M,fdata)

# count number of detected genes and scale.
cdr2 <-colSums(assay(sca)>0)
colData(sca)$cngeneson <- scale(cdr2)

# Fit a hurdle model, modeling the condition and (centered) ngeneson factor, thus adjusting for 
# the cellular detection rate. In this case we set the reference to be 8cell stage using relevel 
# of the factor.

# takes a while to run, so save to an file so that you do not have to rerun each time
savefile<-"data/mouse_embryo/DE/mast_data.Rdata"
if (file.exists(savefile)){
  load(savefile)
}else{  
  cond<-factor(colData(sca)$Stage)
  cond <- relevel(cond,"X8cell")
  colData(sca)$condition<-cond
  zlmCond <- zlm(~condition + cngeneson, sca)
  #Run likelihood ratio test for the condition coefficient.
  summaryCond <- summary(zlmCond, doLRT='conditionX16cell') 
  save(sca,zlmCond,summaryCond,file=savefile)
}
# get the datatable with all outputs
summaryDt <- summaryCond$datatable

# Make a datatable with all the relevant output - combined Hurdle model
fcHurdle <- merge(summaryDt[contrast=='conditionX16cell' & component=='H',.(primerid, `Pr(>Chisq)`)], 
     summaryDt[contrast=='conditionX16cell' & component=='logFC', 
               .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

# change name of Pr(>Chisq) to fdr and sort by fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdle <- fcHurdle[order(fdr),]
head(fcHurdle)
```

    ##    primerid   Pr(>Chisq)      coef      ci.hi     ci.lo          fdr
    ## 1:     Odc1 2.922235e-09 -1.049094 -0.7220156 -1.376173 6.708867e-05
    ## 2:     Elf3 5.232456e-08  3.590321  4.8000744  2.380568 6.006336e-04
    ## 3:    Fbxo3 1.703398e-07  2.678172  3.4682607  1.888083 9.776652e-04
    ## 4:  Gm15645 1.579087e-07 -1.120809 -0.7641414 -1.477477 9.776652e-04
    ## 5:   Ly6g6e 2.610798e-07 -3.208303 -2.1856130 -4.230992 1.198774e-03
    ## 6:    Apoa1 6.037716e-07  2.789231  3.8658661  1.712596 2.310231e-03

``` r
# write output to a table 
write.table(fcHurdle,file="data/mouse_embryo/DE/mast_8cell_vs_16_cell.tab",sep="\t",quote=F)
```

SC3 Kruskall-Wallis test
========================

In the SC3 package they have implemented non-parametric Kruskal-Wallis test for differential expression.

``` r
library(SC3)

# run their DE test using rpkm values and as labels, M$stage
de <- get_de_genes(R,M$Stage)

# output is simply a p-value with no information of directionality.
de.results <- data.frame(gene=rownames(R),p.value=de)
de.results <- de.results[order(de.results$p.value),]

head(de.results)
```

    ##            gene      p.value
    ## Odc1       Odc1 1.312182e-06
    ## Gm15645 Gm15645 3.052558e-05
    ## Socs3     Socs3 1.893245e-04
    ## Ly6g6e   Ly6g6e 6.820155e-04
    ## Ybx1       Ybx1 8.014639e-04
    ## Klf17     Klf17 1.149117e-03

``` r
write.table(de.results,file="data/mouse_embryo/DE/sc3_kwtest_8cell_vs_16_cell.tab",sep="\t",quote=F)
```

Since many different packages have been loaded, you should continue to [Part 2](Differential_gene_expression_part2) in a new R-session so that you can load Seurat unless you have increased the number of DLLs you are using.

##### Session info

``` r
sessionInfo()
```

    ## R version 3.4.1 (2017-06-30)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Sierra 10.12.6
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] SC3_1.7.7                  MAST_1.4.1                
    ##  [3] SummarizedExperiment_1.8.1 DelayedArray_0.4.1        
    ##  [5] matrixStats_0.53.0         Biobase_2.38.0            
    ##  [7] GenomicRanges_1.30.1       GenomeInfoDb_1.14.0       
    ##  [9] IRanges_2.12.0             S4Vectors_0.16.0          
    ## [11] BiocGenerics_0.24.0        scde_2.6.0                
    ## [13] flexmix_2.3-14             lattice_0.20-35           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-131               bitops_1.0-6              
    ##  [3] pbkrtest_0.4-7             doParallel_1.0.11         
    ##  [5] RColorBrewer_1.1-2         rprojroot_1.3-2           
    ##  [7] tools_3.4.1                backports_1.1.2           
    ##  [9] doRNG_1.6.6                R6_2.2.2                  
    ## [11] KernSmooth_2.23-15         lazyeval_0.2.1            
    ## [13] mgcv_1.8-23                colorspace_1.3-2          
    ## [15] nnet_7.3-12                compiler_3.4.1            
    ## [17] quantreg_5.34              Cairo_1.5-9               
    ## [19] SparseM_1.77               pkgmaker_0.22             
    ## [21] caTools_1.17.1             scales_0.5.0              
    ## [23] mvtnorm_1.0-7              DEoptimR_1.0-8            
    ## [25] robustbase_0.92-8          Lmoments_1.2-3            
    ## [27] stringr_1.2.0              digest_0.6.15             
    ## [29] minqa_1.2.4                rmarkdown_1.8             
    ## [31] distillery_1.0-4           XVector_0.18.0            
    ## [33] rrcov_1.4-3                htmltools_0.3.6           
    ## [35] lme4_1.1-15                WriteXLS_4.0.0            
    ## [37] extRemes_2.0-8             limma_3.34.8              
    ## [39] rlang_0.1.6                shiny_1.0.5               
    ## [41] BiocParallel_1.12.0        gtools_3.5.0              
    ## [43] car_2.1-6                  RCurl_1.95-4.10           
    ## [45] magrittr_1.5               modeltools_0.2-21         
    ## [47] GenomeInfoDbData_1.0.0     Matrix_1.2-12             
    ## [49] Rcpp_0.12.15               munsell_0.4.3             
    ## [51] abind_1.4-5                stringi_1.1.6             
    ## [53] yaml_2.1.16                edgeR_3.20.8              
    ## [55] MASS_7.3-48                zlibbioc_1.24.0           
    ## [57] gplots_3.0.1               plyr_1.8.4                
    ## [59] grid_3.4.1                 gdata_2.18.0              
    ## [61] splines_3.4.1              locfit_1.5-9.1            
    ## [63] knitr_1.19                 pillar_1.1.0              
    ## [65] rjson_0.2.15               rngtools_1.2.4            
    ## [67] reshape2_1.4.3             codetools_0.2-15          
    ## [69] evaluate_0.10.1            RcppArmadillo_0.8.300.1.0 
    ## [71] pcaMethods_1.70.0          data.table_1.10.4-3       
    ## [73] httpuv_1.3.5               nloptr_1.0.4              
    ## [75] foreach_1.4.4              MatrixModels_0.4-1        
    ## [77] gtable_0.2.0               RMTstat_0.3               
    ## [79] ggplot2_2.2.1              mime_0.5                  
    ## [81] xtable_1.8-2               e1071_1.6-8               
    ## [83] pcaPP_1.9-73               class_7.3-14              
    ## [85] pheatmap_1.0.8             SingleCellExperiment_1.0.0
    ## [87] tibble_1.4.2               iterators_1.0.9           
    ## [89] registry_0.5               cluster_2.0.6             
    ## [91] Rook_1.1-1                 brew_1.0-6                
    ## [93] ROCR_1.0-7
