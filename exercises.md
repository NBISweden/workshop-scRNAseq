# Welcome to NBIS scRNA-seq tutorial packages

This page contains links to different tutorials that are used in the scRNA-seq analysis course. We see this as a "smörgåsbord" of tutorials, we do not expect that all course participants will have time to run through all of the tutorials, so we suggest that you pick-and-mix the analyses/packages that you find most interesting. 

For many of the packages we recommend that you follow the tutorials supplied by the different packages. But if you have brought your own data to the tutorials, we have included some examples on how to run these tutorials based on your own dataset and also some suggestions on which parts of the tutorials to focus on.

We have run the tutorials with SmartSeq2 data from mouse embryonic development and human innate lymphoid cells. All data and intermediate files are available at uppmax folder: `/proj/b2013006/nobackup/scrnaseq_course/data/`

But be aware, some steps of the tutorials will take long to run with large dataset, so it may be a good idea to use several R-sessions in paralell and work on other stuff while you are waiting. 

As you run into problems, we will try to fill in the [FAQ](FAQ.md) with common quiestons.

If you want to load all the code into R directly, you can access the R-markdown documents by clicking on the `Source` link to the left, or by cloning the whole repository.

### Pipeline for mapping reads, QC and expression estimates

Snakemake pipeline for processing SmartSeq2 data.

*	[Pipeline tutorial](Pipeline_exercise.md) 

### Converting gene names to Ensembl IDs

For those not familiar with working with biomaRt, we suggest that you have a look at this example code for how to convert between different formats using biomaRt. 
 
*	[Tutorial for biomaRt](biomart.md) 
### PCA and clustering

Basic PCA and clustering using base R on mouse embryonic development data.

*	[Tutorial for PCA and clustering](exercises/PCA_and_clustering.md)

### scater package

Tutorial with the scater package for QC of scRNA-seq data

*	[Tutorial for scater](scater_ilc.md)

### Batch analysis

En example on how to check for batch effect and identify genes driven by the batch effect

*	[Tutorial for Batch analysis](batch_analysis.md)

### Normalization

A tutorial showing different methods for data normalizations, mainly using Scran package.

*	[Tutorial for Data normalization](norm_analysis.md)

### SC3 package

Tutorial with the SC3 consensus clustering package

*	[Tutorial for SC3](sc3_ilc.md)

### Pagoda package

Pagoda patway wPCA for clustering of cells. OBS! several steps in this tutorial takes hours to run if you work with your own dataset, a good suggestion is to start with the first steps, knn.error.model, pagoda.varnorm and pagoda.pathway.wPCA and let it run while working on other tutorials. You can also run it with more than one core to speed things up.
 
*	[Tutorial for Pagoda](pagoda_ilc.md)

### Seurat package

Tutorial for Seurat package with normalization, dimensionality reduction and clustering.

*	[Tutorial for Seurat](seurat_analysis.md)

### Differential expression

For this tutorial we have included several different methods for differential expression tests on single cell data, including SCDE, MAST, SC3 and Seurat. The exercise has been split into 2 parts with evaluation of all results in the second part. 

*	[Tutorial for DE detection 1](Differential_gene_expression.md)
*	[Tutorial for DE detection 2](Differential_gene_expression_part2.md)

### Monocle Pseudotime analysis

A tutorial with mouse embryonic data using the Monocle package for pseudotime analysis.

*	[Tutorial for Monocle](monocle_analysis.md)   

### Estimating Batch-Effects

A tutorial for estimating genome-wide and individual genes batch-effects.

*	[Tutorial for Batch-Effects](batch_analysis.md)   

### Comparing Normalization Methods

A tutorial for comparison scRNAseq and bulk RNAseq normalization strategies.

*	[Tutorial for Normalization](norm_analysis.md)   

### UPPMAX Sbatch
 
One example of a sbatch script
 
*	[sbatch scripts](sbatchScript.md)   
 
## Caveat

We will try to keep these tutorials up to date. If you find any errors or things that you think should be updated please contact Asa (asa.bjorklund@scilifelab.se) 
  		
