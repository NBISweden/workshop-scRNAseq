---
layout: default
title: Exercises - scRNAseq course
---

# Welcome to NBIS scRNA-seq tutorials

This page contains links to different tutorials that are used in the scRNA-seq analysis course. We see this as a "smörgåsbord" of tutorials, we do not expect that all course participants will have time to run through all of the tutorials, so we suggest that you pick-and-mix the analyses/packages that you find most interesting. 

For many of the packages we recommend that you follow the tutorials supplied by the different packages. We have included some examples on how to run these tutorials based on your own dataset and also some suggestions on which parts of the tutorials to focus on. 

We have run the tutorials with SmartSeq2 data from:

* Mouse embryonic development from [Deng *et al.*](http://science.sciencemag.org/content/343/6167/193.long) 
* Human innate lymphoid cells from [Björklund *et al.*](https://www.nature.com/articles/ni.3368)

All raw data, and some intermediate files that take long time to compute, are available at uppmax folder: `/proj/uppstore2017171/courses/scrnaseq_course/data/`

You should be able to run these tutorials with your own data. But be aware, some steps of the tutorials will take long to run with large dataset, so it may be a good idea to use several R-sessions in paralell and work on other stuff while you are waiting. 

As you run into problems, we will try to fill in the [FAQ](FAQ) with common quiestons.

If you want to load all the code into R directly, you can access the R-markdown documents at our [github site](https://github.com/NBISweden/workshop-scRNAseq/tree/master/labs)

#### Pipeline for mapping reads, QC and expression estimates

Snakemake pipeline for processing SmartSeq2 data.

*	[Pipeline tutorial](labs/Pipeline_exercise) 

#### Biomart - Convert gene names and fetch annotation

For those not familiar with working with biomaRt, we suggest that you have a look at this example code for how to convert between different formats using biomaRt. 
 
*	[Tutorial for biomaRt](labs/biomart) 

#### PCA, tSNE and clustering

Basic PCA, tSNE and clustering using base R on mouse embryonic development data.

*	[Tutorial for PCA and clustering](labs/PCA_and_clustering)


#### Constructing knn graphs

Construction of graphs from cell-cell similiarities using igraph.  

*       [Tutorial for KNN graphs](labs/igraph)		  

#### scater package

Tutorial with the scater package for QC of scRNA-seq data

*	[Tutorial for scater](labs/scater_ilc)


#### Estimating Batch-Effects

A tutorial for estimating genome-wide and individual genes batch-effects.

*       [Tutorial for Batch-Effects](labs/batch_analysis)

#### Comparing Normalization Methods

A tutorial for comparison scRNAseq and bulk RNAseq normalization strategies.

*       [Tutorial for Normalization](labs/norm_analysis)


#### SC3 package

Tutorial with the SC3 consensus clustering package

*	[Tutorial for SC3](labs/sc3_ilc)

#### Pagoda package

Pagoda patway wPCA for clustering of cells. OBS! several steps in this tutorial takes hours to run if you work with your own dataset, a good suggestion is to start with the first steps, knn.error.model, pagoda.varnorm and pagoda.pathway.wPCA and let it run while working on other tutorials. You can also run it with more than one core to speed things up.
 
*	[Tutorial for Pagoda](labs/pagoda_ilc)

#### Seurat package

Tutorial for Seurat package with normalization, dimensionality reduction and clustering.

*	[Tutorial for Seurat](labs/seurat_analysis)

#### Differential expression

For this tutorial we have included several different methods for differential expression tests on single cell data, including SCDE, MAST, SC3 and Seurat. The exercise has been split into 2 parts with evaluation of all results in the second part. 

*	[Tutorial for DE detection 1](labs/Differential_gene_expression)
*	[Tutorial for DE detection 2](labs/Differential_gene_expression_part2)

#### Monocle Pseudotime analysis

A tutorial with mouse embryonic data using the Monocle package for pseudotime analysis.

*	[Tutorial for Monocle](labs/monocle_analysis)   

#### UPPMAX Sbatch
 
One example of a sbatch script
 
*	[sbatch scripts](labs/sbatchScript)   
 
#### Caveat

We will try to keep these tutorials up to date. If you find any errors or things that you think should be updated please contact Asa (asa.bjorklund@scilifelab.se) 
  		
