---
layout: default
title: Exercises - scRNAseq course
---

##### <img border="0" src="https://www.svgrepo.com/show/6672/exercise.svg" width="40" height="40"> Exercises
***

Here we provide short tutorials on the different steps of scRNAseq analysis using either of the 3 commonly used scRNAseq analysis pipelines, [Seurat](https://satijalab.org/seurat/), [Scran](https://bioconductor.org/packages/release/bioc/html/scran.html) and [Scanpy](https://scanpy.readthedocs.io/en/stable/). It is up to you which one you want to try out, if you finish quickly, you may have time to run several of them or run of the additional labs below. In principle we perform the same steps with all 3 pipelines, but there are some small differences as all different methods are not implemented in all the pipelines.

<br/>

Please me sure you have completed the [<img border="0" src="https://www.svgrepo.com/show/19652/maths-class-materials-cross-of-a-pencil-and-a-ruler.svg" width="20" height="20">**Precourse material**](precourse.md)

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/6672/exercise.svg" width="40" height="40"> MAIN exercises
***



| Tutorial | <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/1448px-R_logo.svg.png" width="20" height="20"> Seurat | <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/1448px-R_logo.svg.png" width="20" height="20"> Scater/Scran | <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Python-logo-notext.svg/1024px-Python-logo-notext.svg.png" width="20" height="20"> Scanpy |
| -------- | ---------- | ---------------- | --------------- |
| QC | [Seurat_qc](labs/compiled/seurat/seurat_01_qc_compiled.md) | [Scater_qc](labs/compiled/scater/scater_01_qc_compiled.md) | [ScanPY_qc](labs/scanpy/qc_3pbmc.ipynb) |
| Dimensionality reduction | [Seurat_dr](labs/compiled/seurat/seurat_02_dim_reduction_compiled.md) | [Scater_dr](labs/compiled/scater/scater_02_dim_reduction_compiled.md) | [ScanPY_dr](labs/scanpy/dim_reduction.ipynb) |
| Data integration | [Seurat_integr](labs/compiled/seurat/seurat_03_integration_compiled.md) | [Scater_integr](labs/compiled/scater/scater_03_integration_compiled.md) | [ScanPY_integr](labs/scanpy/batch_correction_mnn.ipynb) |
| Clustering | [Seurat_clust](labs/compiled/seurat/lab_seurat.html) | [Scater_clust](labs/compiled/scater/lab_scran.html) | [ScanPY_clust](labs/scanpy/qc_3pbmc.ipynb) |
| Differential expression | [Seurat_dge](labs/compiled/seurat/lab_seurat.html) | [Scater_dge](labs/compiled/scater/lab_scran.html) | [ScanPY_dge](labs/scanpy/qc_3pbmc.ipynb) |
| Trajectory inference | [Monocle_ti](labs/compiled/monocle/monocle.html) | [Slingshot_ti](labs/compiled/slingshot/slingshot.html) | [PAGA_ti](labs/paga/paga.ipynb) |

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/48895/exercise.svg" width="40" height="40"> BONUS exercises
***

| Name (link) | Description |
| ----------- | ----------- |
| [Read-Pipeline](labs/Pipeline_exercise) | Snakemake pipeline for processing SmartSeq2 data, mapping reads, QC and expression estimates|
| [biomaRt](labs/biomart) | For those not familiar with working with biomaRt, we suggest that you have a look at this example code for how to convert between different formats using biomaRt|
| [PCA, tSNE and clustering](labs/PCA_and_clustering) | Basic PCA, tSNE and clustering using base R on mouse embryonic development data. | 
| [KNN-graphs](labs/igraph) | Construction of graphs from cell-cell similiarities using igraph|
| [Estimating Batch-Effects](https://bitbucket.org/scilifelab-lts/scrnaseq-labs/src/a228442debe7f8eff28cfdba875349025db9b7a3/batch_analysis.md?fileviewer=file-view-default) | A tutorial for estimating genome-wide and individual genes batch-effects |
| [Normalization Comparison](labs/norm_analysis_v2) | A tutorial for comparison scRNAseq and bulk RNAseq normalisation strategies. | [Tutorial for Normalisation](labs/norm_analysis_v2)  |
| [SC3 package](labs/sc3_R35) | Tutorial with the SC3 consensus clustering package |
| [Trajectory with Monocle2](labs/monocle_analysis) | A tutorial with mouse embryonic data using the Monocle package for pseudotime analysis |
| [Differential expression_2](labs/Differential_gene_expression) [Differential expression_2](labs/Differential_gene_expression) | OBS! This old tutorial uses Seurat v2! For this tutorial we have included several different methods for differential expression tests on single cell data, including SCDE, MAST, SC3 and Seurat. The exercise has been split into 2 parts with evaluation of all results in the second part |
| [UPPMAX Sbatch](labs/sbatchScript) | One example of a sbatch script |
| [Pagoda](labs/pagoda_ilc) | Pagoda patway wPCA for clustering of cells. OBS! several steps in this tutorial takes hours to run if you work with your own dataset, a good suggestion is to start with the first steps, knn.error.model, pagoda.varnorm and pagoda.pathway.wPCA and let it run while working on other tutorials. You can also run it with more than one core to speed things up |


We will try to keep these tutorials up to date. If you find any errors or things that you think should be updated please contact Åsa (asa.bjorklund@scilifelab.se) 
  	
	
<br/>

##### <img border="0" src="https://www.svgrepo.com/show/83019/faq-button.svg" width="40" height="40"> FAQ
***

As you run into problems, we will try to fill in the [FAQ](labs/FAQ) with common questions.

<br/>

	
