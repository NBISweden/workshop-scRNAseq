---
layout: default
title: Exercises - scRNAseq course
---

#### <img border="0" src="https://www.svgrepo.com/show/6672/exercise.svg" width="50" height="50"> Exercises
***

Here we provide short tutorials on the different steps of scRNAseq analysis using either of the 3 commonly used scRNAseq analysis pipelines, [Seurat](https://satijalab.org/seurat/), [Scran](https://bioconductor.org/packages/release/bioc/html/scran.html) and [Scanpy](https://scanpy.readthedocs.io/en/stable/). It is up to you which one you want to try out, if you finish quickly, you may have time to run several of them or run of the additional labs below. In principle we perform the same steps with all 3 pipelines, but there are some small differences as all different methods are not implemented in all the pipelines.

<br/>


##### <img border="0" src="https://www.svgrepo.com/show/6672/exercise.svg" width="40" height="40"> MAIN exercises
***

All scripts (Rmarkdown or ipython notebooks) can be found at our github repo in folder [labs/compiled](https://github.com/NBISweden/workshop-scRNAseq/tree/master/labs/compiled)

Please me sure you have completed the [<img border="0" src="https://www.svgrepo.com/show/19652/maths-class-materials-cross-of-a-pencil-and-a-ruler.svg" width="20" height="20">**Precourse material**](precourse.md) as well as the [<img border="0" src="https://hackernoon.com/hn-images/1*rW03Wtue71AKfxnx6XN_iQ.png" width="20" height="20">**Conda instructions**](conda_instructions.md)

| Tutorial | <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/1448px-R_logo.svg.png" width="20" height="20"> Seurat | <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/R_logo.svg/1448px-R_logo.svg.png" width="20" height="20"> Scater/Scran | <img border="0" src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/c3/Python-logo-notext.svg/1024px-Python-logo-notext.svg.png" width="20" height="20"> Scanpy |
| -------- | ---------- | ---------------- | --------------- |
| <img border="0" src="https://cdn0.iconfinder.com/data/icons/business-and-finance-4-5/68/188-512.png" width="30" height="30"> Quality Control | [Seurat_qc](labs/compiled/seurat/seurat_01_qc.md) | [Scater_qc](labs/compiled/scater/scater_01_qc.md) | [ScanPY_qc](labs/compiled/scanpy/scanpy_01_qc.html) |
| <img border="0" src="https://static.thenounproject.com/png/1551503-200.png" width="30" height="30"> Dimensionality reduction | [Seurat_dr](labs/compiled/seurat/seurat_02_dim_reduction.md) | [Scater_dr](labs/compiled/scater/scater_02_dim_reduction.md) | [Scanpy_dr](labs/compiled/scanpy/scanpy_02_dim_reduction.html) |
| <img border="0" src="http://s16574.pcdn.co/wp-content/uploads/2018/05/cluster-icon.png" width="30" height="30"> Data integration | [Seurat_integr](labs/compiled/seurat/seurat_03_integration.md) | [Scater_integr](labs/compiled/scater/scater_03_integration.md) | [Scanpy_integr](labs/compiled/scanpy/scanpy_03_integration.html) |
| <img border="0" src="https://cdn0.iconfinder.com/data/icons/network-analysis-7/64/cluster-positive-negative-group-collection-512.png" width="30" height="30"> Clustering | [Seurat_clust](labs/compiled/seurat/seurat_04_clustering.html) | [Scater_clust](labs/compiled/scater/scater_04_clustering.html) | [Scanpy_clust](labs/compiled/scanpy/scanpy_04_clustering.html) |
| <img border="0" src="https://static.thenounproject.com/png/1517975-200.png" width="30" height="30"> Differential expression | [Seurat_dge](labs/compiled/seurat/seurat_05_dge.html) | [Scater_dge](labs/compiled/scater/scater_05_dge.html) | [Scanpy_dge](labs/compiled/scanpy/scanpy_05_dge.html) |
| <img border="0" src="https://cdn2.vectorstock.com/i/1000x1000/49/51/route-location-icon-vector-16394951.jpg" width="30" height="30"> Trajectory inference | [Slingshot_ti](labs/compiled/slingshot/slingshot.html) | [Slingshot_ti](labs/compiled/slingshot/slingshot.html) | [PAGA_ti](labs/paga/paga.ipynb) |

<br/>

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/759/exercise.svg" width="40" height="40"> BONUS exercises
***

Please note that the exercises listed below belong to past courses and might not be completely updated. Nonetheless, they provide even more details and options to analyse single cell data.

| Name (link) | Description |
| ----------- | ----------- |
| [Read-Pipeline](oldlabs/Pipeline_exercise) | Snakemake pipeline for processing SmartSeq2 data, mapping reads, QC and expression estimates|
| [biomaRt](oldlabs/biomart) | For those not familiar with working with biomaRt, we suggest that you have a look at this example code for how to convert between different formats using biomaRt|
| [PCA, tSNE and clustering](oldlabs/PCA_and_clustering) | Basic PCA, tSNE and clustering using base R on mouse embryonic development data. |
| [KNN-graphs](oldlabs/igraph) | Construction of graphs from cell-cell similiarities using igraph|
| [Estimating Batch-Effects](https://bitbucket.org/scilifelab-lts/scrnaseq-labs/src/a228442debe7f8eff28cfdba875349025db9b7a3/batch_analysis.md?fileviewer=file-view-default) | A tutorial for estimating genome-wide and individual genes batch-effects |
| [Normalization Comparison](oldlabs/norm_analysis_v2) | A tutorial for comparison scRNAseq and bulk RNAseq normalisation strategies. |
| [SC3 package](oldlabs/sc3_R35) | Tutorial with the SC3 consensus clustering package |
| [Trajectory with Monocle2](oldlabs/monocle_analysis) | A tutorial with mouse embryonic data using the Monocle package for pseudotime analysis |
| [Differential expression_2](oldlabs/Differential_gene_expression) | OBS! This old tutorial uses Seurat v2! For this tutorial we have included several different methods for differential expression tests on single cell data, including SCDE, MAST, SC3 and Seurat. The exercise has been split into 2 parts with evaluation of all results in the second part |
| [UPPMAX Sbatch](oldlabs/sbatchScript) | One example of a sbatch script |
| [Pagoda](oldlabs/pagoda_ilc) | Pagoda pathway wPCA for clustering of cells. OBS! several steps in this tutorial takes hours to run if you work with your own dataset, a good suggestion is to start with the first steps, knn.error.model, pagoda.varnorm and pagoda.pathway.wPCA and let it run while working on other tutorials. You can also run it with more than one core to speed things up |


We will try to keep these tutorials up to date. If you find any errors or things that you think should be updated please contact Åsa (asa.bjorklund@scilifelab.se) or Paulo (paulo.czarnewski@scilifelab.se)

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/83019/faq-button.svg" width="40" height="40"> FAQ
***

As you run into problems, we will try to fill in the [FAQ](labs/FAQ) with common questions.

<br/>

<div style="text-align: right; font-size: 5px"> Icons are provided from [www.svgrepo.com](www.svgrepo.com) </div>
