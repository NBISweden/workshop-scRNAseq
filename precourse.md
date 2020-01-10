---
layout: default
title:  'Precourse Material - scRNAseq course'
---

#### <img border="0" src="https://www.svgrepo.com/show/19652/maths-class-materials-cross-of-a-pencil-and-a-ruler.svg" width="40" height="40"> Precourse material
***

##### <img border="0" src="https://www.svgrepo.com/show/4795/installation-symbol.svg" width="20" height="20"> Installations

We have conda recipies for all R and Python packages in one [file](labs/environment_r.yml). If you have never worked with conda before, please read the [conda instructions](conda_instructions.md).

OBS! Need to fix some paths in instruction.
Also info on Docker?

<br/>

##### <img border="0" src="https://www.svgrepo.com/show/20109/database.svg" width="20" height="20"> Dataset

We will run all tutorials with a set of 3 PBMC 10x datasets from the 10X Genomics website, with different types of library preps. One dataset was done with 10x version2, one with 10x version3 and one with 10x version 3 combined with CITE-seq protein detection. 

These can be fetched using commands (which are also included in the initial labs):

```
mkdir data  
cd data
curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v2/pbmc_1k_v2_filtered_feature_bc_matrix.h5
curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5
curl -O http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5
```

##### <img border="0" src="https://www.svgrepo.com/show/20109/database.svg" width="20" height="20"> Code

All code is also available as R-markdown scripts at the course Github space. If you want a local copy of all course material you can simply clone it with git using:

````
git clone https://github.com/NBISweden/workshop-scRNAseq.git
```

Or download manually from the github site https://github.com/NBISweden/workshop-scRNAseq. 

<br/>


##### <img border="0" src="https://www.svgrepo.com/show/17086/server-client-exchange.svg" width="20" height="20"> Uppmax

**Attention**: This step is no longer required for the course. It is only used for the optional pipeline exercise.


1.   If you do not already have an Uppmax account, create and Uppmax account following these [instructions](files/Apply_for_Uppmax_account.pdf). OBS! It may take a few days to recieve the account, so proceed with this point as soon as possible.

2.   Log in to SUPR and request membership in the project g2019002. Account approval requires manual confirmation from the course organizers, so it may not happen immediately.

3.   Make sure you can connect to Rackham at Uppmax using a terminal. If you use a pc we recommend MobaXterm (http://mobaxterm.mobatek.net) or Windows 10 Bash for Linux.

4.   If you still feel uncertain how to work in a terminal please take time to do the first three parts in the “Unix tutorial for beginners” that can be found here http://www.ee.surrey.ac.uk/Teaching/Unix/ before the course starts. Otherwise you will not be able to take in the practical parts.  

5.   Make sure that you can read and write in the course folder by creating a file with your uppmax user name in the `/proj/g2019002/completed` folder. If you cannot write to the folder, the most likely reason is that you have not requested access to the course project via [SUPR](https://supr.snic.se/), see point 2. OBS! It may take an hour or so from the time your access is approved before you can actually write to the folder. We will check before the course that all students have logged in and done this, so do not forget!


