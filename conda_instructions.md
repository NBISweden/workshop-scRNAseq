---
layout: default
title:  'Schedule - scRNAseq course'
---

<br/>

#### <img border="0" src="https://hackernoon.com/hn-images/1*rW03Wtue71AKfxnx6XN_iQ.png" width="40" height="40"> Conda Instructions
***

In this workshop you will use the classroom computers, so you’ll have all the necessary software installed and ready to run everything.

Briefly, conda allows users to run pieces of code just like the developer created them.


<img border="0" src="https://nbisweden.github.io/excelerate-scRNAseq/logos/conda_illustration.png" width="400">


If you want to install the same software on your own laptop after the course, you can follow the
instructions below to do so. These rely on [Conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html), a self-contained directory that
you can use in order to reproduce all your results. Two of the required software are not available as Conda packages, please see the separate [instructions for installing SingleR and CHETAH](https://raw.githubusercontent.com/NBISweden/excelerate-scRNAseq/master/notes_installation.txt).

Briefly, you need to:  

1. Install Conda and download the `.yaml` file
2. Create and activate the environment
3. Deactivate the environment after running your analyses

You can [read more](https://nbis-reproducible-research.readthedocs.io/en/latest/conda/) about Conda environments and other important concepts to help you make your research reproducible.

<br/>

##### Install Conda and download the environment file
***

You should start by installing Conda. We suggest installing either Miniconda, or Anaconda if storage is
not an issue. After [installing Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html),
download the course Conda file and put it in your working folder. Note that there is a separate [environment.yaml for Mac](https://raw.githubusercontent.com/NBISweden/excelerate-scRNAseq/master/environment.yaml) and for [environment_centos7.yaml for Linux](https://raw.githubusercontent.com/NBISweden/excelerate-scRNAseq/master/environment_centos7.yaml).

<br/>

##### Create an environment from file
***

In terminal, `cd` to your working folder and create an environment called `scRNAseq2020` from the
`environment.yaml` file.

```
conda env create -p scRNAseq2020 -f environment.yaml
```

Several messages will show up on your screen and will tell you about the installation process. This may take a few minutes depending on how many packages are to be installed.

```
##Collecting package metadata: done
##Solving environment: done
##
##Downloading and Extracting Packages
##libcblas-3.8.0       | 6 KB      | ############################################################################# | 100%
##liblapack-3.8.0      | 6 KB      | ############################################################################# | 100%
##...
##Preparing transaction: done
##Verifying transaction: done
##Executing transaction: done
```

<br/>

##### Activate the environment
***

To activate an environment type:

```
source activate ./scRNAseq2020
```

From this point on you can run any of the contents from the course. For instance, you can directly launch RStudio by
typing `rstudio`. Here it is important to add the `&` symbol in the end to be able to use the command line at the same time if needed. You can open other files from Rstudio as well.

```
rstudio ./labs/compiled/my_script.Rmd &
```

<br/>

##### Deactivate the environment
***

After you've ran all your analyses, deactivate the environment by typing:

```
conda deactivate
```
