# Docker Images

The Docker images are based on [rocker-org/rocker-versioned2](https://github.com/rocker-org/rocker-versioned2) images, and follow the same folder structure.

The R images are based on the `rocker/tidyverse` with a Dockerfile per toolkit and the corresponding install scripts in the `scripts/` directory. The Python image is based on the `jupyter/minimal-notebook`.

> **NOTE:** All images contain `quarto` command line utility.

## Description

Here we use three different toolkits, namely `Seurat` (R/RStudio), `Bioconductor` (R/RStudio) and `Scanpy` (Python/Jupyter) to perform scRNAseq analysis provided in the NBIS scRNAseq workshop. 

The different images are differentiated by the following `registry/username/image:tag` convention:

```
ghcr.io/nbisweden/workshop-scrnaseq:<YEAR>-<TOOLKIT>-<LANGUAGE-VERSION>
```

Each image contains the required packages so that, for each toolkit, the following analysis steps can be performed:
* QC
* Dimensionality reduction
* Data integration
* Clustering
* Differential expression
* Celltype prediction
* Trajectory analysis

## Working with Python from R

For the R based images, a Python environment is created using the `reticulate` R package.

> **NOTE:** It is important to know that depending on the package version and it's dependencies, some packages are installed in a `conda` environment, while others are installed in a `venv` environment. 

For `Bioconductor`, the Python environment includes `scanorama`. To use this environment, add the following lines in your notebook:

```R
library(reticulate)
use_virtualenv("/opt/venv/scanorama")
py_discover_config()
```

Additionally, a `conda` environment is available and can be activated by adding the following lines to your notebook:

```R
library(reticulate)
use_condaenv("bioc", conda = "/opt/conda/bin/conda")
py_discover_config()
```

For `Seurat`, the `conda` environment includes `umap-learn`. To use this environment, add the following lines in your notebook:

```R
library(reticulate)
use_condaenv("seurat", conda = "/opt/conda/bin/conda")
py_discover_config()
``` 

The `scanorama` package is available in the Python environment and can be activated by adding the following lines:

```R
library(reticulate)
use_virtualenv("/opt/venv/scanorama")
py_discover_config()
```

## Installing additional R packages

> **NOTE:** These packages will be installed in the *running container*, and are not persisted when the container is terminated.

If you need to install additional R packages from `Cran`, `Bioconductor` or `GitHub`, the following scripts are available for convenience. From the _Terminal_, run the corresponding command below:

* CRAN  

  ```bash
  install2.r --error --skipinstalled -n "$NCPUS" <PACKAGE-NAME>
  ```

* Bioconductor

  ```bash
  installBioc.r --error --skipinstalled -n "$NCPUS" <PACKAGE-NAME>
  ```

* GitHub

  ```bash
  installGithub.r <GITHUB-OWNER>/<GITHUB-REPOSITORY>
  ```

## Installing additional Python packages

> **NOTE:** These packages will be installed in the *running container*, and are not persisted when the container is terminated.

For the Python based image, you can install packages using `mamba`, `pip`, or `conda` (`mamba` is recommended).

## How To Run

### R Based Images

```bash
docker pull ghcr.io/nbisweden/workshop-scrnaseq:<TAG>
docker run --rm -ti -p 8787:8787 -e PASSWORD=scrnaseq -v /path/to/workdir:/home/rstudio/workdir ghcr.io/nbisweden/workshop-scrnaseq:<TAG>
```

In the browser, go to [localhost:8787](localhost:8787).  
Use the following credentials to log in to the RStudio Server:  
> user: rstudio  
> password: scrnaseq

### Python Based Image

```bash
docker pull ghcr.io/nbisweden/workshop-scrnaseq:<TAG>
docker run --rm -ti -p 8888:8888 -v /path/to/workdir:/home/jovyan/workshop-scRNAseq ghcr.io/nbisweden/workshop-scrnaseq:<TAG>
```

In the browser, go to [localhost:8888](localhost:8888).

## How To Build

From the project root directory, run the following command:

```bash
docker build -t ghcr.io/nbisweden/workshop-scrnaseq:<TOOLKIT> --file containers/dockerfiles/<TOOLKIT>.Dockerfile ./containers
```

## Push to GitHub Container Registry

Each toolkit image comes with an associated GitHub Action that builds and pushes the image to GitHub Container Registry (`ghcr.io`). Each workflow is triggered only on changes to files related to corresponding toolkit image.

## Build Singularity Images (Uppmax)

To use these images on Uppmax, you first need to build them as Singularity images. A set of minimal Singularity definition files (`.def`) are provided in the `uppmax/` folder. A set of launch scripts, one for JupyterLab based images and one for RStudio based images are also provided in the same folder.

To build singularity images, use this bash script;

```
#!/bin/bash
## Build singularity sif containers

#SBATCH -A naiss2023-22-1345
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:30:00
#SBATCH -J bld-sif

# start date time
starttime=`date +%s`

export APPTAINER_CACHEDIR=${PWD}
export SINGULARITY_CACHEDIR=${PWD}

singularity cache clean --all

singularity build --force 2024-seurat-r4.3.0.sif singularity_seurat.def
singularity build --force 2024-bioconductor-r4.3.0.sif singularity_bioconductor.def
singularity build --force 2024-scanpy-py3.10.sif singularity_scanpy.def

singularity build --force 2024-seurat_spatial-r4.3.0.sif singularity_seurat_spatial.def
singularity build --force 2024-bioconductor_spatial-r4.3.0.sif singularity_bioconductor_spatial.def
singularity build --force 2024-scanpy_spatial-py3.10.sif singularity_scanpy_spatial.def

singularity build --force 2024-site-r4.3.0.sif singularity_site.def

#end date time
endtime=`date +%s`
echo "End of Script. Script took $(($endtime-$starttime)) seconds."
exit 0
```

## Build site

To render all qmds files and thereby the website as well as generate compiled files, clone the repo, then run this in the root of the repo:

```
## assumes singularity images are at /sw/courses/scrnaseq/singularity/
sbatch scripts/render-singularity.sh
```
