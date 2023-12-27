#!/bin/bash

## RENDER QMD TO HTML
## Runs qmd files and generates html files into the docs/ directory
## Run this script in the root of the repo

## fail fast
set -e

## define docker images
docker_seurat="ghcr.io/nbisweden/workshop-scrnaseq:2023-seurat-r4.3.0"
docker_bioc="ghcr.io/nbisweden/workshop-scrnaseq:2023-bioconductor-r4.3.0"
docker_scanpy="ghcr.io/nbisweden/workshop-scrnaseq:2023-scanpy-py3.10"

# check if in the root of the repo
if [ ! -f "_quarto.yml" ]; then
    echo "Error: Are you in the root of the repo? _quarto.yml is missing."
    exit 1
fi

## site
# docker run --rm --platform=linux/amd64 -v $PWD:/home/rstudio/workdir $docker_image quarto render /home/rstudio/workdir/index.qmd
# docker run --rm --platform=linux/amd64 -v $PWD:/home/rstudio/workdir $docker_image quarto render /home/rstudio/workdir/home_contents.qmd
# docker run --rm --platform=linux/amd64 -v $PWD:/home/rstudio/workdir $docker_image quarto render /home/rstudio/workdir/home_info.qmd
# docker run --rm --platform=linux/amd64 -v $PWD:/home/rstudio/workdir $docker_image quarto render /home/rstudio/workdir/home_precourse.qmd
# docker run --rm --platform=linux/amd64 -v $PWD:/home/rstudio/workdir $docker_image quarto render /home/rstudio/workdir/home_schedule.qmd
# docker run --rm --platform=linux/amd64 -v $PWD:/home/rstudio/workdir $docker_image quarto render /home/rstudio/workdir/home_syllabus.qmd

## seurat
echo "Rendering Seurat files..."
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_01_qc.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_02_dimred.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_03_integration.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_04_clustering.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_05_dge.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_06_celltyping.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_07_spatial.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_seurat quarto render /home/rstudio/workdir/labs/seurat/seurat_08_trajectory.qmd

## bioconductor
echo "Rendering Bioconductor files..."
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_01_qc.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_02_dimred.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_03_integration.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_04_clustering.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_05_dge.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_06_celltyping.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir $docker_bioc quarto render /home/rstudio/workdir/labs/bioc/bioc_07_spatial.qmd

## scanpy
echo "Rendering Scanpy files..."
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_01_qc.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_02_dimred.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_03_integration.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_04_clustering.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_05_dge.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_06_celltyping.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_07_spatial.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir $docker_scanpy quarto render /home/jovyan/workdir/labs/scanpy/scanpy_08_trajectory.qmd

echo "All files rendered successfully."
exit 0
