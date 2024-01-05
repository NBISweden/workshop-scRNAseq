#!/bin/sh

## RENDER QMD TO HTML
##
## Description
## Runs qmd files and generates html files into the docs/ directory
## 
## Usage
## Run this script in the root of the repo
## bash ./scripts/render.sh

## fail fast
set -e

## define docker images
docker_seurat="ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat-r4.3.0"
docker_bioc="ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor-r4.3.0"
docker_scanpy="ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy-py3.10"

docker_seurat_spatial="ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat_spatial-r4.3.0"
docker_bioc_spatial="ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor_spatial-r4.3.0"
docker_scanpy_spatial="ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy_spatial-py3.10"

docker_site="ghcr.io/nbisweden/workshop-scrnaseq:2024-site-r4.3.0"

# check if in the root of the repo
if [ ! -f "_quarto.yml" ]; then
    echo "Error: Are you in the root of the repo? _quarto.yml is missing."
    exit 1
fi

# start time for whole script
start=$(date +%s.%N)

## seurat
echo "Rendering Seurat files..."
start_seurat=$(date +%s.%N)
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_01_qc.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_02_dimred.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_03_integration.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_04_clustering.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_05_dge.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_06_celltyping.qmd
# docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat_spatial quarto render /work/labs/seurat/seurat_07_spatial.qmd
# docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_seurat quarto render /work/labs/seurat/seurat_08_trajectory.qmd
duration_seurat=$(echo "$(date +%s.%N) - $start_seurat" | bc) && echo "Seurat time elapsed: $duration_seurat seconds"

## bioconductor
echo "Rendering Bioconductor files..."
start_bioc=$(date +%s.%N)
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc quarto render /work/labs/bioc/bioc_01_qc.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc quarto render /work/labs/bioc/bioc_02_dimred.qmd
# docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc quarto render /work/labs/bioc/bioc_03_integration.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc quarto render /work/labs/bioc/bioc_04_clustering.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc quarto render /work/labs/bioc/bioc_05_dge.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc quarto render /work/labs/bioc/bioc_06_celltyping.qmd
docker run --rm --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work $docker_bioc_spatial quarto render /work/labs/bioc/bioc_07_spatial.qmd
duration_bioc=$(echo "$(date +%s.%N) - $start_bioc" | bc) && echo "Bioc time elapsed: $duration_bioc seconds"

## scanpy
echo "Rendering Scanpy files..."
start_scanpy=$(date +%s.%N)
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_01_qc.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_02_dimred.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_03_integration.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_04_clustering.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_05_dge.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_06_celltyping.qmd
#docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy_spatial quarto render /work/labs/scanpy/scanpy_07_spatial.qmd
docker run --rm --platform=linux/amd64 -p 8888:8888 -v $PWD:/work $docker_scanpy quarto render /work/labs/scanpy/scanpy_08_trajectory.qmd
duration_scanpy=$(echo "$(date +%s.%N) - $start_scanpy" | bc) && echo "Scanpy time elapsed: $duration_scanpy seconds"

## site
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/index.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/home_contents.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/home_info.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/home_precourse.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/home_schedule.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/home_syllabus.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/other/uppmax.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/other/docker.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/other/containers.qmd
docker run --rm --platform=linux/amd64 -v $PWD:/work $docker_site quarto render /work/other/faq.qmd

echo "Seurat time elapsed: $duration_seurat seconds"
echo "Bioc time elapsed: $duration_bioc seconds"
echo "Scanpy time elapsed: $duration_scanpy seconds"
duration=$(echo "$(date +%s.%N) - $start" | bc) && echo "Total time elapsed: $duration seconds"

echo "All files rendered successfully."
exit 0
