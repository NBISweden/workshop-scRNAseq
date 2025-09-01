#!/bin/sh

## RENDER QMD TO HTML
##
## Description
## Runs qmd files and generates html files into the docs/ directory
## 
## Usage
## Run this script in the root of the repo
## bash ./scripts/render.sh option

## input argument can  be either of:
## all - run all the steps
## seurat - render all seurat labs
## bioc - render all bioc labs
## scanpy - render all scanpy labs
## spatial - render all 3 spatial labs.
## site - render all site stuff.
## compile - compile labs into Rmd/ipynb


## fail fast
set -e

## define docker images
docker_r="ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311"
docker_scanpy="ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256"

## old images for the r spatial labs.
docker_seurat_spatial="ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat_spatial-r4.3.0"
docker_bioc_spatial="ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor_spatial-r4.3.0"

# old site image still works
docker_site="ghcr.io/nbisweden/workshop-scrnaseq:2024-site-r4.3.0"

# check if in the root of the repo
if [ ! -f "_quarto.yml" ]; then
    echo "Error: Are you in the root of the repo? _quarto.yml is missing."
    exit 1
fi

# start time for whole script
start=$(date +%s.%N)

# -u 1000:1000 is useful on linux

## seurat
## OBS! now running the containers with the conda env as entrypoint, then run -n seurat = refers to the conda environment named "seurat"
if [[ "$@" =~ 'seurat' ]] ||  [[ "$@" =~ 'all' ]]
then
    echo "Rendering Seurat files..."
    start_seurat=$(date +%s.%N)
    docker run --rm -it -u root --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_01_qc.qmd
    docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_02_dimred.qmd
    docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_03_integration.qmd
    docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_04_clustering.qmd
    docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_05_dge.qmd
    docker run --rm -it -u root --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_06_celltyping.qmd
    docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/seurat/seurat_07_trajectory.qmd
duration_seurat=$(echo "$(date +%s.%N) - $start_seurat" | bc) && echo "Seurat time elapsed: $duration_seurat seconds"
echo "time elapsed: $duration_seurat seconds"
fi 

if [[ "$@" =~ 'bioc' ]]  ||  [[ "$@" =~ 'all' ]]
then
	## bioconductor, same conda env as seurat.
	echo "Rendering Bioconductor files..."
	start_bioc=$(date +%s.%N)
	docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/bioc/bioc_01_qc.qmd
	docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/bioc/bioc_02_dimred.qmd
	docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/bioc/bioc_03_integration.qmd
	docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/bioc/bioc_04_clustering.qmd
	docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/bioc/bioc_05_dge.qmd
	docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/labs/bioc/bioc_06_celltyping.qmd
	duration_bioc=$(echo "$(date +%s.%N) - $start_bioc" | bc) && echo "Bioc time elapsed: $duration_bioc seconds"
	echo "Bioc time elapsed: $duration_bioc seconds"
fi

## scanpy
if [[ "$@" =~ 'scanpy' ]]  ||  [[ "$@" =~ 'all' ]]
then    
    echo "Rendering Scanpy files..."
    start_scanpy=$(date +%s.%N)
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_01_qc.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_02_dimred.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_03_integration.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_04_clustering.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_05_dge.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_06_celltyping.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_07_trajectory.qmd
    duration_scanpy=$(echo "$(date +%s.%N) - $start_scanpy" | bc) && echo "Scanpy time elapsed: $duration_scanpy seconds"
    echo "Scanpy time elapsed: $duration_scanpy seconds"
fi


if [[ "$@" =~ 'spatial' ]]  ||  [[ "$@" =~ 'all' ]]
then
    echo "Rendering Spatial files..."
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_seurat_spatial quarto render /work/labs/seurat/seurat_08_spatial.qmd
    duration_seurat=$(echo "$(date +%s.%N) - $start_seurat" | bc) && echo "Seurat time elapsed: $duration_seurat seconds"
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_bioc_spatial quarto render /work/labs/bioc/bioc_08_spatial.qmd    
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint "/opt/conda/bin/conda" $docker_scanpy  run -n scanpy quarto render /work/labs/scanpy/scanpy_08_spatial.qmd
fi


if [[ "$@" =~ 'site' ]]  ||  [[ "$@" =~ 'all' ]]
then
    echo "Rendering lectures.."

    
    ## lectures, only 2 that are created with qmd.

    # dge requires ggpubr so this is currently being rendered interactively. ggpubr should be added to the container for next year
    # docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/lectures/dge/index.qmd
    docker run --rm -it --platform=linux/amd64 -v ${PWD}:/home/jovyan/work --entrypoint "/usr/local/conda/bin/conda" $docker_r run -n seurat quarto render /home/jovyan/work/lectures/gsa/index.qmd    

    ## site
    echo "Rendering site files..."
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/index.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/home_contents.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/home_info.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/home_precourse.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/home_schedule.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/home_syllabus.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/uppmax.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/scilifelab-serve.qmd    
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/docker.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/containers.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/containers-spatial.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/faq.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/other/data.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/404.qmd
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work $docker_site quarto render /work/labs/index.qmd

    echo "All of site rendered."
fi




# build compiled files
if [[ "$@" =~ 'compile' ]]  ||  [[ "$@" =~ 'all' ]]
then
    bash ./scripts/compile.sh
    echo "All labs compiled successfully."
fi    


duration=$(echo "$(date +%s.%N) - $start" | bc) && echo "Total time elapsed: $duration seconds"

echo "All files rendered successfully."
exit 0
