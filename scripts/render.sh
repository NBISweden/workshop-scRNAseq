#!/bin/bash

## RENDER QMD TO HTML
##
## Description
## Runs qmd files and generates html files into the docs/ directory
## 
## Usage
## Run this script in the root of the repo
## bash ./scripts/render.sh [all|seurat|bioc|scanpy|spatial|site|compile]
## Optionally: DOCKER_R=your/image:tag bash ./scripts/render.sh [option]

set -euo pipefail

# Docker image variables (can be overridden via env)
DOCKER_R="${DOCKER_R:-ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311}"
DOCKER_SCANPY="${DOCKER_SCANPY:-ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256}"
DOCKER_SEURAT_SPATIAL="${DOCKER_SEURAT_SPATIAL:-ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat_spatial-r4.3.0}"
DOCKER_BIOC_SPATIAL="${DOCKER_BIOC_SPATIAL:-ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor_spatial-r4.3.0}"
DOCKER_SITE="${DOCKER_SITE:-ghcr.io/nbisweden/workshop-scrnaseq:2024-site-r4.3.0}"

# Entry points for conda envs
ENTRYPOINT_R="/usr/local/conda/bin/conda"
ENTRYPOINT_SCANPY="/opt/conda/bin/conda"

LAB_DIR="docs/labs"
LECTURE_DIR="docs/lectures"
SITE_DIR="docs"
OTHER_DIR="docs/other"

# Argument parsing
TOOLKIT="${1:-}"

usage() {
    echo "Usage: $0 [seurat|bioc|scanpy|spatial|site|compile|all]"
    echo "Optionally: DOCKER_<TOOLKIT>=your/image:tag $0 [option]"
    exit 1
}

# Timer function
timer_start() {
    date +%s.%N
}

timer_report() {
    local start=$1
    local end
    end=$(date +%s.%N)
    duration=$(echo "$end - $start" | bc)
    echo "Time elapsed: $duration seconds"
}

# Render a list of files using docker/conda
render_files_r() {
    local files=("$@")
    local start
    start=$(timer_start)
    for file in "${files[@]}"; do
        echo "Rendering $file ..."
        docker run --rm -it --platform=linux/amd64 -u root -v "${PWD}:/home/jovyan/work" \
            --entrypoint "$ENTRYPOINT_R" "$DOCKER_R" run -n seurat quarto render "/home/jovyan/work/$file"
    done
    timer_report "$start"
}

render_files_scanpy() {
    local files=("$@")
    local start
    start=$(timer_start)
    for file in "${files[@]}"; do
        echo "Rendering $file ..."
        docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" \
            --entrypoint "$ENTRYPOINT_SCANPY" "$DOCKER_SCANPY" run -n scanpy quarto render "/work/$file"
    done
    timer_report "$start"
}

render_files_site() {
    local files=("$@")
    local start
    start=$(timer_start)
    for file in "${files[@]}"; do
        echo "Rendering $file ..."
        docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" \
            "$DOCKER_SITE" quarto render "/work/$file"
    done
    timer_report "$start"
}

render_files_spatial() {
    echo ""
    echo "NOTICE: Spatial labs are no longer included in the workshop!"
    echo ""
    # local seurat_file="$LAB_DIR/seurat/seurat_08_spatial.qmd"
    # local bioc_file="$LAB_DIR/bioc/bioc_08_spatial.qmd"
    # local scanpy_file="$LAB_DIR/scanpy/scanpy_08_spatial.qmd"
    # local start
    # start=$(timer_start)
    # echo "Rendering $seurat_file ..."
    # docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" \
    #     "$DOCKER_SEURAT_SPATIAL" quarto render "/work/$seurat_file"
    # echo "Rendering $bioc_file ..."
    # docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" \
    #     "$DOCKER_BIOC_SPATIAL" quarto render "/work/$bioc_file"
    # echo "Rendering $scanpy_file ..."
    # docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" \
    #     --entrypoint "$ENTRYPOINT_SCANPY" "$DOCKER_SCANPY" run -n scanpy quarto render "/work/$scanpy_file"
    # timer_report "$start"
}

render_files_lectures() {
    local lecture_files=("$LECTURE_DIR/gsa/index.qmd")
    local start
    start=$(timer_start)
    for file in "${lecture_files[@]}"; do
        echo "Rendering $file ..."
        docker run --rm -it --platform=linux/amd64 -v "${PWD}:/home/jovyan/work" \
            --entrypoint "$ENTRYPOINT_R" "$DOCKER_R" run -n seurat quarto render "/home/jovyan/work/$file"
    done
    timer_report "$start"
}

render_files_site_pages() {
    local site_files=(
        "$LAB_DIR/index.qmd"
        "$OTHER_DIR/containers-spatial.qmd"
        "$OTHER_DIR/containers.qmd"
        "$OTHER_DIR/data.qmd"
        "$OTHER_DIR/docker.qmd"
        "$OTHER_DIR/faq.qmd"
        "$OTHER_DIR/scilifelab-serve.qmd"
        "$OTHER_DIR/uppmax.qmd"
        "$SITE_DIR/404.qmd"
        "$SITE_DIR/home_contents.qmd"
        "$SITE_DIR/home_info.qmd"
        "$SITE_DIR/home_precourse.qmd"
        "$SITE_DIR/home_schedule.qmd"
        "$SITE_DIR/home_syllabus.qmd"
        "$SITE_DIR/index.qmd"
    )
    render_files_site "${site_files[@]}"
}

render_seurat() {
    local files=(
        "$LAB_DIR/seurat/seurat_01_qc.qmd"
        "$LAB_DIR/seurat/seurat_02_dimred.qmd"
        "$LAB_DIR/seurat/seurat_03_integration.qmd"
        "$LAB_DIR/seurat/seurat_04_clustering.qmd"
        # "$LAB_DIR/seurat/seurat_05_dge.qmd"
        "$LAB_DIR/seurat/seurat_06_celltyping.qmd"
        "$LAB_DIR/seurat/seurat_07_trajectory.qmd"
    )
    echo "Rendering Seurat files..."
    render_files_r "${files[@]}"
}

render_bioc() {
    local files=(
        "$LAB_DIR/bioc/bioc_01_qc.qmd"
        "$LAB_DIR/bioc/bioc_02_dimred.qmd"
        "$LAB_DIR/bioc/bioc_03_integration.qmd"
        "$LAB_DIR/bioc/bioc_04_clustering.qmd"
        # "$LAB_DIR/bioc/bioc_05_dge.qmd"
        "$LAB_DIR/bioc/bioc_06_celltyping.qmd"
    )
    echo "Rendering Bioconductor files..."
    render_files_r "${files[@]}"
}

render_scanpy() {
    local files=(
        "$LAB_DIR/scanpy/scanpy_01_qc.qmd"
        "$LAB_DIR/scanpy/scanpy_02_dimred.qmd"
        "$LAB_DIR/scanpy/scanpy_03_integration.qmd"
        "$LAB_DIR/scanpy/scanpy_04_clustering.qmd"
        # "$LAB_DIR/scanpy/scanpy_05_dge.qmd"
        "$LAB_DIR/scanpy/scanpy_06_celltyping.qmd"
        "$LAB_DIR/scanpy/scanpy_07_trajectory.qmd"
    )
    echo "Rendering Scanpy files..."
    render_files_scanpy "${files[@]}"
}

compile_all() {
    if [ ! -f "./scripts/compile.sh" ]; then
        echo "Error: Cannot find ./scripts/compile.sh for compile step."
        exit 1
    fi
    echo "Compiling labs ..."
    bash ./scripts/compile.sh "all"
    echo "All labs compiled successfully."
}

main() {
    local start
    start=$(timer_start)

    if [ -z "${TOOLKIT}" ]; then
        usage
    fi

    case "${TOOLKIT}" in
        seurat)
            render_seurat
            ;;
        bioc)
            render_bioc
            ;;
        scanpy)
            render_scanpy
            ;;
        spatial)
            render_files_spatial
            ;;
        site)
            render_files_lectures
            render_files_site_pages
            ;;
        compile)
            compile_all
            ;;
        all)
            compile_all
            render_seurat
            render_bioc
            render_scanpy
            render_files_spatial
            render_files_lectures
            render_files_site_pages
            ;;
        *)
            echo "Unknown option '${TOOLKIT}'."
            usage
            ;;
    esac

    timer_report "$start"
    echo "All files rendered successfully."
}

main
