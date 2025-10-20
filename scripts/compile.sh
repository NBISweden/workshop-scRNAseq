#!/bin/bash

# fail on error and unset variables, encourage good scripting
set -euo pipefail

## COMPILE QMD FILES
##
## Description
## Evaluates and replaces meta variable shortcuts
## .qmd files are converted to .qmd files
## .qmd files with jupyter engine are converted to .ipynb
## The following directories are expected:
## docs/labs/seurat, docs/labs/bioc, docs/labs/scanpy
##
## Usage
## Run this script in the root of repo. It takes about 1 min to run.
## bash ./scripts/compile.sh [all|seurat|bioc|scanpy]
## Optionally: DOCKER_SITE=your/image:tag bash ./scripts/compile.sh [option]

DOCKER_SITE="${DOCKER_SITE:-ghcr.io/nbisweden/workshop-scrnaseq:2024-site-r4.3.0}"
OUTPUT_DIR="compiled"
LAB_DIR="docs/labs"
TOOLKIT="${1:-}"

usage() {
    echo "Usage: $0 [seurat|bioc|scanpy|all]"
    echo "Optionally: DOCKER_SITE=your/image:tag $0 [option]"
    exit 1
}

# check if input directory exists
check_input_dir() {
    if [ ! -d "$LAB_DIR/$1" ]; then
        echo "Error: Directory $LAB_DIR/$1 does not exist."
        exit 1
    fi
}

# remove output subdirectory if it exists
check_output_dir() {
    if [ -d "${OUTPUT_DIR}/labs/$1" ]; then
        echo "Directory ${OUTPUT_DIR}/labs/$1 exists. Removing it"
        rm -r "${OUTPUT_DIR}/labs/$1"
    fi
}

# create compiled versions of qmd using profile "compiled"
quarto_compile() {
    echo "Compiling $1 labs ..."
    docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" "$DOCKER_SITE" quarto render --profile compile "/work/$LAB_DIR/$1"/*.qmd --to markdown-header_attributes --metadata engine:markdown --log-level warning,error
}

# Read an md/qmd, remove unnecessary lines from yaml frontmatter
_slim_md_frontmatter() {
    local temp_file
    temp_file=$(mktemp)
    awk '
        /^---$/ && !in_yaml {in_yaml=1; print; next} 
        /^---$/ && in_yaml {in_yaml=0; print; next} 
        !in_yaml {print} 
        in_yaml {
            if (/^(title:|subtitle:|description:)/) {
                print;
                continue_capture=1;
            } else if (continue_capture && !(/^[a-zA-Z0-9_-]+:/)) {
                print;
            } else {
                continue_capture=0;
            }
        }
    ' "$1" > "$temp_file"
    mv "$temp_file" "$1"
}

slim_md_frontmatter() {
    echo "Slimming frontmatter yaml across all $1 .md files"
    find "${OUTPUT_DIR}/labs/$1" -type f -name "*.md" -print0 | while IFS= read -r -d '' file; do
        _slim_md_frontmatter "$file"
    done
}

md_to_qmd() {
    echo "Converting $1 .md files to .qmd"
    shopt -s nullglob
    for file in "${OUTPUT_DIR}/labs/$1"/*.md; do
        mv "$file" "${file%.md}.qmd"
    done
    shopt -u nullglob
}

qmd_to_ipynb() {
    echo "Converting $1 .qmd files to .ipynb"
    shopt -s nullglob
    for file in "${OUTPUT_DIR}/labs/$1"/*.qmd; do
        docker run --rm --platform=linux/amd64 -u 1000:1000 -v "${PWD}:/work" "$DOCKER_SITE" quarto convert "/work/${file}"
        rm -f "$file"
    done
    shopt -u nullglob
}

# Compile one lab
compile_labs() {
    check_input_dir "$1"
    check_output_dir "$1"
    quarto_compile "$1"
    slim_md_frontmatter "$1"
    md_to_qmd "$1"
    if [ "$1" = "scanpy" ]; then
        qmd_to_ipynb "$1"
    fi
}

# main logic, argument parsing
main() {
    if [ -z "${TOOLKIT}" ]; then
        usage
    fi

    case "${TOOLKIT}" in
        seurat|bioc|scanpy)
            compile_labs "${TOOLKIT}"
            ;;
        all)
            for lab in seurat bioc scanpy; do
                compile_labs "$lab"
            done
            ;;
        *)
            echo "Unknown option '${TOOLKIT}'."
            usage
            ;;
    esac
    echo "All files compiled successfully."
}

main
