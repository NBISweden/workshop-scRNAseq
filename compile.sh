#!/bin/bash

## COMPILE QMD FILES
##
## Description
## Evaluates and replaces meta variable shortcuts
## .qmd files are converted to .md files
## .md files with jupyter engine is converted to .ipynb
## The following directories are exprected:
## labs/seurat, labs/bioc, labs/scanpy
##
## Usage
## run this script in the root of repo
## bash ./compile.sh
## takes about 1 min to run

# set output directory for compiled output
output_dir="compiled"
# set input directories with qmd files to compile
input_dirs=("labs/seurat" "labs/bioc" "labs/scanpy")

# fail on error
set -e

# check if in the root of the repo
if [ ! -f "_quarto.yml" ]; then
    echo "Error: Are you in the root of the repo? _quarto.yml is missing."
    exit 1
fi

# check if these directories exist
error=false
for dir in "${input_dirs[@]}"; do
    if [ ! -d "$dir" ]; then
        echo "Error: Directory '$dir' does not exist."
        error=true
    fi
done

if [ "$error" = true ]; then
    exit 1  # Exit with an error code
fi

# if output directory exists, remove it
if [ -d "$output_dir" ]; then
    echo "Directory '$output_dir' exists. Removing it ..."
    rm -r "$output_dir"
fi

# create compiled versions of qmd to using profile "compiled"
echo "Compiling .qmd to .md ..."

quarto render --profile compile labs/seurat/*.qmd --to markdown-header_attributes --metadata engine:markdown --log-level warning,error
quarto render --profile compile labs/bioc/*.qmd --to markdown-header_attributes --metadata engine:markdown --log-level warning,error
quarto render --profile compile labs/scanpy/*.qmd --to markdown-header_attributes --metadata engine:markdown --log-level warning,error


# Read an md/qmd, remove unnecessary lines from yaml, and write to the original file
echo "Slimming yaml across all .md files ..."
slim_yaml() {
    awk '/^---$/ && !in_yaml {in_yaml=1; print; next} 
         /^---$/ && in_yaml {in_yaml=0; print; next} 
         !in_yaml {print} 
         in_yaml {if (/^(title:|subtitle:|description:)/) {print}}' "$1" > "tmp"

    mv "tmp" "$1"
}

find "${output_dir}" -type f -name "*.md" -print0 | while IFS= read -r -d '' file; do
    ## echo "Slimming yaml: $file ..."
    slim_yaml "$file"
done

# converting scanpy md files to ipynb
echo "Converting scanpy .md files to .ipynb ..."
for file in "${output_dir}"/labs/scanpy/*.md; do
    quarto convert "$file" --output "${file%.md}.ipynb"
done

echo "Completed successfully."
