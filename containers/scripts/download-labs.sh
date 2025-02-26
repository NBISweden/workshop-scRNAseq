#! /bin/bash

## Example usage:
#   ./download-labs.sh scanpy
#   ./download-labs.sh seurat

orgurl="https://github.com/NBISweden"
reponame="workshop-scRNAseq-devel"
repodir="compiled/labs"
localdir="labs"
toolkit="$1"

function git_sparse_clone() (
    mkdir -p ${localdir}
    git clone -n --depth=1 --filter=tree:0 ${orgurl}/${reponame} > /dev/null 2>&1
    cd ${reponame}
    git sparse-checkout set --no-cone ${repodir}/${toolkit} > /dev/null 2>&1
    git checkout > /dev/null 2>&1
    cd - > /dev/null 2>&1
    find . -type f -name '*.ipynb' -exec mv -n {} ./${localdir}/ \;
    rm -rf ./${reponame}
)


function select_kernel() (
    notebooks=( $(find . -name "*.ipynb" -print) )
    for nb in "${notebooks[@]}"; do
        jq '.metadata.kernelspec = {"display_name": "scanpy", "language": "python", "name": "scanpy"}' ${nb} > tmp.$$.json && mv tmp.$$.json ${nb}
    done
)


function main() (
    echo "downloading files from ${orgurl}/${reponame}/${repodir}/${toolkit} into ${localdir}/..."
    git_sparse_clone
    echo "making 'scanpy' default kernel..."
    select_kernel
)


main
