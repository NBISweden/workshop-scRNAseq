#! /bin/bash

## Example usage:
#   ./download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq-devel" "compiled/labs" "scanpy" "labs"

orgurl="$1"
reponame="$2"
repodir="$3"
toolkit="$4"
localdir="$5"

function select_kernel() (
    for nb in $(find $1 -name "*.ipynb"); do
        jq '.metadata.kernelspec = {"display_name": "scanpy", "language": "python", "name": "scanpy"}' ${nb} > tmp.$$.json && mv tmp.$$.json ${nb}
    done
)

function main() (
    echo "downloading files from ${orgurl}/${reponame}/${repodir}/${toolkit} into ${localdir}/..."
    mkdir -p ${localdir}
    git clone -n --depth=1 --filter=tree:0 ${orgurl}/${reponame} > /dev/null 2>&1
    cd ${reponame}
    git sparse-checkout set --no-cone ${repodir}/${toolkit} > /dev/null 2>&1
    git checkout > /dev/null 2>&1
    cd - > /dev/null 2>&1
    if [[ "${toolkit}" == "scanpy" ]]
    then
        echo "processing .ipynb"
        find ${reponame} -type f -name "*.ipynb" -exec mv -n {} ./${localdir}/ \;
        echo "making 'scanpy' default kernel..."
        select_kernel ${localdir}
    else
        echo "processing .qmd"
        find ${reponame} -type f -name "*.qmd" -exec mv -n {} ./${localdir}/ \;
    fi
    rm -rf ./${reponame}
)

main
