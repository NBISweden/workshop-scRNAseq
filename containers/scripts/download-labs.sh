#! /bin/bash

## Example usage:
#   ./download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "seurat" "work/labs"
#   ./download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "scanpy" "work/labs"

ORG_URL="$1"
REPO_NAME="$2"
REPO_DIR="$3"
TOOLKIT="$4"
LOCAL_DIR="$5"

HOME="/home/jovyan"
WORK="${HOME}/work"
KERNEL="scanpy"

## Jupyter functions
function make_kernel() (
    mkdir -p ${HOME}/.local/share/jupyter/kernels/${KERNEL}
    cat <<EOF > ${HOME}/.local/share/jupyter/kernels/${KERNEL}/kernel.json
{
  "argv": [
    "pixi", "run", "--frozen", "--manifest-path", "${WORK}/pixi.toml", "python",
    "-Xfrozen_modules=off",
    "-m",
    "ipykernel_launcher",
    "-f",
    "{connection_file}"
  ],
  "display_name": "${KERNEL}",
  "language": "python",
  "metadata": {
    "debugger": true
  },
  "kernel_protocol_version": "5.5"
}
EOF
)

function select_kernel() (
    for nb in $(find $1 -name "*.ipynb"); do
        jq --arg k "$KERNEL" '.metadata.kernelspec = {"display_name": $k, "language": "python", "name": $k}' ${nb} > tmp.$$.json && mv tmp.$$.json ${nb}
    done
)

function main() (
    echo "downloading files from ${ORG_URL}/${REPO_NAME}/${REPO_DIR}/${TOOLKIT} into ${LOCAL_DIR}/..."
    mkdir -p ${LOCAL_DIR}
    git clone -n --depth=1 --filter=tree:0 ${ORG_URL}/${REPO_NAME} > /dev/null 2>&1
    cd ${REPO_NAME}
    git sparse-checkout set --no-cone "${REPO_DIR}/${TOOLKIT}" "${REPO_DIR}/figs" > /dev/null 2>&1
    git checkout > /dev/null 2>&1
    cd - > /dev/null 2>&1
    
    if [ -d "${REPO_NAME}/${REPO_DIR}/figs" ]; then
        echo "moving figs directory..."
        cp -r "${REPO_NAME}/${REPO_DIR}/figs" "./"
    fi
    
    
    if [[ "${TOOLKIT}" == "scanpy" ]]
    then
        echo "processing .ipynb"
        find ${REPO_NAME} -type f -name "*.ipynb" -exec mv -n {} ./${LOCAL_DIR}/ \;
        echo "define '${KERNEL}' kernel..."
        make_kernel
        echo "making '${KERNEL}' default kernel..."
        select_kernel ${LOCAL_DIR}
    else
        echo "processing .qmd"
        find ${REPO_NAME} -type f -name "*.qmd" -exec mv -n {} ./${LOCAL_DIR}/ \;
    fi
    rm -rf ./${REPO_NAME}
)

main
