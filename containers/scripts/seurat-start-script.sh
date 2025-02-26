#!/bin/bash

/usr/local/conda/bin/conda init bash > /dev/null
source ${HOME}/.bashrc

conda activate seurat

export R_HOME=/usr/local/conda/envs/seurat/lib/R
echo "Using R: $(which R)"

echo "Starting server..."

/usr/lib/rstudio-server/bin/rserver \
  --auth-none=${DISABLE_AUTH} \
  --rsession-ld-library-path=${CONDA_PREFIX}/lib \
  --rsession-which-r=$(which R) \
  --server-daemonize=0 \
  --server-user=${USER} \
  --server-working-dir=${HOME} \
  --www-address=0.0.0.0 \
  --www-port=8787
