#!/bin/bash

set -e

## Build ARGs
NCPUS=${NCPUS:--1}
QUARTO_VERSION="1.3.450"

## Function to install apt packages only if they are not installed
function apt_install() {
    if ! dpkg -s "$@" >/dev/null 2>&1; then
        if [ "$(find /var/lib/apt/lists/* | wc -l)" = "0" ]; then
            apt-get update
        fi
        apt-get install -y --no-install-recommends "$@"
    fi
}

apt_install \
    libhdf5-dev \
    libglpk-dev \
    libxt6 \
    patch \
    vim \
    curl

## Install quarto cli
curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
    && apt-get install -y ./quarto-linux-amd64.deb \
    && rm -rf ./quarto-linux-amd64.deb \

## Install Miniconda
wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/rstudio/miniconda.sh
/bin/bash /home/rstudio/miniconda.sh -b -p /opt/conda
rm -rf /home/rstudio/miniconda.sh

## Init conda for root and rstudio users
/opt/conda/bin/conda init bash
su - rstudio -c "/opt/conda/bin/conda init bash"

## Create convenience symlink for installBioc.r
ln -sf \
    "${R_HOME}/site-library/littler/examples/installBioc.r" \
    /usr/local/bin/installBioc.r

## Install packages from CRAN
install2.r --error --skipinstalled -n "$NCPUS" \
    hdf5r \
    umap \
    Seurat

## Install packages with BiocManager (https://stackoverflow.com/a/62456026)
installBioc.r --error --skipinstalled -n "$NCPUS" \
    Biobase \
    glmGamPoi

## Install packages from GitHub
installGithub.r --update FALSE \
    https://github.com/renozao/xbioc/tree/1354168bd7e64be4ee1f9f74e971a61556d75003 \
    https://github.com/meichendong/SCDC/tree/890c604eebd7fffa4a08d7344fbd516df6efcf8d \
    https://github.com/satijalab/seurat-data/tree/d6a8ce61ccb21a3b204f194d07009772c822791d

## Get brain data
wget seurat.nygenome.org/src/contrib/stxBrain.SeuratData_0.1.1.tar.gz \
    && R CMD INSTALL stxBrain.SeuratData_0.1.1.tar.gz \
    && rm -rf stxBrain.SeuratData_0.1.1.tar.gz
    
## Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## Strip binary installed libraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
strip /usr/local/lib/R/site-library/*/libs/*.so
