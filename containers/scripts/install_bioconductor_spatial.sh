#!/bin/bash

set -e

## Build ARGs
NCPUS=${NCPUS:--1}

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
    libbz2-dev \
    libglpk-dev \
    libxt-dev \
    libhdf5-dev \
    patch \
    vim  

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

## Install packages with BiocManager (https://stackoverflow.com/a/62456026)
installBioc.r --error --skipinstalled -n "$NCPUS" \
    batchelor \
    biomaRt \
    scater \
    scran \
    SingleCellExperiment \
    SingleR \
    Spaniel 

## Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## Strip binary installed libraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
strip /usr/local/lib/R/site-library/*/libs/*.so
