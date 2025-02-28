#!/bin/bash

set -euo pipefail

## Build ARGs
NCPUS=${NCPUS:--1}
QUARTO_VERSION="1.3.450"
TINI_VERSION="v0.19.0"

## Install apt packages
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends \
    binutils \
    ca-certificates \
    curl \
    libglpk-dev \
    libhdf5-dev \
    libodbc1 \
    libxml2 \
    libxt6 \
    locales \
    patch \
    software-properties-common \
    vim

apt-get update && apt-get install -y git

echo "en_US.UTF-8 UTF-8" > /etc/locale.gen
locale-gen

add-apt-repository ppa:ubuntu-toolchain-r/test
apt-get -yq update
apt-get -yq install --no-install-recommends libstdc++6

## Install tini
curl -o /tini -sL https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini
chmod +x /tini

## Install quarto cli
curl -o quarto-linux-amd64.deb -sL https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
    && apt-get install -y ./quarto-linux-amd64.deb \
    && rm -rf ./quarto-linux-amd64.deb \

## Install Miniconda
curl -o ${HOME}/miniconda.sh -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ${HOME}/miniconda.sh -bfp /usr/local/conda
rm -rf ${HOME}/miniconda.sh

## Init conda for root and jovyan users
/usr/local/conda/bin/conda init bash
su - ${USER} -c "/usr/local/conda/bin/conda init bash"
source ${HOME}/.bashrc

## Install mamba
conda install --yes -n base -c conda-forge mamba
eval "$(mamba shell hook --shell bash)"
mamba shell init --shell bash
source ${HOME}/.bashrc

## Install conda-lock
mamba install --yes --channel conda-forge conda-lock

## Clean up
apt-get clean
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages
