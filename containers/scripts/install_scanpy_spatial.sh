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
    curl \
    g++ \
    gcc \
    gdebi \
    libfmt-dev

## Install quarto cli
curl -o quarto-linux-${TARGETARCH}.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH}.deb
gdebi --non-interactive quarto-linux-${TARGETARCH}.deb
rm -rf quarto-linux-${TARGETARCH}.deb

## Install Quarto Jupyter extension
python3 -m pip install jupyterlab-quarto

## Clean up
apt-get clean
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages
