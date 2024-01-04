#!/bin/bash

set -e

if [ -z "$QUARTO_VERSION" ]; then
    echo "Error: QUARTO_VERSION is not set or is null."
    exit 1
fi

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
    libgdal-dev \
    libudunits2-dev \
    curl

## Install Quarto
curl -o quarto-linux-amd64.deb -L https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb \
    && apt-get install -y ./quarto-linux-amd64.deb \
    && rm -rf ./quarto-linux-amd64.deb \

## Install R packages from CRAN
install2.r --error --skipinstalled -n "$NCPUS" \
    dplyr \
    here \
    htmlTable \
    leaflet \
    lubridate \
    readxl \
    remotes \
    stringr \
    tidyr \
    writexl \
    yaml

## Clean up
rm -rf /var/lib/apt/lists/*
rm -rf /tmp/downloaded_packages

## Strip binary installed libraries from RSPM
## https://github.com/rocker-org/rocker-versioned2/issues/340
strip /usr/local/lib/R/site-library/*/libs/*.so
