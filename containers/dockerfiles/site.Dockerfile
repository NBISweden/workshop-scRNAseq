# DOCKER FILE FOR WORKSHOP-SCRNASEQ SITE

FROM rocker/r-base:4.3.0

LABEL Description="Docker image for NBIS workshop-scrnaseq"
LABEL Maintainer="roy.francis@nbis.se"
LABEL org.opencontainers.image.source="https://github.com/NBISweden/workshop-scrnaseq"

ARG QUARTO_VERSION="1.3.450"
ARG NCPUS=${NCPUS:--1}

COPY scripts/install_site.sh /rocker_scripts/

RUN /rocker_scripts/install_site.sh

CMD quarto render

# build and run
# docker build --platform=linux/amd64 -t ghcr.io/nbisweden/workshop-scrnaseq:2024-site-r4.3.0 -f dockerfiles/site.Dockerfile .
# docker run --rm --platform=linux/amd64 -v $PWD:/qmd ghcr.io/nbisweden/workshop-scrnaseq:2024-site-r4.3.0 quarto render /qmd/index.qmd
