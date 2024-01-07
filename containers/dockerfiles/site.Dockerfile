FROM rocker/r-base:4.3.0

LABEL Description="Docker image for NBIS workshop-scrnaseq"
LABEL Maintainer="roy.francis@nbis.se"
LABEL org.opencontainers.image.source="https://github.com/NBISweden/workshop-scrnaseq"

ARG QUARTO_VERSION="1.3.450"
ARG NCPUS=${NCPUS:--1}

COPY scripts/install_site.sh /rocker_scripts/

RUN /rocker_scripts/install_site.sh

CMD ["quarto", "render"]
