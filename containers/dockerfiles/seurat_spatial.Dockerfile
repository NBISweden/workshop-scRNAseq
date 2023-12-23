FROM rocker/tidyverse:4.3.2

ARG NCPUS=${NCPUS:--1}

COPY scripts/install_seurat_spatial.sh \
    scripts/reticulate_install_seurat_spatial.R \
    /rocker_scripts/

RUN /rocker_scripts/install_seurat_spatial.sh

RUN Rscript /rocker_scripts/reticulate_install_seurat_spatial.R

## Make user rstudio owner of the Python environment dirs
RUN chgrp -R rstudio /opt/conda && \
    sudo chmod 770 -R /opt/conda
