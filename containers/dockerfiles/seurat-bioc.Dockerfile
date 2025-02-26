FROM ghcr.io/scilifelabdatacentre/serve-rstudio:231030-1146

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG NCPUS=${NCPUS:--1}

USER root

COPY scripts/install_seurat_bioc.sh /tmp/scripts/

RUN /tmp/scripts/install_seurat_bioc.sh

COPY conda/seurat-bioc-environment-2025.yaml conda/seurat-bioc-conda-lock.yaml /tmp/conda/

RUN /usr/local/conda/bin/conda-lock install \
    --conda /usr/local/conda/bin/conda \
    --prefix /usr/local/conda/envs/seurat \
    /tmp/conda/seurat-bioc-conda-lock.yaml \
    && /usr/local/conda/bin/conda clean --all -f -y

# Configure container start
COPY --chown="${NB_UID}:${NB_GID}" --chmod=0755 scripts/seurat-start-script.sh ${HOME}/start-script.sh
COPY --chown="${NB_UID}:${NB_GID}" --chmod=0755 scripts/download-labs.sh ${HOME}/download-labs.sh

USER ${NB_UID}

WORKDIR ${HOME}

EXPOSE 8787

ENTRYPOINT ["/tini", "-g", "--", "./start-script.sh"]
