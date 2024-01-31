FROM jupyter/minimal-notebook:python-3.10

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG TARGETARCH

USER root

COPY scripts/install_scanpy.sh /scripts/

RUN /scripts/install_scanpy.sh

RUN mamba install --yes \
    --channel conda-forge \
    --channel bioconda \
    --channel pytorch \
    --channel anaconda \
    bbknn=1.6.0 \
    fa2=0.3.5 \
    gseapy=1.0.6 \
    matplotlib-venn=0.11.9 \
    openpyxl=3.1.2 \
    pybiomart=0.2.0 \
    scanorama=1.7.4 \
    scanpy=1.9.6 \
    scrublet=0.2.3 \
    scvi-tools=1.0.4 \
    seaborn=0.12.2 && \
    mamba clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    leidenalg==0.10.1 \
    louvain==0.8.1 \
    spatialde==1.1.3

USER ${NB_UID}
