FROM ghcr.io/scilifelabdatacentre/serve-jupyterlab:231124-1427

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG TARGETARCH

# Configure  environment
USER root

COPY scripts/install_scanpy.sh /tmp/scripts/

RUN /tmp/scripts/install_scanpy.sh && \
    mamba install --yes --channel conda-forge conda-lock

COPY conda/scanpy-environment-2025.yaml conda/scanpy-conda-lock.yaml /tmp/conda/

RUN conda-lock install --name scanpy /tmp/conda/scanpy-conda-lock.yaml && \
    mamba run -n scanpy python -m ipykernel install --name scanpy --display-name "scanpy" && \
    mamba clean --all -f -y && \
    chown -R ${NB_UID}:${NB_GID} "${CONDA_DIR}" && \
    chown -R ${NB_UID}:${NB_GID} "${CONDA_DIR}/envs" && \
    chown -R ${NB_UID}:${NB_GID} "/home/${NB_USER}"

# Configure container start
COPY --chown="${NB_UID}:${NB_GID}" --chmod=0755 scripts/scanpy-start-script.sh ${HOME}/start-script.sh
COPY --chown="${NB_UID}:${NB_GID}" --chmod=0755 scripts/download-labs.sh ${HOME}/download-labs.sh

USER ${NB_UID}

EXPOSE 8888

ENTRYPOINT ["tini", "-g", "--", "./start-script.sh"]
