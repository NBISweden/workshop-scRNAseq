FROM ghcr.io/scilifelabdatacentre/serve-jupyterlab:231124-1427

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG TARGETARCH

ENV HOME=/home/${NB_USER}

# Configure  environment
USER root

COPY scripts/install_scanpy.sh /tmp/scripts/

RUN /tmp/scripts/install_scanpy.sh && \
    mamba install --yes --channel conda-forge conda-lock

COPY conda/scanpy-environment-2025.yaml conda/scanpy-conda-lock.yaml /tmp/conda/

RUN conda-lock install --name scanpy /tmp/conda/scanpy-conda-lock.yaml && \
    mamba run -n scanpy python -m ipykernel install --name scanpy --display-name "scanpy" && \
    mamba clean --all -f -y && \
    chown -R ${NB_UID}:${NB_GID} ${CONDA_DIR} && \
    chown -R ${NB_UID}:${NB_GID} ${CONDA_DIR}/envs && \
    chown -R ${NB_UID}:${NB_GID} ${HOME}

RUN mkdir -p ${HOME}/work \
    && chown -R ${NB_UID}:${NB_GID} ${HOME} \
    && cp ~/.bashrc ~/.bash_profile

# Configure container start
COPY scripts/scanpy-start-script.sh ${HOME}/work/start-script.sh
COPY scripts/download-labs.sh ${HOME}/work/download-labs.sh

RUN chown -R ${NB_UID}:${NB_GID} ${HOME} \
    && chmod a+x ${HOME}/work/start-script.sh \
    && chmod a+x ${HOME}/work/download-labs.sh

USER ${NB_USER}

WORKDIR ${HOME}/work

EXPOSE 8888

ENTRYPOINT ["tini", "-g", "--", "./start-script.sh"]
