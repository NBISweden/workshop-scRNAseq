FROM ghcr.io/scilifelabdatacentre/serve-jupyterlab-minimal:250213-1113

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG TARGETARCH

ENV HOME=/home/${NB_USER}

# Configure  environment
USER root

COPY scripts/install_scanpy.sh /tmp/scripts/

RUN /tmp/scripts/install_scanpy.sh

ENV PATH="${HOME}/.pixi/bin:$PATH"

COPY envs/scanpy/pixi.toml envs/scanpy/pixi.lock /home/jovyan/

# Configure container start
COPY scripts/scanpy-start-script.sh ${HOME}/start-script.sh
COPY scripts/download-labs.sh ${HOME}/download-labs.sh
COPY scripts/wait-for-env.sh /usr/local/bin/wait-for-env

RUN chmod a+x ${HOME}/start-script.sh && \
    chmod a+x ${HOME}/download-labs.sh && \
    chmod +x /usr/local/bin/wait-for-env && \
	mkdir -p /var/log/pixi && \
	chown -R jovyan:users /var/log/pixi

RUN mkdir -p /home/jovyan/.local/share/jupyter && \
    chown -R jovyan:users /home/jovyan

ENV PIXI_LOCKED=true
ENV PIXI_CACHE_DIR=$HOME/work/.cache/rattler/cache

USER ${NB_USER}

WORKDIR ${HOME}

EXPOSE 8888

ENTRYPOINT ["tini", "-g", "--", "./start-script.sh"]
