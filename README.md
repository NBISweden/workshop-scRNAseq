# workshop-scRNAseq

![](https://github.com/NBISweden/workshop-scRNAseq/actions/workflows/docker-publish-site.yaml/badge.svg) ![](https://github.com/NBISweden/workshop-scRNAseq/actions/workflows/docker-publish-toolkits.yaml/badge.svg)

This repo contains the course material for NBIS workshop **Single Cell RNA-Seq Data Analyses**. The rendered view of this repo is available [here](https://nbisweden.github.io/workshop-scRNAseq/).

## Environment

```
# for seurat and bioconductor labs
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311`

# for scanpy labs
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20260323-2301
```

## Run labs interactively (locally)

> To run the labs locally follow these [instructions](https://nbisweden.github.io/workshop-scRNAseq/other/docker.html) to install Docker Desktop / Colima, depending on your operating system.

> **IMPORTANT:** If you are using an Apple Silicon (M-chip) you need to follow the Colima [instructions](https://nbisweden.github.io/workshop-scRNAseq/other/docker.html#running-linux-x86_64-containers-on-apple-silicon)!

- Create a new directory and `cd` into it. You will mount this directory to `/home/jovyan/work` in your container so that you can save your work locally.

- To run Seurat or Bioconductor labs in RStudio

```
docker run --rm --platform=linux/amd64 -p 8787:8787 -v ${PWD}:/home/jovyan/work ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311
```

Open in browser: `http://localhost:8787/`

- To run Python labs in JupyterLab

```
docker run --rm --platform=linux/amd64 -p 8888:8888 -v ${PWD}:/home/jovyan/work ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20260323-2301

# Apple Silicon with Colima
docker run --rm --platform=linux/amd64 -p 8888:8888 -v scanpy-labs:/home/jovyan/work ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20260323-2301
```

Open in browser: `http://localhost:8888/lab` and use password `scrnaseq`

- In the container, start a terminal and run the command below to activate the respective environment.

```
# for seurat/bioconductor
conda activate seurat
```

For Scanpy we are using Pixi as environment manager and you do not need to activate the environment. If you are running a command from `/home/jovyan` you just need to prepend any command with `pixi run <CMD>`. In any other directory, you need to tell Pixi which manifest to use as shown below.

```
pixi run --frozen --manifest-path /home/jovyan/pixi.toml <CMD>
```

- To download the compiled labs for the respective toolkit, run the `download-labs.sh` command below provided in the container.

```
# for seurat
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "seurat" "work/labs"

# for bioconductor
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "bioc" "work/labs"

# for scanpy
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "scanpy" "work/labs"
```

---

**2026** • NBIS • SciLifeLab

