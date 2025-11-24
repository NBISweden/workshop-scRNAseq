# workshop-scRNAseq

![](https://github.com/NBISweden/workshop-scRNAseq/actions/workflows/docker-publish-site.yaml/badge.svg) ![](https://github.com/NBISweden/workshop-scRNAseq/actions/workflows/docker-publish-toolkits.yaml/badge.svg)

This repo contains the course material for NBIS workshop **Single Cell RNA-Seq Data Analyses**. The rendered view of this repo is available [here](https://nbisweden.github.io/workshop-scRNAseq/).

## Environment

```
# for seurat and bioconductor labs
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311`

# for scanpy labs
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256

# for optional spatial labs
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat_spatial-r4.3.0
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor_spatial-r4.3.0
docker pull --platform=linux/amd64 ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy_spatial-py3.10
```

## Run labs interactively (locally)

- Create a new directory and `cd` into it. You will mount this directory to `/home/jovyan/work` in your container so that you can save your work locally.

- To run Seurat or Bioconductor labs in RStudio

```
docker run --rm --platform=linux/amd64 -p 8787:8787 -v ${PWD}:/home/jovyan/work ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311
```

Open in browser: `http://localhost:8787/`

- To run Python labs in JupyterLab

```
docker run --rm --platform=linux/amd64 -p 8888:8888 -v ${PWD}:/home/jovyan/work/work ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256
```

Open in browser: `http://localhost:8888/lab` and use password `scrnaseq`

- In the container, start a terminal and run the command below to activate the respective conda environment.

```
# for seurat/bioconductor
conda activate seurat

# for scanpy
conda activate scanpy
```

- To download the compiled labs for the respective toolkit, run the `download-labs.sh` command below provided in the container.

```
# for seurat
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "seurat" "work/labs"

# for bioconductor
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "bioc" "work/labs"

# for scanpy
~/work/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "scanpy" "work/labs"
```

---

**2025** • NBIS • SciLifeLab

