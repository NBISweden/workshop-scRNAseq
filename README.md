# workshop-scRNAseq

![](https://github.com/NBISweden/workshop-scRNAseq/actions/workflows/docker-publish-site.yaml/badge.svg) ![](https://github.com/NBISweden/workshop-scRNAseq/actions/workflows/docker-publish-toolkits.yaml/badge.svg)

This repo contains the course material for NBIS workshop **Single Cell RNA-Seq Data Analyses**. The rendered view of this repo is available [here](https://nbisweden.github.io/workshop-scRNAseq/).

## Contributing

To add or update contents of this repo (for collaborators), first clone the repo, create a new branch, make changes/updates as needed, stage the changes, commit it and push the new branch to GitHub. Then, on GitHub, send a pull request to master.

```
git clone --depth 1 --single-branch --branch master https://github.com/nbisweden/workshop-scrnaseq.git
git checkout -b <branch-name>
git add .
git commit -m "I did this and that"
git push -u origin <branch_name>
```

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
docker run --rm --platform=linux/amd64 -p 8888:8888 -v ${PWD}:/home/jovyan/work ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256
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
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "seurat" "~/work/labs"

# for bioconductor
~/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "bioc" "~/work/labs"

# for scanpy
~/work/download-labs.sh "https://github.com/NBISweden" "workshop-scRNAseq" "compiled/labs" "scanpy" "~/work/labs"
```

## Render labs

Instructions to render the `.qmd` files to `.html`.

- For Seurat labs

```
# r/seurat
docker run --rm -ti --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311 quarto render /work/labs/seurat/seurat_01_qc.qmd

# r/bioc
docker run --rm -ti --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311 quarto render /work/labs/bioc/bioc_01_qc.qmd

# python/scanpy
docker run --rm -ti --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256 quarto render /work/labs/scanpy/scanpy_01_qc.qmd
```

- Successfully rendered outputs are moved to `docs` folder and chunks are cached under `_freeze`.

## Scripts

To render all qmd files in the repo to `docs/` as html output, run

```
bash scripts/render.sh all
```

The `render.sh` script has options to render all of the docs/scripts that are needed, but also to render different parts. The options are:

* all - run all the steps
* seurat - render all seurat labs
* bioc - render all bioc labs
* scanpy - render all scanpy labs
* spatial - render all 3 spatial labs.
* site - render all site stuff like lectures, contents, schedule etc.
* compile - compile labs into Rmd/ipynb


To compile all qmds into `compiled/labs` as qmds and ipynb with evaluated meta variables, can also be run directly with the compile script using:

```
bash scripts/compile.sh
```

---

**2025** • NBIS • SciLifeLab

