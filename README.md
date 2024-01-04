# workshop-scRNAseq

This repo contains the course material for NBIS workshop **Single Cell RNA-Seq Data Analyses**. The rendered view of this repo is available [here](https://nbisweden.github.io/workshop-scrnaseq/).

## Contributing

To add or update contents of this repo (for collaborators), first clone the repo.

```
git clone --depth 1 --single-branch --branch master https://github.com/nbisweden/workshop-scrnaseq.git
```

Make changes/updates as needed. Add the changed files. Commit it. Then push the repo back.

```
git add .
git commit -m "I did this and that"
git push origin
```

## Environment

```
# for seurat labs
docker pull ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat-r4.3.0
docker pull ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat_spatial-r4.3.0

# for bioconductor labs
docker pull ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor-r4.3.0
docker pull ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor_spatial-r4.3.0

# for python labs
docker pull ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy-py3.10
docker pull ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy_spatial-py3.10
```

## Run labs interactively (locally)

- Launch docker container in the project's root folder
- To run Seurat or Bioconductor labs in RStudio

```
docker run --rm -ti --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat-r4.3.0
docker run --rm -ti --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor-r4.3.0
```

- Open in browser: `http://localhost:8787/`, login: rstudio, pass: scrnaseq
- Navigate to `/home/rstudio/workdir/labs` and open qmd files

- To run Python labs in JupyterLab

```
docker run --rm -ti --platform=linux/amd64 -p 8888:8888 -v $PWD:/home/jovyan/workdir ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy-py3.10
```

- Open in browser: `http://127.0.0.1:8888/lab?token=xxxx` (Use exact token from terminal on launch)
- Navigate to `/home/jovyan/workdir/compiled/scanpy` and open .ipynb files

## Render labs

Instructions to render the `.qmd` files to `.html`.

- For Seurat labs

```
# r/seurat
docker run --rm -ti --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work ghcr.io/nbisweden/workshop-scrnaseq:2024-seurat-r4.3.0 quarto render /work/labs/seurat/seurat_01_qc.qmd

# r/bioc
docker run --rm -ti --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/work ghcr.io/nbisweden/workshop-scrnaseq:2024-bioconductor-r4.3.0 quarto render /work/labs/bioc/bioc_01_qc.qmd

# python/scanpy
docker run --rm -ti --platform=linux/amd64 -p 8888:8888 -v $PWD:/work ghcr.io/nbisweden/workshop-scrnaseq:2024-scanpy-py3.10 quarto render /work/labs/scanpy/scanpy_01_qc.qmd
```

- Successfully rendered outputs are moved to `docs` folder and chunks are cached under `_freeze`.

---

**2024** • NBIS • SciLifeLab
