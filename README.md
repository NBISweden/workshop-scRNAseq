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
# for r labs
docker pull susrei/workshop-scrnaseq:2023-bioconductor-r4.3.0-conda

# for python labs
docker pull susrei/workshop-scrnaseq:2023-scanpy-py3.10
```

## Run labs

- Launch docker container in root of the repo

- To run R labs in RStudio

```
docker run --rm -ti --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir susrei/workshop-scrnaseq:2023-bioconductor-r4.3.0-conda
```

- Open in browser: `http://localhost:8787/`, login: rstudio, pass: scrnaseq
- Navigate to `/home/rstudio/workdir/labs` and open qmd files

- To run Python labs in Jupyter notebook

```
docker run --rm -ti --platform=linux/amd64 -p 8888:8888 -e PASSWORD=scrnaseq -v $PWD:/home/jovyan/workdir susrei/workshop-scrnaseq:2023-scanpy-py3.10
```

- Open in browser: `http://127.0.0.1:8888/lab?token=xxxx` (Use exact token from terminal on launch)
- Navigate to `/home/jovyan/workdir/compiled/scanpy` and open .ipynb files

## Render

Instructions to render the `.qmd` files to `.html`.

- For R labs

```
docker run --rm -ti --platform=linux/amd64 -p 8787:8787 -e PASSWORD=scrnaseq -v $PWD:/home/rstudio/workdir susrei/workshop-scrnaseq:2023-bioconductor-r4.3.0-conda quarto render /home/rstudio/workdir/labs/seurat/seurat_01_qc.qmd
```

- For Python labs

```
docker run --rm -ti --platform=linux/amd64 -p 8888:8888 -e PASSWORD=scrnaseq -v $PWD:/home/jovyan/workdir susrei/workshop-scrnaseq:2023-scanpy-py3.10 quarto render /home/jovyan/workdir/labs/scanpy/scanpy_01_qc.qmd
```

- Successfully rendered outputs are moved to `docs` folder and chunks are cached under `_freeze`.

---

**2023** • NBIS • SciLifeLab
