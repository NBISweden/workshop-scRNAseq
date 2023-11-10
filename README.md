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

## Rendering

Be in suitable conda environment or use a docker container with all the tools necessary. Then,

- Build the whole website

```
quarto render
```

- Build individual files

```
quarto render index.qmd
quarto render labs/seurat/seurat_01_qc.qmd
```

Successfully rendered outputs are moved to `docs` folder.

---

**2023** • NBIS • SciLifeLab