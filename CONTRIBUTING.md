## Contributing

To add or update contents of this repo, first clone the repo, create a new branch, make changes/updates as needed, stage the changes, commit it and push the new branch to GitHub. Then, on GitHub, send a pull request to master.

```
git clone --depth 1 --single-branch --branch master https://github.com/nbisweden/workshop-scrnaseq.git
git checkout -b <branch-name>

# make your changes

git add <changed_files>
git commit -m "short descriptive message"
git push -u origin <branch_name>
```

The labs use Quarto's _shortcodes_ for templating, with metadata located in `docs/labs/_metadata.yml`. These are snippets of text identified by a _tag_ of the form `{{< meta >}}` and allow to insert the same content to the different toolkits. The lab templates are stored in `docs/labs/<toolkit>`. To get the final version of the labs you need to compile the labs and the results are stored in `compiled/labs/<toolkit>`.

To make changes to the labs, the workflow is as follows:
- Add/modify metedata snippets.
- Add the new snippet tag to all labs in `docs/labs`.
- Add/modify code blocks in all toolkits.
- Compile labs.
- Test your changes.

Once all changes are tested, you need to render the site. This step takes `.qmd` files and renders `.html` files for the site. The site is located in `docs/_site`.

## Compile labs

To compile all `.qmd` into `compiled/labs` as `.qmd` and `.ipynb` with evaluated meta variables, can be run directly with the compile script using:

```
bash scripts/compile.sh [seurat|bioc|scanpy|all]
```

## Render site

To render all `.qmd` files (site files and labs) in the repo to `docs/_site` as `.html` output, run the command below, choosing the corresponding argument, depending on your changes.

> **WARNING:** Rendering the labs takes several minutes because it executes all the code cells.

```
bash scripts/render.sh [seurat|bioc|scanpy|spatial|site|compile|all]
```

If you made changes to only selected labs, you can render them individually by running with the following commands:

```
# r/seurat
docker run --rm -ti --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/home/jovyan/work --entrypoint /usr/local/conda/bin/conda ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311 run -n seurat quarto render /home/jovyan/work/docs/labs/seurat/seurat_<lab_name>.qmd

# r/bioc
docker run --rm -ti --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/home/jovyan/work --entrypoint /usr/local/conda/bin/conda ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311 run -n seurat quarto render /home/jovyan/work/docs/labs/bioc/bioc_<lab_name>.qmd

# python/scanpy
docker run --rm -ti --platform=linux/amd64 -u 1000:1000 -v ${PWD}:/work --entrypoint /opt/conda/bin/conda ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256 run -u scanpy quarto render /work/labs/scanpy/scanpy_<lab_name>.qmd
```

Successfully rendered outputs are moved to `docs/_site` folder and chunks are cached under `docs/_freeze`. These folders are gitignored.

## Publish site

TBA

