# Contributing

To add or update contents of this repo, first clone the repo, create a new branch (name the branch according to the changes you are going to introduce), make changes/updates as needed, stage the changes, commit and push the new branch to GitHub. Then, on GitHub, open a pull request towards `main`.

```
git clone --depth 1 --single-branch --branch main https://github.com/nbisweden/workshop-scrnaseq.git
git checkout -b <branch-name>

# make your changes

git add <changed_files>
git commit -m "short descriptive message"
git push -u origin <branch_name>
```

The labs use Quarto's _shortcodes_ for templating, with metadata located in `docs/labs/_metadata.yml`. These are snippets of text identified by a _tag_ of the form `{{< meta >}}` that allow to insert the same content to the different toolkits. The lab templates are stored in `docs/labs/<toolkit>`. To get the final version of the labs you need to compile the labs and the results are stored in `compiled/labs/<toolkit>`.

To make changes to the labs, the workflow is as follows:
- Add/modify metadata snippets.
- Add the new snippet tag to all labs in `docs/labs`.
- Add/modify code blocks in each toolkit.
- Compile the labs.
- Test your changes on the compiled labs.
- Render the site and serve it locally.
- (_Course leaders only!_) Publish the site.

Once all changes are tested, you need to [render the site](#render-site). This step takes `.qmd` files and renders `.html` files for the site. The site is located in `docs/_site`.

> **IMPORTANT:** Before rendering any files, make sure `docs/_site` is _synced_ with the contents of `gh-pages` branch. See [site](#site) section for more details.

`docs/_site` folder is ignored by git. This means that it will be visible _as is_ in all your branches. This is important, because when the site is published, it is published from this location. In other words, this should reflect the content of `gh-pages` branch.

## Site

This section contains guidelines on how to make changes to site-related files (_i.e._ `compiled/` and `docs/` folders). Read first the [suggested workflow](#suggested-workflow) section.

> **IMPORTANT:** If you made changes to the labs (`docs/labs`) you need to [compile the labs](#compile-labs) first, before rendering the site.

### Suggested workflow

- Make sure your feature branch is clean (no uncommited changes).
- Switch to `gh-pages` branch and run `git pull` to make sure it is up to date.
- Switch back to the feature branch.
- If you don't have a folder `docs/_site` you should create one.
  ```sh
  mkdir -p docs/_site
  ```
- If you already have a folder `docs/_site`, remove all its content (but not the folder itself).
  ```sh
  rm -rf docs/_site/.
  ```
- Update the folder `docs/_site` with the content of `gh-pages` branch.
  ```sh
  git archive gh-pages | tar -x -C docs/_site
  ```
- Verify that `docs/_site` on your feature branch is now up to date.
- [Render](#render-site) any code files that were changed.
- Serve the site locally and check the changes by running the command below. Then in the browser, navigate to `localhost:8000`, and check the site content.
  ```sh
  cd docs/_site
  python3 -m http.server 8000
  ```
- (_Course leaders only!_) From the `docs` folder, publish the content of `docs/_site`.
  ```sh
  cd ..  # $PROJECT_DIR/docs
  quarto publish gh-pages --no-render
  ```

> **IMPORTANT:** This workflow works with `archive` as well, since we tar the `gh-pages` branch content, and as long as `gh-pages:archive` folder exists, it will be included in the tar, and copied to `docs/_site`. Rendering will not interfere with this folder.

### Compile labs

To compile all `.qmd` into `compiled/labs` as `.qmd` or `.ipynb` with evaluated meta variables, can be achieved directly with the compile script:

```
bash scripts/compile.sh [seurat|bioc|scanpy|all]
```

### Render site

To render `.qmd` files (site files as well as labs) in the repo to `docs/_site` as `.html` output, run the command below, choosing the corresponding argument, depending on your changes.

> **WARNING:** Rendering the _labs_ takes several minutes because it executes all the code cells.

```
bash scripts/render.sh [seurat|bioc|scanpy|spatial|site|compile|all]
```

If you made changes to only selected labs, you can render them individually as shwon below:

```
# r/seurat
docker run --rm -ti --platform=linux/amd64 -u root -v ${PWD}:/home/jovyan/work --entrypoint /usr/local/conda/bin/conda ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311 run -n seurat quarto render /home/jovyan/work/docs/labs/seurat/seurat_<lab_name>.qmd

# r/bioc
docker run --rm -ti --platform=linux/amd64 -u root -v ${PWD}:/home/jovyan/work --entrypoint /usr/local/conda/bin/conda ghcr.io/nbisweden/workshop-scrnaseq-seurat:20250320-2311 run -n seurat quarto render /home/jovyan/work/docs/labs/bioc/bioc_<lab_name>.qmd

# python/scanpy
docker run --rm -ti --platform=linux/amd64 -u jovyan -v ${PWD}:/work --entrypoint /opt/conda/bin/conda ghcr.io/nbisweden/workshop-scrnaseq-scanpy:20250325-2256 run -u scanpy quarto render /work/labs/scanpy/scanpy_<lab_name>.qmd
```

Successfully rendered outputs are moved to the `docs/_site` folder and chunks are cached under `docs/_freeze`. These folders are gitignored.

### Publish site

The site content to be published is in `docs/_site`. Verify the content locally by running
  ```sh
  cd docs/_site
  python3 -m http.server 8000
  ```
and in the browser, navigate to `localhost:8000`.

> **IMPORTANT:** `quarto publish` _erases_ the content of `gh-pages` branch and replaces it with the content of `docs/_site` folder. Make sure you updated it as described in the [suggested workflow](#suggested-workflow). In particular, double-check that the `archive/` folder is synced.

To serve the site from the feature branch, run
  ```sh
  cd docs
  quarto publish gh-pages --no-render
  ```

