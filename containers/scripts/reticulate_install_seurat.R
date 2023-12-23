library(reticulate)

conda_create(
    envname = "seurat",
    python_version = "3.8",
    packages = c("umap-learn==0.5.4"),
    channel = c("conda-forge", "bioconda", "anaconda"),
    conda = "/opt/conda/bin/conda"
)

virtualenv_create("/opt/venv/scanorama")
use_virtualenv("/opt/venv/scanorama")
py_install("scanorama==1.7.3")
