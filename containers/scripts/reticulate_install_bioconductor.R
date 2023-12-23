library(reticulate)

conda_create(
    envname = "bioc",
    python_version = "3.8",
    channel = c("conda-forge", "bioconda", "anaconda"),
    conda = "/opt/conda/bin/conda"
)

virtualenv_create("/opt/venv/scanorama")
use_virtualenv("/opt/venv/scanorama")
py_install("scanorama==1.7.3")
