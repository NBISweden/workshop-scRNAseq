library(reticulate)

conda_create(
    envname = "bioc-spatial",
    python_version = "3.8",
    channel = c("conda-forge", "bioconda", "anaconda"),
    conda = "/opt/conda/bin/conda"
)