library(reticulate)

conda_create(
    envname = "seurat-spatial",
    python_version = "3.10",
    channel = c("conda-forge", "bioconda", "anaconda"),
    conda = "/opt/conda/bin/conda"
)
