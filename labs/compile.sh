#!/bin/bash

# COMPILE QMD/IPYNB
# Replaces meta variables with text
# .qmd files with knitr engine remain as .qmd files
# .qmd files with jupyter engine is converted to .ipynb
# run in root of repo

pwd
cd labs/compiled/seurat
quarto render ../../seurat/seurat_01_qc.qmd --to markdown-header_attributes --output seurat_01_qc.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_02_dimred.qmd --to markdown-header_attributes --output seurat_02_dimred.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_03_integration.qmd --to markdown-header_attributes --output seurat_03_integration.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_04_clustering.qmd --to markdown-header_attributes --output seurat_04_clustering.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_05_dge.qmd --to markdown-header_attributes --output seurat_05_dge.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_06_celltyping.qmd --to markdown-header_attributes --output seurat_06_celltyping.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_07_spatial.qmd --to markdown-header_attributes --output seurat_07_spatial.qmd --metadata engine:markdown
quarto render ../../seurat/seurat_08_trajectory.qmd --to markdown-header_attributes --output seurat_08_trajectory.qmd --metadata engine:markdown

cd ../bioc
quarto render ../../bioc/bioc_01_qc.qmd --to markdown-header_attributes --output bioc_01_qc.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_02_dimred.qmd --to markdown-header_attributes --output bioc_02_dimred.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_03_integration.qmd --to markdown-header_attributes --output bioc_03_integration.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_04_clustering.qmd --to markdown-header_attributes --output bioc_04_clustering.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_05_dge.qmd --to markdown-header_attributes --output bioc_05_dge.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_06_celltyping.qmd --to markdown-header_attributes --output bioc_06_celltyping.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_07_spatial.qmd --to markdown-header_attributes --output bioc_07_spatial.qmd --metadata engine:markdown
quarto render ../../bioc/bioc_08_trajectory.qmd --to markdown-header_attributes --output bioc_08_trajectory.qmd --metadata engine:markdown

cd ../scanpy
quarto render ../../scanpy/scanpy_01_qc.qmd --to markdown-header_attributes --output scanpy_01_qc.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_02_dimred.qmd --to markdown-header_attributes --output scanpy_02_dimred.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_03_integration.qmd --to markdown-header_attributes --output scanpy_03_integration.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_04_clustering.qmd --to markdown-header_attributes --output scanpy_04_clustering.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_05_dge.qmd --to markdown-header_attributes --output scanpy_05_dge.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_06_celltyping.qmd --to markdown-header_attributes --output scanpy_06_celltyping.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_07_spatial.qmd --to markdown-header_attributes --output scanpy_07_spatial.qmd --metadata engine:markdown
quarto render ../../scanpy/scanpy_08_trajectory.qmd --to markdown-header_attributes --output scanpy_08_trajectory.qmd --metadata engine:markdown

cd ../../
quarto convert labs/scanpy/scanpy_01_qc.qmd --output docs/compiled/scanpy/scanpy_01_qc.ipynb
quarto convert labs/scanpy/scanpy_02_dimred.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_02_dimred.qmd
quarto convert labs/scanpy/scanpy_03_integration.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_03_integration.qmd
quarto convert labs/scanpy/scanpy_04_clustering.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_04_clustering.qmd
quarto convert labs/scanpy/scanpy_05_dge.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_05_dge.qmd
quarto convert labs/scanpy/scanpy_06_celltyping.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_06_celltyping.qmd
quarto convert labs/scanpy/scanpy_07_spatial.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_07_spatial.qmd
quarto convert labs/scanpy/scanpy_08_trajectory.qmd --to markdown-header_attributes --output docs/compiled/scanpy/scanpy_08_trajectory.qmd
