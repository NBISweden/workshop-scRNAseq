
cd /Users/asbj/local/misc/course_git/workshop-scRNAseq/labs/compiled/scanpy

conda activate scRNAseq2023_python

echo "Convert to notebook..."
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_01_qc.ipynb
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_02_dim_reduction.ipynb
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_03_integration.ipynb
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_04_clustering.ipynb
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_05_dge.ipynb
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_06_celltype.ipynb
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace scanpy_07_spatial.ipynb


conda activate nbconvert6
echo "Convert to html..."
jupyter nbconvert --to html_toc scanpy_01_qc.ipynb
jupyter nbconvert --to html_toc scanpy_02_dim_reduction.ipynb
jupyter nbconvert --to html_toc scanpy_03_integration.ipynb
jupyter nbconvert --to html_toc scanpy_04_clustering.ipynb
jupyter nbconvert --to html_toc scanpy_05_dge.ipynb
jupyter nbconvert --to html_toc scanpy_06_celltype.ipynb
jupyter nbconvert --to html_toc scanpy_07_spatial.ipynb


