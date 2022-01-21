nb="compiled/scanpy/scanpy_07_spatial.ipynb"

conda activate python_spatial

echo "Convert to notebook..."
jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=1000 --inplace $nb

conda activate scRNAseq2022

echo "Convert to html..."
jupyter nbconvert  --to html_toc --ExecutePreprocessor.timeout=1000  $nb

