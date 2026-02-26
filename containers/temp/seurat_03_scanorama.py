import pyarrow.ipc as ipc
import pyarrow as pa
import numpy as np
import os
import scanorama

# The directory where R saved our Arrow files
data_dir = "/home/jovyan/work/labs/data/scanorama_arrow_data"

datasets = []
genes_list = []

# We have 8 datasets, numbered 1 to 8 based on the R loop
num_datasets = 8

print("Loading data from Arrow IPC files...")

for i in range(1, num_datasets + 1):
    assay_file = os.path.join(data_dir, f"assay_{i}.arrow")
    gene_file = os.path.join(data_dir, f"genes_{i}.arrow")
    
    # 1. Read the Dense Matrix
    # Memory mapping makes this read almost instantly
    with pa.memory_map(assay_file, 'r') as source:
        arrow_table = ipc.RecordBatchStreamReader(source).read_all()
        # Convert to a 2D numpy array. Going through pandas is highly optimized in Arrow.
        np_matrix = arrow_table.to_pandas().to_numpy()
        
    # 2. Read the Gene List
    with pa.memory_map(gene_file, 'r') as source:
        gene_table = ipc.RecordBatchStreamReader(source).read_all()
        # Extract the 'genes' column and convert it directly to a Python list
        genes = gene_table.column('genes').to_pylist()
        
    # IMPORTANT: Scanorama expects matrices in the shape (Cells, Genes).
    # If your R matrices were (Genes, Cells) - which is common in Seurat - 
    # you must transpose them here by uncommenting the next line:
    # np_matrix = np_matrix.T 
        
    datasets.append(np_matrix)
    genes_list.append(genes)
    
    print(f"  Loaded dataset {i}: Matrix shape {np_matrix.shape}, Genes count {len(genes)}")

print("\nStarting Scanorama integration...")

# 3. Run Scanorama
# By providing both the datasets and the genes_list, Scanorama will automatically 
# figure out the intersecting genes across all 8 batches before integrating.
integrated_matrices, corrected_genes = scanorama.integrate(datasets, genes_list)

print("Integration complete!")

# 4. Save the results back to Arrow so R can read them
print("Saving integrated embeddings back to Arrow...")
for i, int_matrix in enumerate(integrated_matrices):
    # Convert the resulting numpy array back to an Arrow Table
    # (Since Scanorama outputs embeddings, usually 100 dimensions per cell)
    col_names = [f"dim_{d+1}" for d in range(int_matrix.shape[1])]
    arrays = [pa.array(int_matrix[:, d]) for d in range(int_matrix.shape[1])]
    result_table = pa.Table.from_arrays(arrays, names=col_names)
    
    output_file = os.path.join(data_dir, f"integrated_{i+1}.arrow")
    with pa.OSFile(output_file, 'wb') as sink:
        with pa.RecordBatchStreamWriter(sink, result_table.schema) as writer:
            writer.write_table(result_table)

print("All done. You can now load the integrated datasets back into R.")
