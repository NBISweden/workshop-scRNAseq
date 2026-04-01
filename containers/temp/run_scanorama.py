import pyarrow.ipc as ipc
import pyarrow as pa
import scanorama
import os
import numpy as np
import argparse

# 1. Set up command-line arguments
parser = argparse.ArgumentParser(description="Run Scanorama via Arrow IPC")
parser.add_argument("--basedir", type=str, required=True, help="Data folder path")
parser.add_argument("--prefix", type=str, required=True, help="Data folder prefix (e.g., 'counts')")
parser.add_argument("--num_datasets", type=int, default=8, help="Number of batches")
args = parser.parse_args()

data_dir = os.path.join(args.basedir, f"scanorama_data_{args.prefix}")
datasets = []
genes_list = []

print(f"--- Starting Scanorama for: {args.prefix} ---")

# 2. Load Data
for i in range(1, args.num_datasets + 1):
    with pa.memory_map(os.path.join(data_dir, f"assay_{i}.arrow"), 'r') as s:
        datasets.append(ipc.RecordBatchStreamReader(s).read_all().to_pandas().to_numpy())
    with pa.memory_map(os.path.join(data_dir, f"genes_{i}.arrow"), 'r') as s:
        genes_list.append(ipc.RecordBatchStreamReader(s).read_all().column('genes').to_pylist())

# 3. Integrate
integrated_matrices, _ = scanorama.integrate(datasets, genes_list)

# 4. Save Results
for i, int_matrix in enumerate(integrated_matrices):
    arrays = [pa.array(int_matrix[:, d]) for d in range(int_matrix.shape[1])]
    names = [f"dim_{d+1}" for d in range(int_matrix.shape[1])]
    table = pa.Table.from_arrays(arrays, names=names)
    
    with pa.OSFile(os.path.join(data_dir, f"integrated_{i+1}.arrow"), 'wb') as sink:
        with pa.RecordBatchStreamWriter(sink, table.schema) as writer:
            writer.write_table(table)

print(f"--- Finished {args.prefix} ---")
