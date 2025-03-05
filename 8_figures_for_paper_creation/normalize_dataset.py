import scanpy as sc
import numpy as np
import gc
import scipy.sparse as sp
from datetime import datetime

timestamp = datetime.now().strftime('%d%m%Y_%H%M%S')

# Parameters
chunk_size = 50000  # Reduced chunk size
target_sum = 1e6
file_path = 'data/gut_data/gut_hs_all_datasets_full_annotated_AM_03032025_141544_raw.h5ad'
output_path = f'data/gut_data/gut_hs_all_datasets_full_annotated_AM_{timestamp}_log.h5ad'

print(f"Reading original file: {file_path}")
print(f"Will write to: {output_path}")

# First, create a temporary AnnData in backed mode to hold structure
adata = sc.read_h5ad(file_path, backed='r')
print(f"Dataset shape: {adata.shape}")
print(f"Data type of X: {type(adata.X)}")

# Process in chunks and collect results
all_processed_chunks = []

for i in range(0, adata.shape[0], chunk_size):
    end = min(i + chunk_size, adata.shape[0])
    print(f"Processing chunk {i}-{end} ({i/adata.shape[0]*100:.1f}% complete)")
    
    # Load chunk into memory
    chunk = adata[i:end].to_memory()
    
    # Make sure X is in the right format (CSR matrix)
    if not sp.issparse(chunk.X):
        chunk.X = sp.csr_matrix(chunk.X)
    else:
        chunk.X = chunk.X.tocsr()
    
    # Normalize and log transform
    sc.pp.normalize_total(chunk, target_sum=target_sum, exclude_highly_expressed=True)
    sc.pp.log1p(chunk)
    
    # Store the processed chunk
    all_processed_chunks.append(chunk)
    
    # Clean up to free memory
    gc.collect()

print("All chunks processed. Combining results...")

# Combine all chunks into a single AnnData object
adata_combined = sc.concat(all_processed_chunks, merge="same")

# Write the final result
print(f"Writing final dataset to {output_path}")
adata_combined.write(output_path)

print("Complete!")