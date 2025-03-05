import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scib_metrics.benchmark import Benchmarker
import scipy

timestamp= '18022025_185142'

# Read data
adata = sc.read_h5ad('data/gut_data/gut_hs_all_datasets_immune_annotated_AM_18022025_185142_raw.h5ad')
if not np.issubdtype(adata.X.dtype, np.floating):
    # If X is a sparse matrix
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.astype(np.float32)
    else:
        adata.X = adata.X.astype(np.float32)

adata_ref = sc.read_h5ad('data/gut_data/Integrated_4_datasets_05042024.h5ad') 
adata_ref.layers["counts"] = adata_ref.X.copy()

sc.pp.highly_variable_genes(
    adata_ref,
    flavor = "seurat_v3",
    n_top_genes = 7000,
    layer = "counts",
    batch_key = "Library_Preparation_Protocol",
    subset = True,
    span = 1
)

genes_to_keep = adata_ref.var_names.intersection(adata.var_names)
adata = adata[:, genes_to_keep].copy()

del adata_ref

# Run benchmarker: batch_key = 'sample_id'
bm = Benchmarker(
    adata,
    batch_key='sample_id',
    label_key='cell_states',
    embedding_obsm_keys=['X_scANVI', 'X_scVI', 'X_pca', 'umap_uncorrected'],
    n_jobs=-1,
)
bm.benchmark()

# Save metrics to DataFrame
metrics_df = pd.DataFrame(bm.get_results())
metrics_df.to_csv(f'/figures/benchmark_metrics_immune_sample_id_{timestamp}.csv')

# Create and save plot
fig = plt.figure(figsize=(10, 6))
bm.plot_results_table(min_max_scale=False)
plt.savefig(f'/figures/benchmark_metrics_immune_sample_id_{timestamp}.png', dpi=300, bbox_inches='tight')
plt.close()