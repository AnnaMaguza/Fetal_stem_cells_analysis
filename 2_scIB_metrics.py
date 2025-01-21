import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scib_metrics.benchmark import Benchmarker
import scipy

timestamp= '20012025_132650'
# Read data
adata = sc.read_h5ad('integration_of_remapped_data/gut_hs_all_datasets_scVI_scANVI_epithelial_cellstates_AM_20012025_132650_raw.h5ad')
if not np.issubdtype(adata.X.dtype, np.floating):
    # If X is a sparse matrix
    if scipy.sparse.issparse(adata.X):
        adata.X = adata.X.astype(np.float32)
    else:
        adata.X = adata.X.astype(np.float32)

adata
# Run benchmarker
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(
    adata,
    flavor = "seurat_v3",
    n_top_genes = 7000,
    layer = "counts",
    batch_key = "library_construnction_and_layout",
    subset = True,
    span = 1
)

bm = Benchmarker(
    adata,
    batch_key='sample_id',
    label_key='cellstates_scANVI',
    embedding_obsm_keys=['X_scANVI', 'X_scVI', 'X_pca'],
    n_jobs=-1,
)
bm.benchmark()

# Save metrics to DataFrame
metrics_df = pd.DataFrame(bm.get_results())
metrics_df.to_csv(f'integration_of_remapped_data/benchmark_metrics_epithelial_sample_id_{timestamp}.csv')

# Create and save plot
fig = plt.figure(figsize=(10, 6))
bm.plot_results_table(min_max_scale=False)
plt.savefig(f'integration_of_remapped_data/benchmark_metrics_epithelial_sample_id_{timestamp}.png', dpi=300, bbox_inches='tight')
plt.close()
