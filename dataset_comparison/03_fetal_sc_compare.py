"""Deep comparison of fetal stem cells between subset (Mar-2025 run) and new big (Oct-2025 run).
Also examines adult stem cells in old (Healthy_colon_adult) vs new big.
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path
import json

OLD_PATH = "/Users/am336941/PhD/data/gut_data/Healthy_colon_adult.h5ad"
NEW_PATH = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
SUBSET_PATH = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
OUT.mkdir(exist_ok=True)
sc.settings.figdir = str(OUT / "figs")
(OUT / "figs").mkdir(exist_ok=True)

# ----------------------------------------------------------------
# 1. SUBSET (Mar 2025 fetal SC) — check leiden clusters
# ----------------------------------------------------------------
print("\n=== Step 1: Loading subset (Mar 2025 fetal SC analysis) ===")
sub = sc.read_h5ad(SUBSET_PATH)
print(f"shape: {sub.shape}, X dtype: {sub.X.dtype}")
print("leiden_cluster value counts:")
print(sub.obs["leiden_cluster"].value_counts(dropna=False).to_string())
print("\ncell_states:", sub.obs["cell_states"].value_counts().to_string())

# Save subset leiden + UMAP for plotting
if "X_umap" in sub.obsm:
    print(f"X_umap available: {sub.obsm['X_umap'].shape}")
    sub.obs[["leiden_cluster"]].assign(
        u1=sub.obsm["X_umap"][:,0], u2=sub.obsm["X_umap"][:,1]
    ).to_parquet(OUT / "subset_umap.parquet")

# Plot subset UMAP coloured by leiden
fig = sc.pl.umap(sub, color=["leiden_cluster", "Study_name", "age_group", "donor_id"],
                 ncols=2, save="_subset_fetalSC_overview.png", show=False, return_fig=True)

# ----------------------------------------------------------------
# 2. NEW BIG → extract fetal stem cells (the same 7308 cells)
# ----------------------------------------------------------------
print("\n=== Step 2: Loading new big and extracting fetal stem cells ===")
new = sc.read_h5ad(NEW_PATH)
print(f"new shape: {new.shape}")
print(f"new layers: {list(new.layers.keys())}")
print(f"new X dtype: {new.X.dtype}")
print(f"new layers['counts'] dtype:", new.layers["counts"].dtype if "counts" in new.layers else "n/a")

# stem cells overall (per cellstates_scANVI)
stem_mask = new.obs["cellstates_scANVI"] == "Stem cells"
fetal_mask = new.obs["age_group"].isin(["first trimester", "second trimester"])

print(f"\nNew stem cells: {stem_mask.sum()}")
print(f"  ↳ fetal: {(stem_mask & fetal_mask).sum()}")
print(f"  ↳ adult: {(stem_mask & (new.obs['age_group']=='adult')).sum()}")
print(f"  ↳ organoid (cell culture model): {(stem_mask & (new.obs['age_group']=='cell culture model')).sum()}")
print(f"  ↳ child stage: {(stem_mask & (new.obs['age_group']=='child stage')).sum()}")

# Subset's cells should all be in new big — verify
common = sub.obs.index.intersection(new.obs.index)
print(f"\nCells from subset present in new big: {len(common)} / {sub.n_obs}")

# ----------------------------------------------------------------
# 3. Compare expression for SHARED cells (subset ↔ new big)
# ----------------------------------------------------------------
print("\n=== Step 3: Compare expression for shared cells ===")
sub_sorted = sub[common, :].copy()
new_fetal_sc = new[common, :].copy()
print(f"sub_sorted: {sub_sorted.shape}, new_fetal_sc: {new_fetal_sc.shape}")

# Gene namespace
print(f"\nsub var index head: {list(sub_sorted.var.index[:5])}")
print(f"new var index head: {list(new_fetal_sc.var.index[:5])}")
gene_overlap = sub_sorted.var.index.intersection(new_fetal_sc.var.index)
print(f"gene overlap (var index): {len(gene_overlap)}")

# Compare total counts per cell
sub_X = sub_sorted.X
new_X = new_fetal_sc.X
if sp.issparse(sub_X):
    sub_total = np.asarray(sub_X.sum(axis=1)).flatten()
else:
    sub_total = sub_X.sum(axis=1)
if sp.issparse(new_X):
    new_total = np.asarray(new_X.sum(axis=1)).flatten()
else:
    new_total = new_X.sum(axis=1)

print(f"\nsubset total counts/cell — mean: {sub_total.mean():.1f}, median: {np.median(sub_total):.1f}")
print(f"new big total counts/cell — mean: {new_total.mean():.1f}, median: {np.median(new_total):.1f}")
print(f"per-cell total-counts correlation: {np.corrcoef(sub_total, new_total)[0,1]:.4f}")
print(f"per-cell ratio (new/old): mean={np.mean(new_total/np.maximum(sub_total,1)):.3f}, median={np.median(new_total/np.maximum(sub_total,1)):.3f}")

# Per-cell n_genes_by_counts
sub_ngenes = (sub_X > 0).sum(axis=1)
if sp.issparse(sub_X):
    sub_ngenes = np.asarray(sub_ngenes).flatten()
new_ngenes = (new_X > 0).sum(axis=1)
if sp.issparse(new_X):
    new_ngenes = np.asarray(new_ngenes).flatten()
print(f"\nsubset genes/cell — mean: {sub_ngenes.mean():.0f}")
print(f"new big genes/cell — mean: {new_ngenes.mean():.0f}")

# Per-cell expression correlation on shared genes
print(f"\n=== Per-cell expression correlation (Spearman on shared {len(gene_overlap)} genes) ===")
sub_shared = sub_sorted[:, gene_overlap].X
new_shared = new_fetal_sc[:, gene_overlap].X
if sp.issparse(sub_shared): sub_shared = sub_shared.toarray()
if sp.issparse(new_shared): new_shared = new_shared.toarray()

# Compute Pearson correlation of log1p(counts) per cell
from scipy.stats import pearsonr
n_sample = 200
idx = np.random.RandomState(0).choice(sub_shared.shape[0], min(n_sample, sub_shared.shape[0]), replace=False)
corrs = []
for i in idx:
    a = np.log1p(sub_shared[i].astype(np.float32))
    b = np.log1p(new_shared[i].astype(np.float32))
    if a.std() > 0 and b.std() > 0:
        c, _ = pearsonr(a, b)
        corrs.append(c)
corrs = np.array(corrs)
print(f"  mean per-cell Pearson r: {corrs.mean():.4f}")
print(f"  median: {np.median(corrs):.4f}")
print(f"  fraction with r > 0.95: {(corrs > 0.95).mean():.3f}")
print(f"  fraction with r > 0.99: {(corrs > 0.99).mean():.3f}")

# Save expression-summary stats
pd.DataFrame({
    "metric": ["total_counts_corr","mean_pearson_r","median_pearson_r",
               "frac_r_gt_0p95","frac_r_gt_0p99",
               "subset_mean_counts","new_mean_counts","ratio_new_over_old"],
    "value": [np.corrcoef(sub_total, new_total)[0,1],
              float(corrs.mean()), float(np.median(corrs)),
              float((corrs>0.95).mean()), float((corrs>0.99).mean()),
              float(sub_total.mean()), float(new_total.mean()),
              float(new_total.mean()/sub_total.mean())]
}).to_csv(OUT / "expression_comparison_stats.csv", index=False)

# ----------------------------------------------------------------
# 4. Re-cluster fetal stem cells from new big (just fetal SCs)
# ----------------------------------------------------------------
print("\n=== Step 4: Re-cluster fetal SC from new big ===")
new_fsc = new[stem_mask & fetal_mask, :].copy()
print(f"fetal SC from new big: {new_fsc.shape}")

# Use raw counts -> normalize -> HVG -> PCA -> neighbors -> umap -> leiden
new_fsc.X = new_fsc.layers["counts"].copy() if "counts" in new_fsc.layers else new_fsc.X.copy()
sc.pp.normalize_total(new_fsc, target_sum=1e4)
sc.pp.log1p(new_fsc)
sc.pp.highly_variable_genes(new_fsc, n_top_genes=3000, flavor="seurat", batch_key="batch")
new_fsc = new_fsc[:, new_fsc.var["highly_variable"]].copy()
sc.pp.scale(new_fsc, max_value=10)
sc.tl.pca(new_fsc, n_comps=30)
sc.pp.neighbors(new_fsc, n_neighbors=15, n_pcs=30)
sc.tl.umap(new_fsc)
sc.tl.leiden(new_fsc, resolution=0.5, key_added="leiden_05")
sc.tl.leiden(new_fsc, resolution=0.3, key_added="leiden_03")
sc.tl.leiden(new_fsc, resolution=0.8, key_added="leiden_08")
print("leiden 0.5:", new_fsc.obs["leiden_05"].value_counts().to_string())
print("leiden 0.3:", new_fsc.obs["leiden_03"].value_counts().to_string())
fig = sc.pl.umap(new_fsc, color=["leiden_03","leiden_05","leiden_08","Study_name","batch","age_group","library_preparation_protocol","donor_id","gut_region","Cell_cycle_phase"],
                 ncols=3, save="_newbig_fetalSC_recluster.png", show=False, return_fig=True)

# Quick: use the precomputed scVI/scANVI from new big (already integrated)
print("\n=== Step 5: Use new big's precomputed integration on fetal SC ===")
new_fsc2 = new[stem_mask & fetal_mask, :].copy()
sc.pp.neighbors(new_fsc2, use_rep="X_scVI", n_neighbors=15, key_added="scvi_nb")
sc.tl.umap(new_fsc2, neighbors_key="scvi_nb", key_added="umap_scvi")
sc.tl.leiden(new_fsc2, resolution=0.3, neighbors_key="scvi_nb", key_added="leiden_scvi_03")
sc.tl.leiden(new_fsc2, resolution=0.5, neighbors_key="scvi_nb", key_added="leiden_scvi_05")
print("leiden_scvi_03:", new_fsc2.obs["leiden_scvi_03"].value_counts().to_string())
print("leiden_scvi_05:", new_fsc2.obs["leiden_scvi_05"].value_counts().to_string())

new_fsc2.obsm["X_umap"] = new_fsc2.obsm["umap_scvi"]
fig = sc.pl.umap(new_fsc2, color=["leiden_scvi_03","leiden_scvi_05","Study_name","batch","age_group","library_preparation_protocol","gut_region","Cell_cycle_phase"],
                 ncols=3, save="_newbig_fetalSC_scvi.png", show=False, return_fig=True)

# Save
new_fsc.write_h5ad(OUT / "newbig_fetalSC_recluster.h5ad")
new_fsc2.write_h5ad(OUT / "newbig_fetalSC_scvi.h5ad")

print("\n=== Done. Outputs in:", OUT)
