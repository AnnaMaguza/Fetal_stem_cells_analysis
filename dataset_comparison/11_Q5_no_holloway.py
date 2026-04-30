"""Q5 (corrected): re-cluster fetal SC excluding Holloway enteroid cells.
The subset's growth_condition is sparse; use Study_name=='Holloway_2021' as the
enteroid identifier (the Holloway 2021 dataset is fetal-enteroid in its entirety).
"""
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

SUBSET = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
NEW_BIG = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

# Quick check: is there any non-enteroid Holloway in the new big atlas?
print("=== Checking Holloway 2021 across new big atlas ===")
big = ad.read_h5ad(NEW_BIG, backed="r")
hol_obs = big.obs[big.obs["Study_name"] == "Holloway_2021"]
print(f"Holloway 2021 total cells: {len(hol_obs)}")
print(f"\ngrowth_condition value counts (top 30):")
print(hol_obs["growth_condition"].astype(str).value_counts().head(30).to_string())
print(f"\nMaterial Type:")
print(hol_obs["Material Type"].astype(str).value_counts().to_string())
print(f"\nage_group:")
print(hol_obs["age_group"].astype(str).value_counts().to_string())
print(f"\nfull_age (top 15):")
print(hol_obs["full_age"].astype(str).value_counts().head(15).to_string())
big.file.close()

# Load subset
print("\nLoading subset ...")
sub = sc.read_h5ad(SUBSET)
print(f"subset: {sub.shape}")

# Identify Holloway cells (these are the enteroids)
hol_mask = sub.obs["Study_name"] == "Holloway_2021"
print(f"\nHolloway 2021 cells (enteroids) in subset: {hol_mask.sum()}")
print(f"Elementaite 2021 cells (primary tissue) in subset: {(~hol_mask).sum()}")

# Subset to Elmentaite only (true primary fetal tissue)
sub_pt = sub[~hol_mask].copy()
print(f"\n=== Reclustering: {sub_pt.n_obs} primary-tissue fetal SC (Elmentaite only) ===")
print(f"  age_group: {sub_pt.obs['age_group'].value_counts().to_dict()}")
print(f"  donor_id (top 15): {sub_pt.obs['donor_id'].value_counts().head(15).to_dict()}")
print(f"  gut_region: {sub_pt.obs['gut_region'].value_counts().to_dict()}")

# Fresh PCA reclustering
sub_pt.X = sub_pt.X.astype(np.float32)
sub_pt.layers["counts"] = sub_pt.X.copy()
sc.pp.normalize_total(sub_pt, target_sum=1e4)
sc.pp.log1p(sub_pt)
sc.pp.highly_variable_genes(sub_pt, n_top_genes=3000, flavor="seurat", batch_key="batch")
sub_pt_h = sub_pt[:, sub_pt.var["highly_variable"]].copy()
sc.pp.scale(sub_pt_h, max_value=10)
sc.tl.pca(sub_pt_h, n_comps=30)
sc.pp.neighbors(sub_pt_h, n_neighbors=15, n_pcs=30)
sc.tl.umap(sub_pt_h)
sc.tl.leiden(sub_pt_h, resolution=0.3, key_added="leiden_03", flavor="igraph", n_iterations=2, directed=False)
sc.tl.leiden(sub_pt_h, resolution=0.5, key_added="leiden_05", flavor="igraph", n_iterations=2, directed=False)
sc.tl.leiden(sub_pt_h, resolution=0.8, key_added="leiden_08", flavor="igraph", n_iterations=2, directed=False)
print(f"\n=== Fresh PCA Leiden 0.3 ({sub_pt_h.obs['leiden_03'].nunique()} clusters): ===")
print(sub_pt_h.obs["leiden_03"].value_counts().to_string())
print(f"\n=== Leiden 0.5 ({sub_pt_h.obs['leiden_05'].nunique()} clusters): ===")
print(sub_pt_h.obs["leiden_05"].value_counts().to_string())

# Use precomputed scVI from focused subset
print("\n=== Using precomputed focused scVI from subset ===")
sc.pp.neighbors(sub_pt, use_rep="X_scVI", n_neighbors=15, key_added="scvi_nb")
sc.tl.umap(sub_pt, neighbors_key="scvi_nb", key_added="umap_scvi")
sc.tl.leiden(sub_pt, resolution=0.3, neighbors_key="scvi_nb", key_added="leiden_scvi_03", flavor="igraph", n_iterations=2, directed=False)
sc.tl.leiden(sub_pt, resolution=0.5, neighbors_key="scvi_nb", key_added="leiden_scvi_05", flavor="igraph", n_iterations=2, directed=False)
print(f"focused scVI Leiden 0.3 ({sub_pt.obs['leiden_scvi_03'].nunique()} clusters):")
print(sub_pt.obs["leiden_scvi_03"].value_counts().to_string())
print(f"\nfocused scVI Leiden 0.5 ({sub_pt.obs['leiden_scvi_05'].nunique()} clusters):")
print(sub_pt.obs["leiden_scvi_05"].value_counts().to_string())

# Plot fresh PCA
fig, axs = plt.subplots(2, 4, figsize=(24, 12))
sc.pl.umap(sub_pt_h, color="leiden_03", ax=axs[0,0], show=False, title=f"Fresh PCA leiden 0.3", legend_loc="on data")
sc.pl.umap(sub_pt_h, color="leiden_05", ax=axs[0,1], show=False, title="leiden 0.5", legend_loc="on data")
sc.pl.umap(sub_pt_h, color="Cell_cycle_phase", ax=axs[0,2], show=False, title="cell cycle")
sc.pl.umap(sub_pt_h, color="age_group", ax=axs[0,3], show=False, title="age group")
sc.pl.umap(sub_pt_h, color="donor_id", ax=axs[1,0], show=False, title="donor", legend_loc="none")
sc.pl.umap(sub_pt_h, color="gut_region", ax=axs[1,1], show=False, title="gut region")
sc.pl.umap(sub_pt_h, color="library_preparation_protocol", ax=axs[1,2], show=False, title="lib protocol")
sc.pl.umap(sub_pt_h, color="batch", ax=axs[1,3], show=False, title="batch", legend_loc="none")
fig.suptitle(f"Q5: Elmentaite-only fetal SC (n={sub_pt.n_obs}) — fresh PCA reclustering", fontsize=14)
fig.savefig(FIGS / "12_Q5_elmentaite_only_freshPCA.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Plot focused scVI
sub_pt.obsm["X_umap"] = sub_pt.obsm["umap_scvi"]
fig, axs = plt.subplots(2, 4, figsize=(24, 12))
sc.pl.umap(sub_pt, color="leiden_scvi_03", ax=axs[0,0], show=False, title="Focused scVI leiden 0.3", legend_loc="on data")
sc.pl.umap(sub_pt, color="leiden_scvi_05", ax=axs[0,1], show=False, title="leiden 0.5", legend_loc="on data")
sc.pl.umap(sub_pt, color="Cell_cycle_phase", ax=axs[0,2], show=False, title="cell cycle")
sc.pl.umap(sub_pt, color="age_group", ax=axs[0,3], show=False, title="age group")
sc.pl.umap(sub_pt, color="donor_id", ax=axs[1,0], show=False, title="donor", legend_loc="none")
sc.pl.umap(sub_pt, color="gut_region", ax=axs[1,1], show=False, title="gut region")
sc.pl.umap(sub_pt, color="full_age", ax=axs[1,2], show=False, title="full_age", legend_loc="none")
sc.pl.umap(sub_pt, color="batch", ax=axs[1,3], show=False, title="batch", legend_loc="none")
fig.suptitle(f"Q5: Elmentaite-only fetal SC (n={sub_pt.n_obs}) — focused scVI", fontsize=14)
fig.savefig(FIGS / "13_Q5_elmentaite_only_scvi.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Markers on fresh PCA
fig, axs = plt.subplots(3, 4, figsize=(22, 15))
markers = ["LGR5","OLFM4","ASCL2","SMOC2","MKI67","PCNA","TOP2A","SOX9","MTRNR2L12","ASS1","FXYD3","CKB"]
for i, m in enumerate(markers):
    if m in sub_pt_h.var.index:
        sc.pl.umap(sub_pt_h, color=m, ax=axs[i//4, i%4], show=False, title=m, use_raw=False)
fig.suptitle("Q5: Elmentaite-only — markers on fresh PCA UMAP")
fig.savefig(FIGS / "14_Q5_markers_freshPCA.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Markers on focused scVI
fig, axs = plt.subplots(3, 4, figsize=(22, 15))
for i, m in enumerate(markers):
    if m in sub_pt.var.index:
        sc.pl.umap(sub_pt, color=m, ax=axs[i//4, i%4], show=False, title=m, use_raw=False)
fig.suptitle("Q5: Elmentaite-only — markers on focused scVI UMAP")
fig.savefig(FIGS / "15_Q5_markers_scvi.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# DE for the cleanest clustering (focused scVI 0.3, since it is integrated and gives few clusters)
print("\n=== DE genes per leiden_scvi_03 cluster (focused scVI; most relevant) ===")
sub_pt_de = sub_pt.copy()
sc.tl.rank_genes_groups(sub_pt_de, "leiden_scvi_03", method="wilcoxon", n_genes=15)
de_dfs = []
for cl in sorted(sub_pt_de.obs["leiden_scvi_03"].unique(), key=lambda x: int(x)):
    n = (sub_pt_de.obs["leiden_scvi_03"] == cl).sum()
    df = sc.get.rank_genes_groups_df(sub_pt_de, group=cl).head(15)
    df["cluster"] = cl; df["n_cells"] = n
    de_dfs.append(df)
    print(f"\n--- Cluster {cl} (n={n}) top 15 markers ---")
    print(df[["names","scores","pvals_adj"]].to_string(index=False))
pd.concat(de_dfs).to_csv(OUT / "Q5_DE_focused_scvi_03.csv", index=False)

# Save
for c in list(sub_pt.obs.columns):
    if sub_pt.obs[c].dtype == 'O':
        sub_pt.obs[c] = sub_pt.obs[c].astype(str)
for c in list(sub_pt_h.obs.columns):
    if sub_pt_h.obs[c].dtype == 'O':
        sub_pt_h.obs[c] = sub_pt_h.obs[c].astype(str)
sub_pt.write_h5ad(OUT / "Q5_fetalSC_elmentaite_only.h5ad")
sub_pt_h.write_h5ad(OUT / "Q5_fetalSC_elmentaite_only_freshPCA.h5ad")
print("\n=== Done ===")
