"""Q5: re-cluster fetal SC without enteroid cells (Holloway primary kept if any)."""
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

SUBSET = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
NEW_BIG = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

print("Loading subset ...")
sub = sc.read_h5ad(SUBSET)
print(f"subset: {sub.shape}")

# Examine Holloway sample structure to find non-enteroid cells
hol = sub[sub.obs["Study_name"] == "Holloway_2021"]
print(f"\nHolloway 2021 cells in subset: {hol.n_obs}")
print("\ngrowth_condition:")
print(hol.obs["growth_condition"].value_counts(dropna=False).to_string())
print("\nMaterial Type:")
print(hol.obs["Material Type"].value_counts(dropna=False).to_string())
print("\ndevelopmental_stage:")
print(hol.obs["developmental_stage"].value_counts(dropna=False).to_string())
print("\nfull_age:")
print(hol.obs["full_age"].value_counts(dropna=False).to_string())
print("\nLibrary_preparation_protocol:")
print(hol.obs["library_preparation_protocol"].value_counts(dropna=False).to_string())
print("\nSample_id top 10:")
print(hol.obs["sample_id"].value_counts(dropna=False).head(10).to_string())

# Identify enteroid cells: anything with 'enteroid' in growth_condition is enteroid
ent_mask = sub.obs["growth_condition"].astype(str).str.contains("enteroid", case=False, na=False)
print(f"\nEnteroid cells (by growth_condition): {ent_mask.sum()}")
print(f"  by Study: {sub.obs.loc[ent_mask,'Study_name'].value_counts().to_dict()}")
print(f"Non-enteroid cells: {(~ent_mask).sum()}")
print(f"  by Study: {sub.obs.loc[~ent_mask,'Study_name'].value_counts().to_dict()}")

# Are there ANY Holloway cells that are NOT enteroid?
hol_non_ent = hol.obs.loc[~hol.obs["growth_condition"].astype(str).str.contains("enteroid", case=False, na=False)]
print(f"\nHolloway 2021 cells NOT labeled enteroid: {len(hol_non_ent)}")
if len(hol_non_ent) > 0:
    print(hol_non_ent["growth_condition"].value_counts().to_string())

# The new big atlas has more Holloway cells including primary if any
print("\n--- Checking Holloway in NEW BIG atlas ---")
import anndata as ad
big = ad.read_h5ad(NEW_BIG, backed="r")
hol_big_mask = big.obs["Study_name"] == "Holloway_2021"
print(f"Holloway in big: {hol_big_mask.sum()}")
hol_big_obs = big.obs[hol_big_mask]
print("Holloway big growth_condition:")
print(hol_big_obs["growth_condition"].value_counts(dropna=False).head(20).to_string())
print("Holloway big age_group:")
print(hol_big_obs["age_group"].value_counts(dropna=False).to_string())
print("Holloway big sample_id top 10:")
print(hol_big_obs["sample_id"].value_counts(dropna=False).head(10).to_string())
big.file.close()

# Re-cluster: keep only primary tissue (non-enteroid) fetal SC
sub_no_ent = sub[~ent_mask].copy()
print(f"\n=== Reclustering: {sub_no_ent.n_obs} primary-tissue fetal SC ===")
print(f"  Study_name: {sub_no_ent.obs['Study_name'].value_counts().to_dict()}")
print(f"  age_group:  {sub_no_ent.obs['age_group'].value_counts().to_dict()}")
print(f"  donor_id (top 10): {sub_no_ent.obs['donor_id'].value_counts().head(10).to_dict()}")

# Re-do clustering
sub_no_ent.X = sub_no_ent.X.astype(np.float32)
# Save raw counts in layer
sub_no_ent.layers["counts"] = sub_no_ent.X.copy()
sc.pp.normalize_total(sub_no_ent, target_sum=1e4)
sc.pp.log1p(sub_no_ent)
sc.pp.highly_variable_genes(sub_no_ent, n_top_genes=3000, flavor="seurat", batch_key="batch")
sub_no_ent_h = sub_no_ent[:, sub_no_ent.var["highly_variable"]].copy()
sc.pp.scale(sub_no_ent_h, max_value=10)
sc.tl.pca(sub_no_ent_h, n_comps=30)
sc.pp.neighbors(sub_no_ent_h, n_neighbors=15, n_pcs=30)
sc.tl.umap(sub_no_ent_h)
sc.tl.leiden(sub_no_ent_h, resolution=0.3, key_added="leiden_03")
sc.tl.leiden(sub_no_ent_h, resolution=0.5, key_added="leiden_05")
sc.tl.leiden(sub_no_ent_h, resolution=0.8, key_added="leiden_08")
print(f"\nFresh PCA Leiden 0.3: {sub_no_ent_h.obs['leiden_03'].value_counts().to_dict()}")
print(f"Fresh PCA Leiden 0.5: {sub_no_ent_h.obs['leiden_05'].value_counts().to_dict()}")

# Use the focused scVI from subset for a second view
if "X_scVI" in sub_no_ent.obsm:
    sc.pp.neighbors(sub_no_ent, use_rep="X_scVI", n_neighbors=15, key_added="scvi_nb")
    sc.tl.umap(sub_no_ent, neighbors_key="scvi_nb", key_added="umap_scvi")
    sc.tl.leiden(sub_no_ent, resolution=0.3, neighbors_key="scvi_nb", key_added="leiden_scvi_03")
    sc.tl.leiden(sub_no_ent, resolution=0.5, neighbors_key="scvi_nb", key_added="leiden_scvi_05")
    print(f"\nfocused scVI Leiden 0.3: {sub_no_ent.obs['leiden_scvi_03'].value_counts().to_dict()}")
    print(f"focused scVI Leiden 0.5: {sub_no_ent.obs['leiden_scvi_05'].value_counts().to_dict()}")

# Plot fresh PCA
fig, axs = plt.subplots(2, 3, figsize=(20, 12))
sc.pl.umap(sub_no_ent_h, color="leiden_03", ax=axs[0,0], show=False, title=f"No enteroids (n={sub_no_ent.n_obs}) — fresh PCA — leiden 0.3", legend_loc="on data")
sc.pl.umap(sub_no_ent_h, color="leiden_05", ax=axs[0,1], show=False, title="leiden 0.5", legend_loc="on data")
sc.pl.umap(sub_no_ent_h, color="leiden_08", ax=axs[0,2], show=False, title="leiden 0.8", legend_loc="on data")
sc.pl.umap(sub_no_ent_h, color="Cell_cycle_phase", ax=axs[1,0], show=False, title="cell cycle")
sc.pl.umap(sub_no_ent_h, color="donor_id", ax=axs[1,1], show=False, title="donor", legend_loc="none")
sc.pl.umap(sub_no_ent_h, color="gut_region", ax=axs[1,2], show=False, title="gut region")
fig.suptitle(f"Q5: fetal SC without enteroids — fresh PCA reclustering (n={sub_no_ent.n_obs})", fontsize=14)
fig.savefig(FIGS / "12_Q5_no_enteroids_freshPCA.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Plot focused scVI
if "leiden_scvi_03" in sub_no_ent.obs.columns:
    sub_no_ent.obsm["X_umap"] = sub_no_ent.obsm["umap_scvi"]
    fig, axs = plt.subplots(2, 3, figsize=(20, 12))
    sc.pl.umap(sub_no_ent, color="leiden_scvi_03", ax=axs[0,0], show=False, title="No enteroids — focused scVI — leiden 0.3", legend_loc="on data")
    sc.pl.umap(sub_no_ent, color="leiden_scvi_05", ax=axs[0,1], show=False, title="leiden 0.5", legend_loc="on data")
    sc.pl.umap(sub_no_ent, color="Cell_cycle_phase", ax=axs[0,2], show=False, title="cell cycle")
    sc.pl.umap(sub_no_ent, color="donor_id", ax=axs[1,0], show=False, title="donor", legend_loc="none")
    sc.pl.umap(sub_no_ent, color="gut_region", ax=axs[1,1], show=False, title="gut region")
    sc.pl.umap(sub_no_ent, color="age_group", ax=axs[1,2], show=False, title="age group")
    fig.suptitle(f"Q5: fetal SC without enteroids — focused scVI reclustering (n={sub_no_ent.n_obs})", fontsize=14)
    fig.savefig(FIGS / "13_Q5_no_enteroids_scvi.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

# Plot key markers
fig, axs = plt.subplots(3, 4, figsize=(22, 15))
markers_to_plot = ["LGR5","OLFM4","ASCL2","SMOC2","MKI67","PCNA","TOP2A","SOX9","KRT20","ALPI","MUC2","CHGA"]
for i, m in enumerate(markers_to_plot):
    if m in sub_no_ent_h.var.index:
        sc.pl.umap(sub_no_ent_h, color=m, ax=axs[i//4, i%4], show=False, title=m, use_raw=False)
fig.suptitle("Markers on no-enteroid reclustering (fresh PCA)")
fig.savefig(FIGS / "14_Q5_markers_freshPCA.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Save (clean obs first to avoid h5ad issue)
for c in list(sub_no_ent.obs.columns):
    if sub_no_ent.obs[c].dtype == 'O':
        sub_no_ent.obs[c] = sub_no_ent.obs[c].astype(str)
for c in list(sub_no_ent_h.obs.columns):
    if sub_no_ent_h.obs[c].dtype == 'O':
        sub_no_ent_h.obs[c] = sub_no_ent_h.obs[c].astype(str)
try:
    sub_no_ent.write_h5ad(OUT / "Q5_fetalSC_no_enteroids.h5ad")
    sub_no_ent_h.write_h5ad(OUT / "Q5_fetalSC_no_enteroids_freshPCA.h5ad")
    print(f"\nSaved Q5 outputs.")
except Exception as e:
    print(f"\nh5ad save warning: {e}")

# Top DE genes between leiden 0.3 clusters (no enteroid)
print("\n=== DE genes per leiden 0.3 cluster (fresh PCA) ===")
sc.tl.rank_genes_groups(sub_no_ent_h, "leiden_03", method="wilcoxon", n_genes=15)
for cl in sorted(sub_no_ent_h.obs["leiden_03"].unique()):
    print(f"\nCluster {cl} (n={(sub_no_ent_h.obs['leiden_03']==cl).sum()}) top markers:")
    df = sc.get.rank_genes_groups_df(sub_no_ent_h, group=cl).head(15)
    print(df.to_string(index=False))

print("\n=== Done ===")
