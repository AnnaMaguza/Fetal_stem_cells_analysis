"""Final visual comparison: subset (3-cluster) vs new big (re-clustered)."""
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

SUBSET_PATH = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"

# ---------- 1. Subset: 3-cluster reference (Mar-2025 analysis) ----------
sub = sc.read_h5ad(SUBSET_PATH)
print(f"subset shape: {sub.shape}")
print("subset obsm keys:", list(sub.obsm.keys()))
print("subset leiden_cluster:")
print(sub.obs["leiden_cluster"].value_counts(dropna=False))

fig, axs = plt.subplots(2, 3, figsize=(20, 12))
sc.pl.umap(sub, color="leiden_cluster", ax=axs[0,0], show=False, title="Subset: leiden_cluster (Mar 2025)", legend_loc="on data")
sc.pl.umap(sub, color="Study_name", ax=axs[0,1], show=False, title="Subset: Study")
sc.pl.umap(sub, color="age_group", ax=axs[0,2], show=False, title="Subset: age_group")
sc.pl.umap(sub, color="donor_id", ax=axs[1,0], show=False, title="Subset: donor", legend_loc="none")
sc.pl.umap(sub, color="library_preparation_protocol", ax=axs[1,1], show=False, title="Subset: lib protocol")
sc.pl.umap(sub, color="Cell_cycle_phase", ax=axs[1,2], show=False, title="Subset: cell cycle")
fig.suptitle("Fetal SC subset (Mar-2025): 3 clean clusters", fontsize=14)
fig.savefig(FIGS / "01_subset_fetalSC_3clusters.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Show key markers in subset
fig, axs = plt.subplots(2, 4, figsize=(22, 10))
for i, m in enumerate(["LGR5","OLFM4","ASCL2","SMOC2","MKI67","PCNA","TOP2A","SOX9"]):
    if m in sub.var.index:
        sc.pl.umap(sub, color=m, ax=axs[i//4, i%4], show=False, title=m, layer="counts" if "counts" in sub.layers else None)
fig.suptitle("Subset fetal SC: stem-cell + cycling markers")
fig.savefig(FIGS / "02_subset_markers.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Composition by cluster
print("\n=== Subset: cluster composition ===")
ct = pd.crosstab(sub.obs["leiden_cluster"], sub.obs["Cell_cycle_phase"])
print(ct.to_string())
ct2 = pd.crosstab(sub.obs["leiden_cluster"], sub.obs["Study_name"])
print(ct2.to_string())
ct3 = pd.crosstab(sub.obs["leiden_cluster"], sub.obs["gut_region"])
print(ct3.to_string())
ct.to_csv(OUT / "subset_cluster_x_cellcycle.csv")
ct2.to_csv(OUT / "subset_cluster_x_study.csv")
ct3.to_csv(OUT / "subset_cluster_x_region.csv")

# ---------- 2. New big fetal SC re-clustered (scVI from full atlas) ----------
new_fsc2 = ad.read_h5ad(OUT / "newbig_fetalSC_scvi.h5ad")
new_fsc2.obsm["X_umap"] = new_fsc2.obsm["umap_scvi"]

fig, axs = plt.subplots(2, 3, figsize=(20, 12))
sc.pl.umap(new_fsc2, color="leiden_scvi_03", ax=axs[0,0], show=False, title=f"NewBig fetal SC: leiden 0.3 (scVI from full atlas)", legend_loc="none")
sc.pl.umap(new_fsc2, color="leiden_scvi_05", ax=axs[0,1], show=False, title="NewBig fetal SC: leiden 0.5", legend_loc="none")
sc.pl.umap(new_fsc2, color="Study_name", ax=axs[0,2], show=False, title="Study")
sc.pl.umap(new_fsc2, color="age_group", ax=axs[1,0], show=False, title="age_group")
sc.pl.umap(new_fsc2, color="donor_id", ax=axs[1,1], show=False, title="donor", legend_loc="none")
sc.pl.umap(new_fsc2, color="Cell_cycle_phase", ax=axs[1,2], show=False, title="cell cycle")
fig.suptitle("NewBig fetal SC (n=7308) — clustered using scVI from FULL 387k-cell atlas", fontsize=14)
fig.savefig(FIGS / "03_newbig_fetalSC_global_scVI.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# ---------- 3. New big fetal SC re-clustered fresh (HVG/PCA, no scVI) ----------
new_fsc1 = ad.read_h5ad(OUT / "newbig_fetalSC_recluster.h5ad")
fig, axs = plt.subplots(2, 3, figsize=(20, 12))
sc.pl.umap(new_fsc1, color="leiden_03", ax=axs[0,0], show=False, title="NewBig fetal SC: leiden 0.3 (fresh PCA)", legend_loc="none")
sc.pl.umap(new_fsc1, color="leiden_05", ax=axs[0,1], show=False, title="leiden 0.5", legend_loc="none")
sc.pl.umap(new_fsc1, color="Study_name", ax=axs[0,2], show=False, title="Study")
sc.pl.umap(new_fsc1, color="batch", ax=axs[1,0], show=False, title="batch", legend_loc="none")
sc.pl.umap(new_fsc1, color="library_preparation_protocol", ax=axs[1,1], show=False, title="lib protocol")
sc.pl.umap(new_fsc1, color="donor_id", ax=axs[1,2], show=False, title="donor", legend_loc="none")
fig.suptitle("NewBig fetal SC — re-clustered FRESH (HVG+PCA, no integration)", fontsize=14)
fig.savefig(FIGS / "04_newbig_fetalSC_fresh_PCA.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# ---------- 4. Side-by-side: subset 3-cluster labels carried over to new big embedding ----------
# Map subset's leiden_cluster onto new big fetal SC and visualize on both UMAPs
subset_leiden = sub.obs["leiden_cluster"].astype(str)
new_fsc1.obs["subset_leiden"] = subset_leiden.reindex(new_fsc1.obs.index).values
new_fsc2.obs["subset_leiden"] = subset_leiden.reindex(new_fsc2.obs.index).values

fig, axs = plt.subplots(1, 3, figsize=(22, 7))
sc.pl.umap(sub, color="leiden_cluster", ax=axs[0], show=False, title="Subset UMAP (focused scVI)\nlabel = subset leiden", legend_loc="on data")
sc.pl.umap(new_fsc2, color="subset_leiden", ax=axs[1], show=False, title="NewBig UMAP (global scVI)\nlabel = subset leiden")
sc.pl.umap(new_fsc1, color="subset_leiden", ax=axs[2], show=False, title="NewBig UMAP (fresh PCA)\nlabel = subset leiden")
fig.suptitle("SAME 7308 cells, SAME labels, 3 different embeddings", fontsize=14)
fig.savefig(FIGS / "05_three_embeddings_compared.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Are the subset's 3 clusters still visible in the new embedding?
ct = pd.crosstab(new_fsc2.obs["subset_leiden"], new_fsc2.obs["leiden_scvi_05"])
ct.to_csv(OUT / "subset_leiden_x_newbig_leiden.csv")
print("\n=== Subset leiden ↔ NewBig leiden_scvi_05 contingency ===")
print(ct.to_string())

print("\nDone — figures in", FIGS)
