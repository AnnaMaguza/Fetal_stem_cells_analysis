"""Match old (pre-remap) fetal SCs to new (post-remap) by (donor + barcode).
Then compare expression on matched cells across the two genome references.
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path
import re
import matplotlib.pyplot as plt

OLD_FSC = "/Users/am336941/PhD/data/gut_data/before_remapping/Fetal_healthy_stem_cells_leiden.h5ad"
SUBSET = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
NEW_BIG = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

# Load old fetal SC (full, in memory — only 158 MB)
print("Loading old fetal SC ...")
old = sc.read_h5ad(OLD_FSC)
print(f"old: {old.shape}, X dtype={old.X.dtype}")
print("old leiden:")
print(old.obs["leiden"].value_counts().to_string())
print("old cluster:")
print(old.obs["cluster"].value_counts().to_string())
print("old Donor_ID head:", old.obs["Donor_ID"].value_counts().head(20).to_string())
print("old Study_name:", old.obs["Study_name"].value_counts().to_string())

# Extract barcode (first 16-bp string) and donor
def bc16(s):
    m = re.match(r"^([ACGT]{16})", s)
    return m.group(1) if m else None
old_obs = old.obs.copy()
old_obs["bc16"] = pd.Series(old.obs.index, index=old.obs.index).map(bc16)
old_obs["match_key"] = old_obs["bc16"].astype(str) + "|" + old_obs["Donor_ID"].astype(str)

# Load new subset (in memory, 7,308 cells)
print("\nLoading new fetal SC subset ...")
new = sc.read_h5ad(SUBSET)
print(f"new: {new.shape}, X dtype={new.X.dtype}")
print("new Study_name:", new.obs["Study_name"].value_counts().to_string())
print("new donor_id top:", new.obs["donor_id"].value_counts().head(20).to_string())

new_obs = new.obs.copy()
new_obs["bc16"] = pd.Series(new.obs.index, index=new.obs.index).map(bc16)
new_obs["match_key"] = new_obs["bc16"].astype(str) + "|" + new_obs["donor_id"].astype(str)

print(f"\nold match_keys unique: {old_obs['match_key'].nunique()} (cells {len(old_obs)})")
print(f"new match_keys unique: {new_obs['match_key'].nunique()} (cells {len(new_obs)})")

common_keys = set(old_obs["match_key"]) & set(new_obs["match_key"])
print(f"\n=== MATCHED cells (barcode + donor) ===")
print(f"matched keys: {len(common_keys)}")

# Take first occurrence per key (in case of dupes)
old_obs2 = old_obs[old_obs["match_key"].isin(common_keys)].drop_duplicates("match_key")
new_obs2 = new_obs[new_obs["match_key"].isin(common_keys)].drop_duplicates("match_key")
print(f"old matched cells: {len(old_obs2)}")
print(f"new matched cells: {len(new_obs2)}")

# Reorder both to same key order
key_order = sorted(common_keys)
old_obs2 = old_obs2.set_index("match_key").reindex(key_order)
new_obs2 = new_obs2.set_index("match_key").reindex(key_order)

# Get matrices
old_X = old[old_obs2["Cell_ID"].values if "Cell_ID" in old_obs2 else old_obs2.index, :]
old_idx = old.obs.reset_index().set_index("match_key" if "match_key" in old.obs.columns else "cell_id")
# Build proper subsets
old_obs_full = old.obs.copy()
old_obs_full["match_key"] = old_obs["match_key"].values
new_obs_full = new.obs.copy()
new_obs_full["match_key"] = new_obs["match_key"].values

# Index lookups
old_to_idx = {k: i for i, k in enumerate(old_obs_full["match_key"]) if k in common_keys}
new_to_idx = {k: i for i, k in enumerate(new_obs_full["match_key"]) if k in common_keys}

old_idx_arr = np.array([old_to_idx[k] for k in key_order])
new_idx_arr = np.array([new_to_idx[k] for k in key_order])

old_sub = old[old_idx_arr].copy()
new_sub = new[new_idx_arr].copy()
print(f"old_sub: {old_sub.shape}, new_sub: {new_sub.shape}")
old_sub.obs["match_key"] = key_order
new_sub.obs["match_key"] = key_order

# Confirm leiden labels are present
print("\nold_sub leiden distribution among matched cells:")
print(old_sub.obs["leiden"].value_counts().to_string())
print("\nold_sub cluster distribution among matched cells:")
print(old_sub.obs["cluster"].value_counts().to_string())

# ----------------------------------------------------------------
# Compare expression for matched cells, on shared gene symbols
# ----------------------------------------------------------------
print("\n=== Expression comparison: same cells, two genome references ===")
print(f"old genes: {old.n_vars} (var index = symbol)")
print(f"new genes: {new.n_vars} (var index mixed symbol/Ensembl)")

old_genes = old.var.index.astype(str)
new_genes = new.var.index.astype(str)
shared_symbols = old_genes.intersection(new_genes)
print(f"shared gene symbols: {len(shared_symbols)}")

old_X_shared = old_sub[:, shared_symbols].X
new_X_shared = new_sub[:, shared_symbols].X
if sp.issparse(old_X_shared): old_X_shared = old_X_shared.toarray()
if sp.issparse(new_X_shared): new_X_shared = new_X_shared.toarray()

# Old X is float32 normalized; new X is int64 raw counts
print(f"old_X_shared dtype: {old_X_shared.dtype}, sample values: {old_X_shared[0,:5]}")
print(f"new_X_shared dtype: {new_X_shared.dtype}, sample values: {new_X_shared[0,:5]}")

# Normalise both to log1p of normalised counts so they're comparable
def normalise(M):
    M = M.astype(np.float32)
    tot = M.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1
    return np.log1p(M / tot * 1e4)

old_norm = normalise(old_X_shared)
new_norm = normalise(new_X_shared)

# Per-cell Pearson correlation
from scipy.stats import pearsonr, spearmanr
n_sample = 500
rng = np.random.RandomState(0)
idx = rng.choice(old_norm.shape[0], min(n_sample, old_norm.shape[0]), replace=False)
pcorrs, scorrs = [], []
for i in idx:
    a = old_norm[i]; b = new_norm[i]
    if a.std() > 0 and b.std() > 0:
        pcorrs.append(pearsonr(a, b)[0])
        scorrs.append(spearmanr(a, b)[0])
pcorrs = np.array(pcorrs); scorrs = np.array(scorrs)

print(f"\nPER-CELL EXPRESSION CORRELATION across genome references (n={len(pcorrs)} sampled cells):")
print(f"  Pearson  mean={pcorrs.mean():.4f}  median={np.median(pcorrs):.4f}  q25={np.quantile(pcorrs,.25):.4f}  q75={np.quantile(pcorrs,.75):.4f}")
print(f"  Spearman mean={scorrs.mean():.4f}  median={np.median(scorrs):.4f}")
print(f"  fraction Pearson > 0.95: {(pcorrs>0.95).mean():.3f}")
print(f"  fraction Pearson > 0.99: {(pcorrs>0.99).mean():.3f}")
print(f"  fraction Pearson > 0.80: {(pcorrs>0.80).mean():.3f}")

# Per-gene comparison: mean expression old vs new on matched cells
print("\n=== Per-gene mean (matched cells) ===")
old_mean = old_norm.mean(axis=0)
new_mean = new_norm.mean(axis=0)
g_pearson = pearsonr(old_mean, new_mean)[0]
print(f"per-gene mean correlation across {len(shared_symbols)} shared genes: r={g_pearson:.4f}")

# Compare per-cell totals
old_total = old_X_shared.sum(axis=1)
new_total = new_X_shared.sum(axis=1)
total_corr = pearsonr(old_total, new_total)[0]
print(f"\nper-cell total counts (raw, on shared genes): old mean={old_total.mean():.1f}, new mean={new_total.mean():.1f}")
print(f"  ratio new/old mean: {new_total.mean()/old_total.mean():.3f}")
print(f"  per-cell total-counts correlation: {total_corr:.4f}")

# Save matched cells with old leiden labels for downstream analyses
matched = pd.DataFrame({
    "match_key": key_order,
    "old_leiden": old_sub.obs["leiden"].values,
    "old_cluster": old_sub.obs["cluster"].values,
    "old_donor": old_sub.obs["Donor_ID"].values,
    "old_age": old_sub.obs["Age"].values,
    "old_region": old_sub.obs["Location"].values,
    "old_total_counts": old_total,
    "new_total_counts": new_total,
    "new_donor": new_sub.obs["donor_id"].values,
    "new_study": new_sub.obs["Study_name"].values,
    "new_age_group": new_sub.obs["age_group"].values,
})
matched.to_csv(OUT / "matched_old_new.csv", index=False)
print(f"\nWrote {OUT/'matched_old_new.csv'} (n={len(matched)})")

# ----------------------------------------------------------------
# Cluster-membership: where did the OLD 3 clusters go in NEW data?
# ----------------------------------------------------------------
# Are ALL the old fetal SCs still labeled "Stem cells" in new big? Need to check NEW BIG for the unmatched ones too.
print("\n=== Tracking unmatched OLD cells in NEW BIG ===")
all_new_big = sc.read_h5ad(NEW_BIG, backed="r")
print(f"new big shape: {all_new_big.shape}")
nbig_obs = all_new_big.obs.copy()
nbig_obs["bc16"] = pd.Series(all_new_big.obs.index, index=all_new_big.obs.index).map(bc16)
nbig_obs["match_key"] = nbig_obs["bc16"].astype(str) + "|" + nbig_obs["donor_id"].astype(str)
all_new_big.file.close()

old_keys = set(old_obs["match_key"])
print(f"\nold fetal SCs total (unique keys): {old_obs['match_key'].nunique()}")
in_newbig = nbig_obs["match_key"].isin(old_keys)
print(f"old fetal SC keys present anywhere in new big: {in_newbig.sum()} cells (over {nbig_obs.loc[in_newbig,'match_key'].nunique()} keys)")

# Their labels in the new big
match_records = nbig_obs[in_newbig].copy()
match_records.to_csv(OUT / "old_fsc_in_newbig.csv")
print("\ncellstates_scANVI distribution for old fetal SCs in new big:")
print(match_records["cellstates_scANVI"].value_counts(dropna=False).head(30).to_string())
print("\ncelltype distribution for old fetal SCs in new big:")
print(match_records["celltype"].value_counts(dropna=False).to_string())
print("\nage_group distribution for old fetal SCs in new big:")
print(match_records["age_group"].value_counts(dropna=False).to_string())
print("\nStudy_name distribution for old fetal SCs in new big:")
print(match_records["Study_name"].value_counts(dropna=False).to_string())

# Plot: subset's UMAP recolored by OLD leiden (where matched)
# but on the OLD UMAP, recolor by where they ended up in new big
# First let's just recover: among matched cells, where do they live in new subset's UMAP?

# Map old leiden onto new subset by match_key
new_obs_with_old = new.obs.copy()
new_obs_with_old["match_key"] = new_obs["match_key"].values
key_to_old_leiden = pd.Series(old_obs["leiden"].values, index=old_obs["match_key"]).to_dict()
key_to_old_cluster = pd.Series(old_obs["cluster"].values, index=old_obs["match_key"]).to_dict()
new.obs["old_leiden"] = new_obs_with_old["match_key"].map(key_to_old_leiden).values
new.obs["old_cluster"] = new_obs_with_old["match_key"].map(key_to_old_cluster).values
print("\n=== old leiden labels propagated to new subset ===")
print(new.obs["old_leiden"].value_counts(dropna=False).to_string())

# Plot
fig, axs = plt.subplots(1, 3, figsize=(22, 7))
sc.pl.umap(old, color="leiden", ax=axs[0], show=False, title="OLD UMAP (pre-remap, 7,817 cells, 3 leiden clusters)", legend_loc="on data")
sc.pl.umap(new, color="old_cluster", ax=axs[1], show=False, title="NEW SUBSET UMAP\nrecoloured by OLD cluster name (matched cells only)")
sc.pl.umap(new, color="leiden_cluster", ax=axs[2], show=False, title="NEW SUBSET UMAP — its own leiden")
fig.suptitle("Old (pre-remap) clusters → tracked into new subset", fontsize=14)
fig.savefig(FIGS / "06_old_clusters_tracked.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Contingency: old leiden × new leiden_cluster (subset)
ct = pd.crosstab(new.obs["old_leiden"], new.obs["leiden_cluster"], dropna=False)
ct.to_csv(OUT / "old_leiden_x_new_leiden.csv")
print("\nContingency old leiden × new leiden_cluster:")
print(ct.to_string())

# Save final old fetal SC obs with updated info
old.obs.to_csv(OUT / "old_fsc_obs_full.csv")
new.obs.to_csv(OUT / "new_fsc_subset_obs_with_old_labels.csv")

print("\n=== Done ===")
