"""Compare RAW counts (int) old vs new for matched cells, on shared gene symbols.
Uses Fetal_healthy_stem_cells_scvi.h5ad which has layers['counts'].
Also tracks where the old fetal SCs ended up (label, study) in the new big.
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path
import re
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

OLD_FSC_SCVI = "/Users/am336941/PhD/data/gut_data/before_remapping/Fetal_healthy_stem_cells_scvi.h5ad"
OLD_FSC_LEIDEN = "/Users/am336941/PhD/data/gut_data/before_remapping/Fetal_healthy_stem_cells_leiden.h5ad"
SUBSET = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

# old: load counts layer (raw int) — but X is HVG-restricted. Use the full file's X if it's raw.
print("Loading old fetal SC scvi (has counts layer) ...")
old_scvi = sc.read_h5ad(OLD_FSC_SCVI)
print(f"old_scvi: {old_scvi.shape}, X dtype={old_scvi.X.dtype}, layers: {list(old_scvi.layers.keys())}")
print(f"counts layer dtype: {old_scvi.layers['counts'].dtype}")
print(f"counts layer sample: {old_scvi.layers['counts'].toarray()[0,:10] if sp.issparse(old_scvi.layers['counts']) else old_scvi.layers['counts'][0,:10]}")
# Let's verify counts layer is int counts
cl = old_scvi.layers["counts"]
if sp.issparse(cl):
    print(f"counts layer first cell sum: {cl[0].sum():.1f}")
else:
    print(f"counts layer first cell sum: {cl[0].sum():.1f}")

# Get old labels from leiden file (full feature set)
print("\nLoading old fetal SC leiden (for cluster labels) ...")
old_full = ad.read_h5ad(OLD_FSC_LEIDEN, backed="r")
old_labels = old_full.obs[["leiden", "cluster", "Donor_ID", "Sample_ID", "Cell States", "Cell States Kong", "Age", "Location"]].copy()
old_full.file.close()

# Add labels to scvi
old_scvi.obs = old_scvi.obs.join(old_labels[["leiden", "cluster"]], how="left")
print("old_scvi leiden:", old_scvi.obs["leiden"].value_counts().to_string())

# new
print("\nLoading new fetal SC subset ...")
new = sc.read_h5ad(SUBSET)
print(f"new: {new.shape}, X dtype={new.X.dtype}")

# Build match keys
def bc16(s):
    m = re.match(r"^([ACGT]{16})", s)
    return m.group(1) if m else None
old_scvi.obs["bc16"] = pd.Series(old_scvi.obs.index, index=old_scvi.obs.index).map(bc16)
old_scvi.obs["match_key"] = old_scvi.obs["bc16"].astype(str) + "|" + old_scvi.obs["Donor_ID"].astype(str)
new.obs["bc16"] = pd.Series(new.obs.index, index=new.obs.index).map(bc16)
new.obs["match_key"] = new.obs["bc16"].astype(str) + "|" + new.obs["donor_id"].astype(str)

common_keys = set(old_scvi.obs["match_key"]) & set(new.obs["match_key"])
print(f"matched keys (donor + 16bp barcode): {len(common_keys)}")

# Slice both — first occurrence per key
def first_idx_by_key(adata, keys):
    seen = {}
    out = []
    keys_set = set(keys)
    for i, k in enumerate(adata.obs["match_key"]):
        if k in keys_set and k not in seen:
            seen[k] = i
            out.append(i)
    return seen
old_keymap = first_idx_by_key(old_scvi, common_keys)
new_keymap = first_idx_by_key(new, common_keys)
key_order = sorted(common_keys)
old_idx = [old_keymap[k] for k in key_order]
new_idx = [new_keymap[k] for k in key_order]
old_sub = old_scvi[old_idx].copy()
new_sub = new[new_idx].copy()
print(f"matched: old_sub {old_sub.shape}, new_sub {new_sub.shape}")

# Use the raw counts layer for old, raw X for new (already int)
old_counts = old_sub.layers["counts"]
new_counts = new_sub.X
if sp.issparse(old_counts): old_counts = old_counts.toarray()
if sp.issparse(new_counts): new_counts = new_counts.toarray()

# Old scvi file is HVG-restricted (7000 genes). Get gene names
old_genes = old_sub.var.index.astype(str)  # 7000 HVGs
new_genes = new_sub.var.index.astype(str)
shared_symbols = old_genes.intersection(new_genes)
print(f"\nold genes (HVG only): {len(old_genes)}, new genes: {len(new_genes)}, shared: {len(shared_symbols)}")

# Subset to shared
old_idx_g = [old_sub.var.index.get_loc(g) for g in shared_symbols]
new_idx_g = [new_sub.var.index.get_loc(g) for g in shared_symbols]
old_X_shared = old_counts[:, old_idx_g].astype(np.float32)
new_X_shared = new_counts[:, new_idx_g].astype(np.float32)

print(f"\nold raw counts: dtype={old_X_shared.dtype}, sum/cell mean={old_X_shared.sum(1).mean():.0f}, median={np.median(old_X_shared.sum(1)):.0f}")
print(f"new raw counts: dtype={new_X_shared.dtype}, sum/cell mean={new_X_shared.sum(1).mean():.0f}, median={np.median(new_X_shared.sum(1)):.0f}")
print(f"per-cell counts ratio (new/old): mean={(new_X_shared.sum(1)/np.maximum(old_X_shared.sum(1),1)).mean():.3f}, median={np.median(new_X_shared.sum(1)/np.maximum(old_X_shared.sum(1),1)):.3f}")

# Per-cell n_genes_by_counts
old_ngenes = (old_X_shared > 0).sum(1)
new_ngenes = (new_X_shared > 0).sum(1)
print(f"\nold genes/cell mean (HVGs only): {old_ngenes.mean():.0f}")
print(f"new genes/cell mean (HVGs only, same shared genes): {new_ngenes.mean():.0f}")

# Normalise both to log1p(CP10K)
def lognormcp10k(M):
    tot = M.sum(1, keepdims=True); tot[tot==0] = 1
    return np.log1p(M/tot*1e4)
old_norm = lognormcp10k(old_X_shared)
new_norm = lognormcp10k(new_X_shared)

# Per-cell Pearson correlation on log1p(CP10K) of shared HVGs
n = 1000
rng = np.random.RandomState(0)
idx = rng.choice(old_norm.shape[0], min(n, old_norm.shape[0]), replace=False)
pcorrs = []; rawcorrs = []
for i in idx:
    a, b = old_norm[i], new_norm[i]
    if a.std()>0 and b.std()>0:
        pcorrs.append(pearsonr(a,b)[0])
    a2, b2 = old_X_shared[i], new_X_shared[i]
    if a2.std()>0 and b2.std()>0:
        rawcorrs.append(pearsonr(a2,b2)[0])
pcorrs = np.array(pcorrs); rawcorrs = np.array(rawcorrs)
print(f"\n=== Per-cell expression correlation (RAW counts → log1p(CP10K)), shared HVGs (n={len(pcorrs)} cells sampled) ===")
print(f"  Pearson  mean={pcorrs.mean():.4f}  median={np.median(pcorrs):.4f}  q25={np.quantile(pcorrs,.25):.4f}  q75={np.quantile(pcorrs,.75):.4f}")
print(f"  fraction >0.95: {(pcorrs>0.95).mean():.3f}")
print(f"  fraction >0.99: {(pcorrs>0.99).mean():.3f}")
print(f"\n  RAW Pearson mean={rawcorrs.mean():.4f}  median={np.median(rawcorrs):.4f}")

# Per-gene comparison
old_mean = old_norm.mean(0)
new_mean = new_norm.mean(0)
g_pearson = pearsonr(old_mean, new_mean)[0]
print(f"\nper-gene log1p(CP10K) mean correlation (n={len(shared_symbols)} HVGs): r={g_pearson:.4f}")
g_pcorr_raw = pearsonr(old_X_shared.sum(0), new_X_shared.sum(0))[0]
print(f"per-gene total raw counts correlation: r={g_pcorr_raw:.4f}")

# Identify genes with biggest fold changes
g_log_diff = new_mean - old_mean   # log1p(CP10K) difference
g_df = pd.DataFrame({
    "gene": shared_symbols,
    "old_mean_log1p_CP10K": old_mean,
    "new_mean_log1p_CP10K": new_mean,
    "delta_new_minus_old": g_log_diff,
    "old_total": old_X_shared.sum(0),
    "new_total": new_X_shared.sum(0),
})
g_df["abs_delta"] = g_df["delta_new_minus_old"].abs()
g_df = g_df.sort_values("abs_delta", ascending=False)
g_df.to_csv(OUT / "gene_expression_old_vs_new_HVG.csv", index=False)
print("\nTop 20 genes with biggest |log1p(CP10K)| shift across genomes:")
print(g_df.head(20).to_string(index=False))

# Distribution plot
fig, axs = plt.subplots(1, 3, figsize=(18, 5))
axs[0].hist(pcorrs, bins=40)
axs[0].set_xlabel("Per-cell Pearson r (log1p CP10K, shared HVGs)"); axs[0].set_ylabel("# cells")
axs[0].axvline(np.median(pcorrs), color="red", linestyle="--", label=f"median={np.median(pcorrs):.3f}")
axs[0].legend(); axs[0].set_title("Per-cell expression correlation\nold vs new (matched 955 cells)")

axs[1].scatter(old_X_shared.sum(1), new_X_shared.sum(1), s=2, alpha=0.4)
mx = max(old_X_shared.sum(1).max(), new_X_shared.sum(1).max())
axs[1].plot([0,mx],[0,mx], 'k--', lw=0.5)
axs[1].set_xlabel("OLD per-cell total counts (HVGs)"); axs[1].set_ylabel("NEW per-cell total counts (HVGs)")
axs[1].set_title(f"per-cell totals\nratio new/old median={np.median(new_X_shared.sum(1)/np.maximum(old_X_shared.sum(1),1)):.2f}")

axs[2].scatter(old_mean, new_mean, s=2, alpha=0.4)
mx = max(old_mean.max(), new_mean.max())
axs[2].plot([0,mx],[0,mx], 'k--', lw=0.5)
axs[2].set_xlabel("OLD mean log1p(CP10K)"); axs[2].set_ylabel("NEW mean log1p(CP10K)")
axs[2].set_title(f"per-gene means (shared HVGs)\nr={g_pearson:.3f}")
fig.suptitle("Genome remapping impact on EXPRESSION (955 matched cells, shared HVGs only)")
fig.savefig(FIGS / "07_remap_expression_impact.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# ---------- Track ALL old fetal SC into the new big ----------
print("\n=== Tracking old fetal SC into new big atlas ===")
NEW_BIG = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
big = ad.read_h5ad(NEW_BIG, backed="r")
big_obs = big.obs.copy()
big_obs["bc16"] = pd.Series(big.obs.index, index=big.obs.index).map(bc16)
big_obs["match_key"] = big_obs["bc16"].astype(str) + "|" + big_obs["donor_id"].astype(str)
big.file.close()

old_keys_set = set(old_scvi.obs["match_key"])
big_obs["old_match"] = big_obs["match_key"].isin(old_keys_set)
hits = big_obs[big_obs["old_match"]].copy()

# Bring in old leiden labels for the matched cells
old_keys_to_leiden = pd.Series(old_scvi.obs["leiden"].values, index=old_scvi.obs["match_key"]).to_dict()
old_keys_to_cluster = pd.Series(old_scvi.obs["cluster"].values, index=old_scvi.obs["match_key"]).to_dict()
hits["old_leiden"] = hits["match_key"].map(old_keys_to_leiden)
hits["old_cluster"] = hits["match_key"].map(old_keys_to_cluster)

print(f"\nold cells matched anywhere in new big: {len(hits)} (over {hits['match_key'].nunique()} unique keys)")
print(f"  → out of total old fetal SC: {len(old_keys_set)} unique keys")
print(f"  → retention rate: {hits['match_key'].nunique() / len(old_keys_set) * 100:.1f}%")

# Cross-tab: old_cluster × cellstates_scANVI in new big
ct = pd.crosstab(hits["old_cluster"], hits["cellstates_scANVI"], dropna=False)
ct["TOTAL"] = ct.sum(axis=1)
ct.to_csv(OUT / "old_cluster_x_newbig_cellstates.csv")
print("\n=== old cluster × new big cellstates_scANVI ===")
print(ct.to_string())

# Per-cluster: how many cells lost / kept-as-stem / reassigned
old_cluster_counts = old_scvi.obs["cluster"].value_counts()
print(f"\nOriginal cluster sizes: {old_cluster_counts.to_dict()}")
for cl in old_cluster_counts.index:
    orig_n = old_cluster_counts[cl]
    matched = hits[hits["old_cluster"] == cl]
    matched_n = matched["match_key"].nunique()
    still_stem = (matched["cellstates_scANVI"] == "Stem cells").sum()
    print(f"\n  Cluster '{cl}': original {orig_n}")
    print(f"    matched in new big: {matched_n} ({matched_n/orig_n*100:.1f}%)")
    print(f"    still 'Stem cells': {still_stem} ({still_stem/orig_n*100:.1f}%)")
    print(f"    reassigned to:")
    print((matched.groupby("cellstates_scANVI").size().sort_values(ascending=False).head(8)).to_string())

# ---------- Save key tables ----------
hits.to_csv(OUT / "old_fsc_tracked_in_newbig.csv", index=True)
print("\n=== Done ===")
