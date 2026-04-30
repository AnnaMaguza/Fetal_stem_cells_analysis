"""Q1 + Q2 + Q3 + Q5: deep follow-up analyses.
Q1: Were lost cells low quality?
Q2: Were reassigned cells correctly relabeled?
Q3: Were newly-recovered Elmentaite stem cells good quality?
Q5: Recluster fetal SC excluding enteroid cells.
"""
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path
import re
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

OLD_FSC = "/Users/am336941/PhD/data/gut_data/before_remapping/Fetal_healthy_stem_cells_leiden.h5ad"
NEW_BIG = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
SUBSET = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

def bc16(s):
    m = re.match(r"^([ACGT]{16})", s)
    return m.group(1) if m else None

# ------------------------------------------------------------------
# Q1: were lost cells low-quality?
# ------------------------------------------------------------------
print("\n=== Q1: lost cells QC ===")
old = sc.read_h5ad(OLD_FSC)
old.obs["bc16"] = pd.Series(old.obs.index, index=old.obs.index).map(bc16)
old.obs["match_key"] = old.obs["bc16"].astype(str) + "|" + old.obs["Donor_ID"].astype(str)
print(f"old: {old.shape}, columns: {list(old.obs.columns)[:20]}")

# Read new big in backed mode for matching
big = ad.read_h5ad(NEW_BIG, backed="r")
big_obs = big.obs.copy()
big_obs["bc16"] = pd.Series(big.obs.index, index=big.obs.index).map(bc16)
big_obs["match_key"] = big_obs["bc16"].astype(str) + "|" + big_obs["donor_id"].astype(str)

# big match keys for any cell
big_keys = set(big_obs["match_key"])
old_keys = set(old.obs["match_key"])

old.obs["in_new_big"] = old.obs["match_key"].isin(big_keys)
print(f"\nold cells found in new big atlas: {old.obs['in_new_big'].sum()} / {len(old.obs)}")
print(f"old cells LOST: {(~old.obs['in_new_big']).sum()}")

# QC metrics distribution: lost vs retained
qc_cols = ["n_genes_by_counts", "total_counts", "pct_counts_mito", "pct_counts_ribo"]
print("\n--- QC metrics: lost vs retained ---")
qc_summary = []
for col in qc_cols:
    ret = old.obs.loc[old.obs["in_new_big"], col].astype(float)
    lost = old.obs.loc[~old.obs["in_new_big"], col].astype(float)
    print(f"\n{col}:")
    print(f"  retained (n={len(ret)}): mean={ret.mean():.1f}, median={ret.median():.1f}, q25={ret.quantile(.25):.1f}, q75={ret.quantile(.75):.1f}")
    print(f"  lost     (n={len(lost)}): mean={lost.mean():.1f}, median={lost.median():.1f}, q25={lost.quantile(.25):.1f}, q75={lost.quantile(.75):.1f}")
    try:
        u, p = mannwhitneyu(ret, lost, alternative="two-sided")
        print(f"  Mann-Whitney p = {p:.2e}")
    except Exception as e:
        print(f"  MW err: {e}")
    qc_summary.append({"metric": col, "retained_median": ret.median(), "lost_median": lost.median(),
                        "retained_mean": ret.mean(), "lost_mean": lost.mean(),
                        "n_retained": len(ret), "n_lost": len(lost), "mw_p": p})
pd.DataFrame(qc_summary).to_csv(OUT / "Q1_lost_vs_retained_qc.csv", index=False)

# Plot QC distributions
fig, axs = plt.subplots(1, 4, figsize=(22, 5))
for ax, col in zip(axs, qc_cols):
    ret = old.obs.loc[old.obs["in_new_big"], col].astype(float)
    lost = old.obs.loc[~old.obs["in_new_big"], col].astype(float)
    bins = np.linspace(min(ret.min(), lost.min()), max(ret.quantile(.99), lost.quantile(.99)), 40)
    ax.hist(ret, bins=bins, alpha=0.55, label=f"retained (n={len(ret)})", density=True)
    ax.hist(lost, bins=bins, alpha=0.55, label=f"lost (n={len(lost)})", density=True)
    ax.axvline(ret.median(), color="C0", linestyle="--", alpha=0.7)
    ax.axvline(lost.median(), color="C1", linestyle="--", alpha=0.7)
    ax.set_xlabel(col); ax.set_ylabel("density")
    ax.set_title(col); ax.legend()
fig.suptitle("Q1: Are lost cells low quality? (old QC metrics, lost vs retained)", fontsize=13)
fig.savefig(FIGS / "08_Q1_lost_vs_retained_qc.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# By cluster: how lost are different clusters?
print("\nLoss rate per old cluster:")
for cl in old.obs["cluster"].unique():
    sub = old.obs[old.obs["cluster"] == cl]
    lost_n = (~sub["in_new_big"]).sum()
    print(f"  {cl}: {lost_n}/{len(sub)} = {lost_n/len(sub)*100:.1f}% lost")

# By donor
print("\nLoss rate per donor (top 15 by old size):")
donor_loss = old.obs.groupby("Donor_ID").agg(n_total=("match_key", "count"), n_lost=("in_new_big", lambda x: (~x).sum()))
donor_loss["pct_lost"] = donor_loss["n_lost"] / donor_loss["n_total"] * 100
print(donor_loss.sort_values("n_total", ascending=False).head(15).to_string())
donor_loss.to_csv(OUT / "Q1_loss_by_donor.csv")

# ------------------------------------------------------------------
# Q2: were reassigned cells correctly relabeled?
# ------------------------------------------------------------------
print("\n\n=== Q2: reassigned cells — biological correctness ===")

# Load new big into memory using a slim subset
big.file.close()
print("Loading new big atlas (this takes a minute)...")
big = sc.read_h5ad(NEW_BIG)
print(f"loaded: {big.shape}")
big.obs["bc16"] = pd.Series(big.obs.index, index=big.obs.index).map(bc16)
big.obs["match_key"] = big.obs["bc16"].astype(str) + "|" + big.obs["donor_id"].astype(str)
big.obs["was_old_fetal_sc"] = big.obs["match_key"].isin(old_keys)

# Build maps
old_keys_to_leiden = pd.Series(old.obs["leiden"].values, index=old.obs["match_key"]).to_dict()
old_keys_to_cluster = pd.Series(old.obs["cluster"].values, index=old.obs["match_key"]).to_dict()
big.obs["old_leiden"] = big.obs["match_key"].map(old_keys_to_leiden)
big.obs["old_cluster"] = big.obs["match_key"].map(old_keys_to_cluster)

# Marker scoring
markers = {
    "stem":          ["LGR5","OLFM4","ASCL2","SMOC2","RNF43","TROY"],
    "TA":            ["MKI67","PCNA","TOP2A","CDK1","STMN1","UBE2C"],
    "enterocyte":    ["KRT20","ALPI","SI","SLC15A1","SLC5A1","MEP1A","APOA4","FABP1","FABP2","RBP2"],
    "proximal_prog": ["SOX9","CDX2","HES1","HOPX","BMP4","ID3"],
    "goblet":        ["MUC2","TFF3","CLCA1","SPINK4","REG4"],
    "BEST4_epi":     ["BEST4","OTOP2","CA7","SPIB"],
    "EEC":           ["CHGA","CHGB","NEUROD1","INSM1","TPH1"],
    "paneth":        ["DEFA5","DEFA6","LYZ","ITLN2","ANG4"],
    "tuft":          ["DCLK1","POU2F3","TRPM5","LRMP"],
    "mesenchymal":   ["COL1A1","COL3A1","ACTA2","DCN","VIM","PDGFRA"],
    "mesoderm":      ["HAND1","HAND2","ZEB2","TBX18","LUM"],
    "endothelial":   ["PECAM1","CDH5","KDR","CLDN5"],
    "smc":           ["MYH11","ACTG2","TAGLN","CNN1"],
}

# Compute scores using sc.tl.score_genes for each marker set
big_norm = big.copy()
sc.pp.normalize_total(big_norm, target_sum=1e4)
sc.pp.log1p(big_norm)
for name, gs in markers.items():
    in_var = [g for g in gs if g in big_norm.var.index]
    if in_var:
        sc.tl.score_genes(big_norm, gene_list=in_var, score_name=f"score_{name}", use_raw=False, random_state=0)
        print(f"  scored {name}: {len(in_var)}/{len(gs)} markers in var")

# Pull only old fetal SC cells from big
old_in_big_mask = big_norm.obs["was_old_fetal_sc"].values
recover = big_norm[old_in_big_mask].copy()
print(f"\nold fetal SCs found in new big: {recover.n_obs}")

# Their new label
print("\ndistribution of cellstates_scANVI:")
print(recover.obs["cellstates_scANVI"].value_counts(dropna=False).head(20).to_string())

# For each new label, compare scores
print("\n--- For each new cellstates_scANVI label, mean marker scores ---")
score_cols = [f"score_{n}" for n in markers]
agg = recover.obs.groupby("cellstates_scANVI", observed=True)[score_cols].mean()
agg["n_cells"] = recover.obs.groupby("cellstates_scANVI", observed=True).size()
agg = agg.sort_values("n_cells", ascending=False)
print(agg.head(10).to_string())
agg.to_csv(OUT / "Q2_reassigned_marker_scores.csv")

# Test specifically: stem-vs-non-stem marker scores for each reassignment label
print("\n--- Q2 hypothesis test: do reassigned cells really express the new-label markers? ---")
key_assignments = {
    "Stem cells":          "stem",
    "Proximal progenitor": "proximal_prog",
    "TA":                  "TA",
    "Enterocyte":          "enterocyte",
    "Mesoderm":            "mesoderm",
    "Myofibroblast":       "mesenchymal",
    "SMC":                 "smc",
    "Mesothelium":         "mesoderm",
    "BEST4+ epithelial":   "BEST4_epi",
    "Goblet cell":         "goblet",
}
for label, marker_key in key_assignments.items():
    cells = recover.obs[recover.obs["cellstates_scANVI"] == label]
    if len(cells) < 5: continue
    s_stem = cells[f"score_stem"].values
    s_self = cells[f"score_{marker_key}"].values
    print(f"\n  '{label}' (n={len(cells)}):")
    print(f"     stem score:        mean={s_stem.mean():.3f}, median={np.median(s_stem):.3f}")
    print(f"     {marker_key:14s} score: mean={s_self.mean():.3f}, median={np.median(s_self):.3f}")
    print(f"     {marker_key} > stem in {(s_self > s_stem).sum()}/{len(cells)} = {(s_self > s_stem).mean()*100:.1f}% of cells")

# Plot: among "Stem→Enterocyte" reassignments, do they really express enterocyte markers?
fig, axs = plt.subplots(2, 4, figsize=(20, 10))
big_stem = big_norm.obs["cellstates_scANVI"] == "Stem cells"
old_to_other = {
    "Enterocyte": ("score_enterocyte","score_stem"),
    "TA": ("score_TA","score_stem"),
    "Proximal progenitor": ("score_proximal_prog","score_stem"),
    "Goblet cell": ("score_goblet","score_stem"),
}
for i, (lbl, (s_self, s_stem)) in enumerate(old_to_other.items()):
    # all stem cells in new big
    mask_all_stem = big_stem.values
    mask_old_lbl  = old_in_big_mask & (big_norm.obs["cellstates_scANVI"] == lbl).values
    ax = axs[0, i]
    ax.scatter(big_norm.obs.loc[mask_all_stem, s_stem], big_norm.obs.loc[mask_all_stem, s_self],
               s=2, alpha=0.3, label=f"all 'Stem cells' (n={mask_all_stem.sum()})")
    ax.scatter(big_norm.obs.loc[mask_old_lbl, s_stem], big_norm.obs.loc[mask_old_lbl, s_self],
               s=12, alpha=0.6, color="red", label=f"old fetal SC → '{lbl}' (n={mask_old_lbl.sum()})")
    ax.set_xlabel(s_stem); ax.set_ylabel(s_self); ax.set_title(f"'{lbl}'"); ax.legend()
    # boxplot
    ax = axs[1, i]
    cats = []; vals = []; labs = []
    for clname, mask in [("all 'Stem cells'", mask_all_stem),
                         (f"all '{lbl}'", (big_norm.obs['cellstates_scANVI']==lbl).values),
                         (f"old SC → '{lbl}'", mask_old_lbl)]:
        if mask.sum() > 0:
            v = big_norm.obs.loc[mask, s_self].values
            vals.append(v); labs.append(f"{clname}\nn={mask.sum()}")
    ax.boxplot(vals, labels=labs, showmeans=True)
    ax.set_ylabel(s_self); ax.set_title(f"score_{s_self.split('_',1)[1]} distributions")
    ax.tick_params(axis='x', rotation=15, labelsize=8)
fig.suptitle("Q2: were reassigned old fetal SCs really expressing the new-label markers?", fontsize=13)
fig.savefig(FIGS / "09_Q2_reassigned_validation.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------
# Q3: newly-recovered Elmentaite stem cells (in new but not in old's 7,817)
# ------------------------------------------------------------------
print("\n\n=== Q3: newly-recovered Elmentaite stem cells ===")
new_fetal_sc_mask = (big_norm.obs["cellstates_scANVI"] == "Stem cells") & big_norm.obs["age_group"].isin(["first trimester","second trimester"])
elm_fetal_sc = big_norm[new_fetal_sc_mask & (big_norm.obs["Study_name"]=="Elementaite_2021")].copy()
print(f"new dataset Elmentaite fetal stem cells: {elm_fetal_sc.n_obs}")
elm_in_old = elm_fetal_sc.obs["was_old_fetal_sc"]
print(f"  ↳ in OLD fetal SC analysis (matched by donor+barcode): {elm_in_old.sum()}")
print(f"  ↳ NEW (not in old):                                     {(~elm_in_old).sum()}")

# QC metrics in new big
print("\nQC of NEW Elmentaite stem cells (not in old):")
new_only = elm_fetal_sc[~elm_in_old].copy()
old_match = elm_fetal_sc[elm_in_old].copy()
new_qc = ["n_genes_by_counts", "total_counts", "percent_mito", "percent_ribo", "doublet_scores", "predicted_doublets", "cell_passed_qc"]
print("\n  metric | old_in_new (median) | new_only (median)")
for col in new_qc:
    if col in elm_fetal_sc.obs.columns:
        try:
            ret = old_match.obs[col].astype(float)
            new = new_only.obs[col].astype(float)
            print(f"  {col}: {ret.median():.2f}  vs  {new.median():.2f}")
        except Exception:
            ret = old_match.obs[col].astype(str).value_counts()
            new = new_only.obs[col].astype(str).value_counts()
            print(f"  {col}: {dict(ret)} vs {dict(new)}")

# Marker score comparison
print("\nMarker scores (old_in_new vs new_only):")
for s in ["score_stem","score_TA","score_proximal_prog","score_enterocyte","score_goblet","score_BEST4_epi","score_EEC"]:
    if s in elm_fetal_sc.obs.columns:
        print(f"  {s}: old={old_match.obs[s].median():.3f}  new={new_only.obs[s].median():.3f}")

# Possible reason: were new_only cells in the old INTEGRATED dataset (557k cells) but not in fetal SC subset? Let's see
# The old integrated atlas (Integrated_4_datasets_05042024.h5ad) has 557,099 cells. Did the user filter to only Elmentaite fetal stem cells in the old workflow?
# Let me also look at donor distribution
print("\nDonor distribution of new-only Elmentaite stem cells:")
print(new_only.obs["donor_id"].value_counts().head(20).to_string())
print("\nDonor distribution of old-matched Elmentaite stem cells:")
print(old_match.obs["donor_id"].value_counts().head(20).to_string())

# Check if new-only donors exist in old fetal SC at all
old_donors = set(old.obs["Donor_ID"].astype(str))
new_only_donors = set(new_only.obs["donor_id"].astype(str))
print(f"\nold fetal SC donors: {sorted(old_donors)}")
print(f"new_only donors:     {sorted(new_only_donors)}")
print(f"new_only donors in old fetal SC: {len(new_only_donors & old_donors)} of {len(new_only_donors)}")

# Save subsets for downstream
new_only.write_h5ad(OUT / "Q3_newly_recovered_elm_stem_cells.h5ad")
old_match.write_h5ad(OUT / "Q3_old_matched_elm_stem_cells.h5ad")

# Plot QC histograms for Q3
fig, axs = plt.subplots(2, 4, figsize=(22, 10))
qc_to_plot = ["n_genes_by_counts","total_counts","percent_mito","percent_ribo","doublet_scores","S_score","G2M_score"]
for i, col in enumerate(qc_to_plot):
    ax = axs[i//4, i%4]
    if col in elm_fetal_sc.obs.columns:
        ret = old_match.obs[col].astype(float)
        new = new_only.obs[col].astype(float)
        bins = np.linspace(min(ret.min(), new.min()), max(ret.quantile(.99), new.quantile(.99)), 40)
        ax.hist(ret, bins=bins, alpha=0.55, label=f"in old (n={len(ret)})", density=True)
        ax.hist(new, bins=bins, alpha=0.55, label=f"new-only (n={len(new)})", density=True)
        ax.set_xlabel(col); ax.legend(); ax.set_title(col)
fig.suptitle("Q3: Are newly-recovered Elmentaite stem cells good quality?", fontsize=13)
fig.savefig(FIGS / "10_Q3_newly_recovered_QC.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Marker score boxplots for Q3
fig, axs = plt.subplots(1, 4, figsize=(20, 5))
for i, s in enumerate(["score_stem","score_TA","score_proximal_prog","score_enterocyte"]):
    if s not in elm_fetal_sc.obs.columns: continue
    ax = axs[i]
    ax.boxplot([old_match.obs[s].values, new_only.obs[s].values],
               labels=[f"in old\n(n={old_match.n_obs})", f"new only\n(n={new_only.n_obs})"], showmeans=True)
    ax.set_ylabel(s); ax.set_title(s)
fig.suptitle("Q3: marker scores — newly-recovered vs already-known", fontsize=13)
fig.savefig(FIGS / "11_Q3_marker_scores.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------
# Q5: recluster fetal SC excluding enteroid cells
# ------------------------------------------------------------------
print("\n\n=== Q5: recluster excluding fetal enteroids ===")
sub = sc.read_h5ad(SUBSET)
print(f"subset: {sub.shape}")
print("growth_condition value counts:")
print(sub.obs["growth_condition"].value_counts(dropna=False).head(15).to_string())
print("\nMaterial Type value counts:")
print(sub.obs["Material Type"].value_counts(dropna=False).head(15).to_string())
print("\nStudy_name × growth_condition:")
gxc = pd.crosstab(sub.obs["Study_name"], sub.obs["growth_condition"], dropna=False)
print(gxc.to_string())
print("\nStudy_name × Material Type:")
sxm = pd.crosstab(sub.obs["Study_name"], sub.obs["Material Type"], dropna=False)
print(sxm.to_string())

# Identify enteroid cells: anything with "enteroid" in growth_condition, OR Material Type != "primary tissue", OR specifically Holloway 2021 if they're all enteroid-derived in vitro
ent_mask = sub.obs["growth_condition"].astype(str).str.contains("enteroid", case=False, na=False)
print(f"\ncells with 'enteroid' in growth_condition: {ent_mask.sum()}")

# Also check by material type
mat_cell = sub.obs["Material Type"].astype(str) == "cell"
mat_pt = sub.obs["Material Type"].astype(str) == "primary tissue"
print(f"Material Type=='cell':            {mat_cell.sum()}")
print(f"Material Type=='primary tissue': {mat_pt.sum()}")

# Confirm: enteroid cells = Holloway entirely?
print(f"\nenteroid cells by Study: {sub.obs.loc[ent_mask, 'Study_name'].value_counts().to_string()}")
print(f"non-enteroid cells by Study: {sub.obs.loc[~ent_mask, 'Study_name'].value_counts().to_string()}")

# Subset: keep non-enteroid only
sub_no_ent = sub[~ent_mask].copy()
print(f"\nfetal SC excluding enteroids: {sub_no_ent.shape}")

# Recluster from raw counts using a fresh PCA + UMAP and fresh leiden
sub_no_ent.X = sub_no_ent.X.astype(np.float32)
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
print("\nleiden 0.3 distribution:")
print(sub_no_ent_h.obs["leiden_03"].value_counts().to_string())
print("\nleiden 0.5 distribution:")
print(sub_no_ent_h.obs["leiden_05"].value_counts().to_string())

# Use the precomputed scVI from new big atlas (focused subset)
sub_no_ent.obsm["X_scVI_focused"] = sub_no_ent.obsm["X_scVI"]   # this came from the original integration
sc.pp.neighbors(sub_no_ent, use_rep="X_scVI_focused", n_neighbors=15, key_added="scvi_nb")
sc.tl.umap(sub_no_ent, neighbors_key="scvi_nb", key_added="umap_scvi_focused")
sc.tl.leiden(sub_no_ent, resolution=0.3, neighbors_key="scvi_nb", key_added="leiden_scvi_03")
sc.tl.leiden(sub_no_ent, resolution=0.5, neighbors_key="scvi_nb", key_added="leiden_scvi_05")
print("\nleiden_scvi_03:", sub_no_ent.obs["leiden_scvi_03"].value_counts().to_string())
print("leiden_scvi_05:", sub_no_ent.obs["leiden_scvi_05"].value_counts().to_string())

# Plot Q5
sub_no_ent.obsm["X_umap"] = sub_no_ent.obsm["umap_scvi_focused"]
fig, axs = plt.subplots(2, 3, figsize=(20, 12))
sc.pl.umap(sub_no_ent_h, color="leiden_03", ax=axs[0,0], show=False, title="No enteroids — fresh PCA — leiden 0.3", legend_loc="on data")
sc.pl.umap(sub_no_ent_h, color="leiden_05", ax=axs[0,1], show=False, title="No enteroids — fresh PCA — leiden 0.5", legend_loc="on data")
sc.pl.umap(sub_no_ent_h, color="Cell_cycle_phase", ax=axs[0,2], show=False, title="cell cycle")
sc.pl.umap(sub_no_ent, color="leiden_scvi_03", ax=axs[1,0], show=False, title="No enteroids — focused scVI — leiden 0.3", legend_loc="on data")
sc.pl.umap(sub_no_ent, color="leiden_scvi_05", ax=axs[1,1], show=False, title="leiden 0.5", legend_loc="on data")
sc.pl.umap(sub_no_ent, color="donor_id", ax=axs[1,2], show=False, title="donor")
fig.suptitle(f"Q5: fetal SC excluding enteroids (n={sub_no_ent.n_obs})", fontsize=14)
fig.savefig(FIGS / "12_Q5_no_enteroids.png", dpi=120, bbox_inches="tight")
plt.close(fig)

# Save
sub_no_ent.write_h5ad(OUT / "Q5_fetalSC_no_enteroids.h5ad")
sub_no_ent_h.write_h5ad(OUT / "Q5_fetalSC_no_enteroids_freshPCA.h5ad")

# Compare donor distribution before/after exclusion
print("\nBefore (all fetal SC subset):")
print(sub.obs["donor_id"].value_counts().head(10).to_string())
print("\nAfter (no enteroids):")
print(sub_no_ent.obs["donor_id"].value_counts().head(10).to_string())

print("\n=== Done ===")
