"""Q3 quick figure: QC + marker score histograms for newly-recovered Elmentaite stem cells."""
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import re
import matplotlib.pyplot as plt

NEW_BIG = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
OLD_FSC = "/Users/am336941/PhD/data/gut_data/before_remapping/Fetal_healthy_stem_cells_leiden.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
FIGS = OUT / "figs"; FIGS.mkdir(exist_ok=True)

def bc16(s):
    m = re.match(r"^([ACGT]{16})", s)
    return m.group(1) if m else None

print("loading old ...")
old = sc.read_h5ad(OLD_FSC, backed="r")
old.obs["bc16"] = pd.Series(old.obs.index, index=old.obs.index).map(bc16)
old.obs["match_key"] = old.obs["bc16"].astype(str) + "|" + old.obs["Donor_ID"].astype(str)
old_keys = set(old.obs["match_key"])
old.file.close()

print("loading new big ...")
big = sc.read_h5ad(NEW_BIG)
big.obs["bc16"] = pd.Series(big.obs.index, index=big.obs.index).map(bc16)
big.obs["match_key"] = big.obs["bc16"].astype(str) + "|" + big.obs["donor_id"].astype(str)
big.obs["was_old_fetal_sc"] = big.obs["match_key"].isin(old_keys)

# Subset to Elmentaite fetal stem cells
mask = (big.obs["cellstates_scANVI"]=="Stem cells") & big.obs["age_group"].isin(["first trimester","second trimester"]) & (big.obs["Study_name"]=="Elementaite_2021")
elm_fsc = big[mask].copy()
print(f"Elmentaite fetal stem cells: {elm_fsc.n_obs}")
print(f"  in old: {elm_fsc.obs['was_old_fetal_sc'].sum()}")
print(f"  new only: {(~elm_fsc.obs['was_old_fetal_sc']).sum()}")

# Score markers
sc.pp.normalize_total(elm_fsc, target_sum=1e4)
sc.pp.log1p(elm_fsc)
markers = {
    "stem":          ["LGR5","OLFM4","ASCL2","SMOC2","RNF43"],
    "TA":            ["MKI67","PCNA","TOP2A","CDK1","STMN1","UBE2C"],
    "enterocyte":    ["KRT20","ALPI","SI","SLC15A1","SLC5A1","MEP1A","APOA4","FABP1","FABP2"],
    "proximal_prog": ["SOX9","CDX2","HES1","HOPX","BMP4","ID3"],
}
for n, gs in markers.items():
    in_var = [g for g in gs if g in elm_fsc.var.index]
    if in_var:
        sc.tl.score_genes(elm_fsc, gene_list=in_var, score_name=f"score_{n}", use_raw=False, random_state=0)

old_match = elm_fsc[elm_fsc.obs["was_old_fetal_sc"]]
new_only  = elm_fsc[~elm_fsc.obs["was_old_fetal_sc"]]

# Plot QC + scores
fig, axs = plt.subplots(2, 4, figsize=(22, 10))
qc_metrics = ["n_genes_by_counts","total_counts","percent_mito","doublet_scores"]
score_metrics = ["score_stem","score_TA","score_enterocyte","score_proximal_prog"]
for i, m in enumerate(qc_metrics):
    if m not in elm_fsc.obs.columns: continue
    ret = old_match.obs[m].astype(float).values
    new = new_only.obs[m].astype(float).values
    bins = np.linspace(min(ret.min(), new.min()), max(np.quantile(ret, .99), np.quantile(new, .99)), 40)
    axs[0,i].hist(ret, bins=bins, alpha=0.55, label=f"in OLD (n={len(ret)})", density=True)
    axs[0,i].hist(new, bins=bins, alpha=0.55, label=f"new-only (n={len(new)})", density=True)
    axs[0,i].axvline(np.median(ret), color="C0", linestyle="--")
    axs[0,i].axvline(np.median(new), color="C1", linestyle="--")
    axs[0,i].set_xlabel(m); axs[0,i].set_title(m); axs[0,i].legend()
for i, m in enumerate(score_metrics):
    if m not in elm_fsc.obs.columns: continue
    ret = old_match.obs[m].astype(float).values
    new = new_only.obs[m].astype(float).values
    axs[1,i].boxplot([ret, new], tick_labels=[f"in OLD\n(n={len(ret)})", f"new-only\n(n={len(new)})"], showmeans=True)
    axs[1,i].set_ylabel(m); axs[1,i].set_title(m)
fig.suptitle("Q3: newly-recovered Elmentaite stem cells — QC and marker scores", fontsize=14)
fig.savefig(FIGS / "10_Q3_newly_recovered_QC.png", dpi=120, bbox_inches="tight")
plt.close(fig)
print("saved Q3 figure.")
