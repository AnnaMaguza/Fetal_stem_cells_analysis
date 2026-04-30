"""Inspect the actual OLD fetal stem-cell files (before remapping)."""
import anndata as ad
import pandas as pd
import json
from pathlib import Path

OLD_DIR = "/Users/am336941/PhD/data/gut_data/before_remapping"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")

files = [
    ("old_fsc_leiden",  f"{OLD_DIR}/Fetal_healthy_stem_cells_leiden.h5ad"),
    ("old_fsc_scvi",    f"{OLD_DIR}/Fetal_healthy_stem_cells_scvi.h5ad"),
    ("old_fsc",         f"{OLD_DIR}/Fetal_healthy_stem_cells.h5ad"),
    ("old_integrated",  f"{OLD_DIR}/Integrated_4_datasets_05042024.h5ad"),
]

for label, path in files:
    print(f"\n========== {label}: {path} ==========")
    a = ad.read_h5ad(path, backed="r")
    print(f"shape: {a.shape}")
    print(f"X dtype: {a.X.dtype}")
    print(f"obsm: {list(a.obsm.keys())}")
    print(f"layers: {list(a.layers.keys())}")
    print(f"obs columns ({len(a.obs.columns)}): {list(a.obs.columns)}")
    print(f"var columns ({len(a.var.columns)}): {list(a.var.columns)}")
    print(f"obs head:\n{a.obs.head(3).to_string()}")
    print(f"\nvar index head: {list(a.var.index[:8])}")
    if "gene_id" in a.var.columns:
        print(f"var gene_id head: {list(a.var['gene_id'][:5])}")
    a.obs.to_parquet(OUT / f"obs_{label}.parquet")
    a.var.to_parquet(OUT / f"var_{label}.parquet")

    # Look for clustering / leiden columns
    for col in ["leiden","Leiden","clusters","cluster","cell_states","cell_state","Cell_states","Cell States","celltype"]:
        if col in a.obs.columns:
            vc = a.obs[col].value_counts(dropna=False).head(20)
            print(f"\n  '{col}' ↦ top values:")
            print(vc.to_string())
    a.file.close()

print("\n=== Done ===")
