"""Inspect metadata of both datasets in backed mode."""
import anndata as ad
import pandas as pd
import json
from pathlib import Path

OLD_PATH = "/Users/am336941/PhD/data/gut_data/Healthy_colon_adult.h5ad"
NEW_PATH = "/Users/am336941/PhD/data/gut_data/gut_hs_all_datasets_full_annotated_AM_30102025_181544_raw.h5ad"
SUBSET_PATH = "/Users/am336941/PhD/data/gut_data/gut_hs_fetalSC_AM_05032025_150941_raw.h5ad"
OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")
OUT.mkdir(exist_ok=True)


def inspect(path, label):
    print(f"\n=== {label}: {path} ===")
    a = ad.read_h5ad(path, backed="r")
    print(f"shape: {a.shape}")
    print(f"obs columns ({len(a.obs.columns)}): {list(a.obs.columns)}")
    print(f"var columns ({len(a.var.columns)}): {list(a.var.columns)}")
    print(f"obsm keys: {list(a.obsm.keys())}")
    print(f"varm keys: {list(a.varm.keys())}")
    print(f"layers: {list(a.layers.keys())}")
    print(f"uns keys: {list(a.uns.keys())}")
    print(f"X dtype: {a.X.dtype}")
    print("\n--- obs head ---")
    print(a.obs.head().to_string())
    print("\n--- var head ---")
    print(a.var.head().to_string())
    # Save full obs to disk for downstream
    safe = label.lower().replace(" ", "_")
    a.obs.to_parquet(OUT / f"obs_{safe}.parquet")
    a.var.to_parquet(OUT / f"var_{safe}.parquet")
    # Per-column unique-value summary for low-cardinality columns
    summary = {}
    for c in a.obs.columns:
        try:
            n_uniq = a.obs[c].nunique(dropna=False)
            summary[c] = {"n_unique": int(n_uniq), "dtype": str(a.obs[c].dtype)}
            if n_uniq <= 60:
                vc = a.obs[c].value_counts(dropna=False).head(60)
                summary[c]["top_values"] = vc.to_dict()
        except Exception as e:
            summary[c] = {"error": str(e)}
    with open(OUT / f"obs_summary_{safe}.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    a.file.close()
    return summary


s_old = inspect(OLD_PATH, "old")
s_new = inspect(NEW_PATH, "new")
s_sub = inspect(SUBSET_PATH, "subset")

print("\n=== Done. Outputs in:", OUT)
