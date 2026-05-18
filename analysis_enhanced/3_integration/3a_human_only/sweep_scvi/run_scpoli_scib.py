"""Compute scIB on every scPoli sweep output (runs in lgr5_enhanced, where
scib-metrics is installed). Outputs:
    sweep_scpoli/<name>/scib_metrics.csv
    sweep_scpoli_summary.csv  (combined)
"""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path
import anndata as ad
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

DATA_OUT = Path("/Users/am336941/PhD/data/Fetal_stem_cells_analysis_enhanced")
SWEEP = DATA_OUT / "sweep_scpoli"


def main() -> int:
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
    bio = BioConservation(
        isolated_labels=True, nmi_ari_cluster_labels_kmeans=True,
        silhouette_label=True, clisi_knn=True,
    )
    batch = BatchCorrection(
        graph_connectivity=True, ilisi_knn=True, kbet_per_label=False,
        pcr_comparison=True, bras=True,
    )
    rows = []
    for run_dir in sorted(p for p in SWEEP.iterdir() if p.is_dir() and (p / "integrated.h5ad").exists()):
        name = run_dir.name
        scib_csv = run_dir / "scib_metrics.csv"
        print(f"[{name}] running scIB on {run_dir/'integrated.h5ad'}")
        a = ad.read_h5ad(run_dir / "integrated.h5ad")
        a = a[a.obs["cell_states"].astype(str) != "unknown"].copy()
        if "X_scpoli" not in a.obsm:
            print(f"[{name}] no X_scpoli in obsm — skip"); continue
        bm = Benchmarker(
            a, batch_key="Study_name", label_key="cell_states",
            embedding_obsm_keys=["X_scpoli"],
            bio_conservation_metrics=bio,
            batch_correction_metrics=batch, n_jobs=2,
        )
        bm.benchmark()
        res = bm.get_results(min_max_scale=False, clean_names=True)
        res.to_csv(scib_csv)
        print(f"[{name}] scIB:\n{res}")
        manifest = json.loads((run_dir / "manifest.json").read_text())
        row = {"name": name, **{f"cfg__{k}": v for k, v in manifest.get("config", {}).items()},
               "wallclock_min": manifest.get("wallclock_min"),
               **{c: res.iloc[0][c] for c in res.columns}}
        rows.append(row)
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(DATA_OUT / "sweep_scpoli_summary.csv", index=False)
        print("wrote", DATA_OUT / "sweep_scpoli_summary.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
