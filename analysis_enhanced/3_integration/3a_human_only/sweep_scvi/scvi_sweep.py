"""Per-config scVI sweep runner on Object A.

Usage:
    uv run --project /Users/am336941/uv_envs/pyenv313 \
        python analysis_enhanced/3_integration/3a_human_only/sweep_scvi/scvi_sweep.py <config_name> [--force]
    uv run --project /Users/am336941/uv_envs/pyenv313 \
        python analysis_enhanced/3_integration/3a_human_only/sweep_scvi/scvi_sweep.py --all [--force]
    uv run --project /Users/am336941/uv_envs/pyenv313 \
        python analysis_enhanced/3_integration/3a_human_only/sweep_scvi/scvi_sweep.py --baseline-only [--force]

Reads `scvi_sweep_configs.json`. Before any scVI config runs, an uncorrected
PCA + UMAP baseline (with scIB metrics) is computed on the same input so each
integration can be judged against it. For each named config:
    1. Load object_a_human.h5ad
    2. HVG-select on raw counts (seurat_v3, batch-aware on Study_name)
    3. Set up scvi-tools AnnData (batch_key=Study_name; per-config covariates)
    4. Train SCVI for max_epochs (no early stopping by default)
    5. Compute UMAP, save figures + scIB-metrics, persist model + integrated h5ad

Output layout:
    figures_root/baseline_uncorrected/   PCA + UMAP + scIB on uncorrected data
    figures_root/<name>/                 training_loss + UMAPs + scIB + manifest
    models_root/<name>/integrated.h5ad   full obj w/ X_scvi + UMAP
    models_root/<name>/model/            scvi model

Idempotent: if `integrated.h5ad` and `scib_metrics.csv` exist, skip unless
`--force` is passed.
"""
from __future__ import annotations

import argparse
import gc
import json
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="scvi")

sc.settings.verbosity = 1
sc.settings.dpi = 120
sc.settings.dpi_save = 300
plt.rcParams.update({
    "savefig.bbox": "tight", "savefig.dpi": 300, "figure.dpi": 120,
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
})

ROOT = Path(__file__).resolve().parent
CFG_PATH = ROOT / "scvi_sweep_configs.json"

COLOR_KEYS = ["Study_name", "cell_states", "age_group", "lgr5_status", "gut_region"]
QC_KEYS = ["total_counts", "n_genes_by_counts"]


def device() -> str:
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def load_configs() -> tuple[dict, list[dict]]:
    blob = json.loads(CFG_PATH.read_text())
    return blob["fixed"], blob["configs"]


def X_is_raw(adata) -> bool:
    return np.array_equal(adata.X.sum(axis=0).astype(int), adata.X.sum(axis=0))


def _prep_for_hvg(adata, fixed) -> None:
    """Apply gene filter and QC metric calculation in-place."""
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, log1p=False)
    min_cells = fixed.get("min_cells_per_gene", 10)
    if min_cells:
        n_before = adata.n_vars
        sc.pp.filter_genes(adata, min_cells=min_cells)
        if "counts" in adata.layers and adata.layers["counts"].shape[1] != adata.n_vars:
            adata.layers["counts"] = adata.X.copy()
        print(f"  gene filter min_cells={min_cells}: {n_before} -> {adata.n_vars}")


def _save_umap_panels(adata, fig_dir: Path, prefix: str) -> None:
    """Save the two standard UMAP composites (metadata + QC)."""
    colors = [c for c in COLOR_KEYS if c in adata.obs.columns]
    with plt.rc_context():
        sc.set_figure_params(dpi=300, figsize=(6, 6))
        sc.pl.umap(adata, color=colors, cmap="magma_r", ncols=3, size=2,
                   frameon=False, show=False)
        plt.savefig(fig_dir / f"{prefix}_metadata.png", bbox_inches="tight")
        plt.close()
    with plt.rc_context():
        sc.set_figure_params(dpi=300, figsize=(6, 6))
        sc.pl.umap(adata, color=QC_KEYS, cmap="magma_r", ncols=2, size=2,
                   frameon=False, show=False)
        plt.savefig(fig_dir / f"{prefix}_qc.png", bbox_inches="tight")
        plt.close()


def _scib_benchmark(adata, embedding_keys: list[str], fixed: dict,
                    out_csv: Path, plot_path: Path) -> pd.DataFrame:
    """Run scIB benchmarker on the given obsm embeddings and save outputs."""
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
    bio = BioConservation(
        isolated_labels=True,
        nmi_ari_cluster_labels_kmeans=True,
        silhouette_label=True,
        clisi_knn=True,
    )
    batch = BatchCorrection(
        graph_connectivity=True,
        ilisi_knn=True,
        kbet_per_label=False,
        pcr_comparison=True,
        bras=True,
    )
    a_lab = adata[adata.obs[fixed["label_key"]].astype(str) != fixed["unlabeled_category"]].copy()
    bm = Benchmarker(
        a_lab,
        batch_key=fixed["batch_key"],
        label_key=fixed["label_key"],
        embedding_obsm_keys=embedding_keys,
        bio_conservation_metrics=bio,
        batch_correction_metrics=batch,
        n_jobs=2,
    )
    bm.benchmark()
    res = bm.get_results(min_max_scale=False, clean_names=True)
    res.to_csv(out_csv)
    plt.figure(figsize=(10, 6))
    bm.plot_results_table(min_max_scale=False, show=False)
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    return res


def run_baseline_uncorrected(fixed: dict, force: bool = False) -> dict | None:
    """Uncorrected log-normalised PCA + UMAP + scIB on the same input as scVI."""
    figures_root = Path(fixed["figures_root"])
    out_dir = figures_root / "baseline_uncorrected"
    out_dir.mkdir(parents=True, exist_ok=True)

    scib_csv = out_dir / "scib_metrics.csv"
    umap_meta = out_dir / "umap_metadata.png"
    manifest = out_dir / "manifest.json"

    if scib_csv.exists() and umap_meta.exists() and not force:
        print("[baseline_uncorrected] already complete — skip (use --force to redo)")
        if manifest.exists():
            return json.loads(manifest.read_text())
        return None

    t0 = time.time()
    print(f"[baseline_uncorrected] start  device=cpu (PCA only)")

    adata = ad.read_h5ad(fixed["object_a_path"])
    print(f"[baseline_uncorrected] loaded {adata.shape}  X_is_raw={X_is_raw(adata)}")

    _prep_for_hvg(adata, fixed)

    hvg_batch = fixed["hvg_batch_key"] if fixed["hvg_batch_key"] in adata.obs.columns else None
    n_hvg = fixed.get("baseline_n_top_hvgs", 2000)
    sc.pp.highly_variable_genes(
        adata, flavor=fixed["hvg_flavor"], n_top_genes=n_hvg,
        layer="counts", batch_key=hvg_batch, subset=True,
    )
    print(f"[baseline_uncorrected] HVG-subset {adata.shape}")

    adata.X = adata.layers["counts"].astype(np.float32).copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.scale(adata, max_value=10, zero_center=True)
    sc.tl.pca(adata, n_comps=50, random_state=fixed["random_state"])
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=50, metric="minkowski")
    sc.tl.umap(adata, random_state=fixed["random_state"], min_dist=0.2, spread=2.0)

    # PCA scatter (metadata-coloured)
    colors = [c for c in COLOR_KEYS if c in adata.obs.columns]
    with plt.rc_context():
        sc.set_figure_params(dpi=300, figsize=(6, 6))
        sc.pl.pca(adata, color=colors, ncols=3, size=2, frameon=False, show=False)
        plt.savefig(out_dir / "pca_metadata.png", bbox_inches="tight")
        plt.close()

    _save_umap_panels(adata, out_dir, prefix="umap")

    # scIB on PCA — direct comparator for X_scvi later
    res = _scib_benchmark(
        adata, embedding_keys=["X_pca"], fixed=fixed,
        out_csv=scib_csv, plot_path=out_dir / "scib_benchmark.png",
    )
    print(f"[baseline_uncorrected] scIB:\n{res}")

    summary = {
        "name": "baseline_uncorrected",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_top_hvgs": int(n_hvg),
        "wallclock_min": float((time.time() - t0) / 60),
        "scib": res.iloc[0].to_dict() if len(res) else {},
    }
    manifest.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[baseline_uncorrected] done in {summary['wallclock_min']:.1f} min")

    del adata
    gc.collect()
    return summary


def run_one(cfg: dict, fixed: dict, force: bool = False) -> dict:
    """Train scVI for a single config; compute scIB; return summary dict."""
    name = cfg["name"]
    models_root = Path(fixed["models_root"])
    figures_root = Path(fixed["figures_root"])

    fig_dir = figures_root / name
    model_dir_root = models_root / name
    fig_dir.mkdir(parents=True, exist_ok=True)
    model_dir_root.mkdir(parents=True, exist_ok=True)

    integrated = model_dir_root / "integrated.h5ad"
    model_dir = model_dir_root / "model"
    scib_csv = fig_dir / "scib_metrics.csv"
    benchmark_plot_path = fig_dir / "scib_benchmark.png"
    manifest = fig_dir / "manifest.json"

    if integrated.exists() and scib_csv.exists() and not force:
        print(f"[{name}] already complete — skip (use --force to redo)")
        if manifest.exists():
            return json.loads(manifest.read_text())
        return {"name": name, "status": "skipped_complete"}

    # Fast resume: model + integrated exist but scIB never ran.
    if (
        not force
        and integrated.exists()
        and not scib_csv.exists()
        and model_dir.exists()
    ):
        print(f"[{name}] resume — integrated.h5ad present, re-running scIB only")
        adata = ad.read_h5ad(integrated)
        res = _scib_benchmark(adata, ["X_scvi"], fixed, scib_csv, benchmark_plot_path)
        summary = _make_manifest(name, cfg, adata, n_epochs_run=None,
                                 t0=time.time(), accel=device(), res=res)
        manifest.write_text(json.dumps(summary, indent=2, default=str))
        return summary

    t0 = time.time()
    print(f"[{name}] start  device={device()}  cfg={cfg}")

    adata = ad.read_h5ad(fixed["object_a_path"])
    print(f"[{name}] loaded {adata.shape}  X_is_raw={X_is_raw(adata)}  layers={list(adata.layers.keys())}")

    _prep_for_hvg(adata, fixed)

    batch_key = fixed["batch_key"]
    cat_keys = list(cfg.get("categorical_covariate_keys", []))
    cont_keys = list(cfg.get("continuous_covariate_keys", []))
    for col in [batch_key, fixed["label_key"], *cat_keys]:
        if col not in adata.obs.columns:
            raise RuntimeError(f"missing obs column {col!r}")

    hvg_batch = fixed["hvg_batch_key"] if fixed["hvg_batch_key"] in adata.obs.columns else None
    sc.pp.highly_variable_genes(
        adata, flavor=fixed["hvg_flavor"], n_top_genes=cfg["n_top_hvgs"],
        layer="counts", batch_key=hvg_batch, subset=True,
    )
    print(f"[{name}] HVG-subset {adata.shape} (batch_key={hvg_batch})")

    # scvi-tools setup — batch_key always Study_name; covariates per config.
    scvi.settings.seed = fixed["random_state"]
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=batch_key,
        labels_key=fixed["label_key"],
        categorical_covariate_keys=cat_keys or None,
        continuous_covariate_keys=cont_keys or None,
    )
    model = scvi.model.SCVI(
        adata,
        n_layers=fixed["n_layers"],
        n_latent=cfg["n_latent"],
        n_hidden=cfg["n_hidden"],
        gene_likelihood=cfg["gene_likelihood"],
        dropout_rate=fixed.get("dropout_rate", 0.1),
        dispersion=fixed.get("dispersion", "gene-batch"),
    )
    print(f"[{name}] model built: latent={cfg['n_latent']} hidden={cfg['n_hidden']} "
          f"layers={fixed['n_layers']} gl={cfg['gene_likelihood']} "
          f"cat={cat_keys} cont={cont_keys}")

    accel = device()
    model.train(
        max_epochs=fixed["max_epochs"],
        early_stopping=fixed.get("early_stopping", True),
        check_val_every_n_epoch=fixed.get("check_val_every_n_epoch", 1),
        batch_size=fixed.get("batch_size", 256),
        enable_progress_bar=True,
        accelerator=accel,
    )
    n_epochs_run = len(model.history["elbo_train"])
    print(f"[{name}] trained {n_epochs_run} epochs in {(time.time()-t0)/60:.1f} min")

    model.save(str(model_dir), overwrite=True)

    # Training-loss plot
    hist_train = model.history["elbo_train"].iloc[:, 0]
    hist_val = model.history["elbo_validation"].iloc[:, 0]
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(hist_train.index, hist_train.values, label="train")
    ax.plot(hist_val.index, hist_val.values, label="validation")
    ax.set_xlabel("epoch"); ax.set_ylabel("ELBO")
    ax.set_title(f"scVI — {name}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "training_loss.png")
    plt.close(fig)

    # Latent + UMAP
    adata.obsm["X_scvi"] = model.get_latent_representation(give_mean=True)
    sc.pp.neighbors(adata, use_rep="X_scvi", n_neighbors=50, metric="minkowski")
    sc.tl.umap(adata, random_state=fixed["random_state"], min_dist=0.2, spread=2.0)
    _save_umap_panels(adata, fig_dir, prefix="umap")

    adata.write_h5ad(integrated, compression="gzip")

    res = _scib_benchmark(adata, ["X_scvi"], fixed, scib_csv, benchmark_plot_path)
    print(f"[{name}] scIB:\n{res}")

    summary = _make_manifest(name, cfg, adata, n_epochs_run, t0, accel, res)
    manifest.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[{name}] done in {summary['wallclock_min']:.1f} min")

    del model, adata
    gc.collect()
    if accel == "cuda":
        torch.cuda.empty_cache()
    elif accel == "mps":
        try: torch.mps.empty_cache()
        except Exception: pass
    return summary


def _make_manifest(name, cfg, adata, n_epochs_run, t0, accel, res) -> dict:
    return {
        "name": name,
        "config": cfg,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_epochs_run": int(n_epochs_run) if n_epochs_run is not None else None,
        "wallclock_min": float((time.time() - t0) / 60),
        "device": accel,
        "scib": res.iloc[0].to_dict() if len(res) else {},
    }


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config_name", nargs="?", help="single config to run; omit if --all")
    ap.add_argument("--all", action="store_true", help="run baseline + every scVI config")
    ap.add_argument("--baseline-only", action="store_true",
                    help="only run uncorrected PCA/UMAP baseline")
    ap.add_argument("--no-baseline", action="store_true",
                    help="skip the uncorrected baseline step")
    ap.add_argument("--force", action="store_true", help="re-run even if outputs exist")
    args = ap.parse_args(argv)

    fixed, configs = load_configs()
    Path(fixed["figures_root"]).mkdir(parents=True, exist_ok=True)
    Path(fixed["models_root"]).mkdir(parents=True, exist_ok=True)

    if args.baseline_only:
        run_baseline_uncorrected(fixed, force=args.force)
        return 0

    if args.all:
        targets = configs
    elif args.config_name:
        targets = [c for c in configs if c["name"] == args.config_name]
        if not targets:
            print(f"unknown config {args.config_name}; valid: {[c['name'] for c in configs]}")
            return 2
    else:
        ap.print_usage()
        return 2

    if not args.no_baseline:
        run_baseline_uncorrected(fixed, force=args.force)

    summaries = []
    for c in targets:
        s = run_one(c, fixed, force=args.force)
        summaries.append(s)

    if args.all:
        sweep_df = pd.json_normalize(summaries, sep="__")
        out = Path(fixed["figures_root"]) / "sweep_summary.csv"
        sweep_df.to_csv(out, index=False)
        print("wrote", out)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
