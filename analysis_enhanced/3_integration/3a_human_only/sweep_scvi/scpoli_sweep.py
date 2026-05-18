"""Per-config scPoli sweep runner on Object A.

Runs in the `scarches_env` (Python 3.9) — invoke with:
    uv run --project /Users/am336941/uv_envs/scarches_env \
        python analysis_enhanced/scpoli_sweep.py <name> [--all] [--force]

scIB scoring is NOT done in this script — the integrated h5ad is written and
scIB is computed afterwards by a follow-up step in `lgr5_enhanced` (where
scib-metrics is installed). See run_scpoli_scib.py.

Per config writes:
    sweep_scpoli/<name>/integrated.h5ad
    sweep_scpoli/<name>/model/
    sweep_scpoli/<name>/training_loss.png
    sweep_scpoli/<name>/umap__study_states_age.png
    sweep_scpoli/<name>/manifest.json
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
import torch

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ROOT = Path(__file__).resolve().parent
CFG_PATH = ROOT / "scpoli_sweep_configs.json"

sc.settings.dpi = 120
sc.settings.dpi_save = 300
plt.rcParams.update({
    "savefig.bbox": "tight", "savefig.dpi": 300, "figure.dpi": 120,
    "font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
})


def device() -> str:
    if torch.cuda.is_available():
        return "cuda"
    if torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def load_configs() -> tuple[dict, list[dict]]:
    blob = json.loads(CFG_PATH.read_text())
    return blob["fixed"], blob["configs"]


def run_one(cfg: dict, fixed: dict, force: bool = False) -> dict:
    from scarches.models.scpoli import scPoli
    name = cfg["name"]
    out_root = Path(fixed["output_root"])
    out_dir = out_root / name
    out_dir.mkdir(parents=True, exist_ok=True)
    integrated = out_dir / "integrated.h5ad"
    manifest_path = out_dir / "manifest.json"
    if integrated.exists() and manifest_path.exists() and not force:
        print(f"[{name}] already complete — skip (--force to redo)")
        return json.loads(manifest_path.read_text())

    t0 = time.time()
    accel = device()
    print(f"[{name}] start  device={accel}  cfg={cfg}")

    adata = ad.read_h5ad(fixed["object_a_path"])
    print(f"[{name}] loaded {adata.shape}")

    # If the input already has the HVG subset (the canonical scVI-winner h5ad), skip
    # HVG selection — the gene space is shared with the scVI/scANVI integrations and
    # this avoids a numpy/skmisc binary incompatibility in scarches_env.
    if not fixed.get("skip_hvg"):
        if fixed.get("min_cells_per_gene"):
            n_before = adata.n_vars
            sc.pp.filter_genes(adata, min_cells=fixed["min_cells_per_gene"])
            if "counts" in adata.layers and adata.layers["counts"].shape[1] != adata.n_vars:
                adata.layers["counts"] = adata.X.copy()
            print(f"[{name}] gene filter min_cells={fixed['min_cells_per_gene']}: {n_before} -> {adata.n_vars}")
        sc.pp.highly_variable_genes(
            adata, flavor=fixed["hvg_flavor"], n_top_genes=fixed["n_top_hvgs"],
            layer="counts", batch_key=fixed["hvg_batch_key"], subset=True,
        )
        print(f"[{name}] HVG-subset {adata.shape}")
    else:
        print(f"[{name}] skip_hvg=True; using pre-HVG'd input {adata.shape}")

    # scPoli expects raw counts in .X
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()

    np.random.seed(fixed["random_state"])
    torch.manual_seed(fixed["random_state"])

    # Build scPoli
    model = scPoli(
        adata,
        condition_keys=list(fixed["condition_keys"]),
        cell_type_keys=list(fixed["cell_type_keys"]),
        unknown_ct_names=list(fixed["unknown_ct_names"]),
        latent_dim=cfg["latent_dim"],
        embedding_dims=cfg["embedding_dims"],
        hidden_layer_sizes=cfg["hidden_layer_sizes"],
        recon_loss=fixed["recon_loss"],
    )
    print(f"[{name}] model built: latent={cfg['latent_dim']} emb={cfg['embedding_dims']}")

    # Train
    model.train(
        n_epochs=fixed["n_epochs"],
        pretraining_epochs=fixed["pretraining_epochs"],
        early_stopping_kwargs=fixed["early_stopping_kwargs"],
        use_early_stopping=True,
    )
    n_epochs_run = len(getattr(model.trainer, "logs", {}).get("epoch_train_loss", [])) or fixed["n_epochs"]
    print(f"[{name}] trained {n_epochs_run} epochs in {(time.time()-t0)/60:.1f} min")

    # Save model + loss
    model_dir = out_dir / "model"
    model.save(str(model_dir), overwrite=True)
    # scPoli stores per-epoch losses in model.trainer.logs (dict[str, list[float]])
    logs = getattr(model.trainer, "logs", {})
    fig, ax = plt.subplots(figsize=(5, 3.5))
    plotted = False
    for k in ("epoch_train_loss", "epoch_val_loss", "epoch_kl_loss", "epoch_recon_loss"):
        if k in logs and len(logs[k]):
            ax.plot(logs[k], label=k.replace("epoch_", ""))
            plotted = True
    if plotted:
        ax.set_xlabel("epoch"); ax.set_ylabel("loss")
        ax.set_title(f"scPoli — {name}")
        ax.legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(out_dir / "training_loss.png")
    plt.close(fig)

    # Latent + UMAP
    adata.obsm["X_scpoli"] = model.get_latent(adata, mean=True)
    sc.pp.neighbors(adata, use_rep="X_scpoli", n_neighbors=15)
    sc.tl.umap(adata, random_state=fixed["random_state"])
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sc.pl.umap(adata, color="Study_name",  ax=axes[0], show=False, frameon=False)
    sc.pl.umap(adata, color="cell_states", ax=axes[1], show=False, frameon=False,
               legend_fontsize=6, legend_loc="right margin")
    sc.pl.umap(adata, color="age_group",   ax=axes[2], show=False, frameon=False)
    fig.suptitle(f"scPoli {name}  (latent={cfg['latent_dim']}, emb={cfg['embedding_dims']})")
    fig.tight_layout()
    fig.savefig(out_dir / "umap__study_states_age.png")
    plt.close(fig)

    adata.write_h5ad(integrated, compression="gzip")

    summary = {
        "name": name,
        "config": cfg,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_epochs_run": int(n_epochs_run),
        "wallclock_min": float((time.time() - t0) / 60),
        "device": accel,
    }
    manifest_path.write_text(json.dumps(summary, indent=2, default=str))
    print(f"[{name}] done in {summary['wallclock_min']:.1f} min")

    del model, adata
    gc.collect()
    if accel == "cuda": torch.cuda.empty_cache()
    elif accel == "mps":
        try: torch.mps.empty_cache()
        except Exception: pass
    return summary


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("config_name", nargs="?")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    fixed, configs = load_configs()
    if args.all:
        targets = configs
    elif args.config_name:
        targets = [c for c in configs if c["name"] == args.config_name]
        if not targets:
            print(f"unknown config; valid: {[c['name'] for c in configs]}")
            return 2
    else:
        ap.print_usage(); return 2
    for c in targets:
        run_one(c, fixed, force=args.force)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
