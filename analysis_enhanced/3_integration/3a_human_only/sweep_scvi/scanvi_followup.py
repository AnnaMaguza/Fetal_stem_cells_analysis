"""Train scANVI from the winning scVI model on Object A.

Reads:
    - DATA_OUT/sweep_scvi_winner.txt          (winner name written by 3a1b)
    - DATA_OUT/sweep_scvi/<winner>/model/     (scVI model artefact)
    - DATA_OUT/sweep_scvi/<winner>/integrated.h5ad   (HVG-subset AnnData used by scVI)

Writes:
    - DATA_OUT/integrated_obj_a_scanvi.h5ad   (latent in X_scanvi + UMAP)
    - DATA_OUT/model_scanvi_obj_a_scvi/       (canonical scANVI model)
    - DATA_OUT/figures/scvi_sweep/scanvi_training_loss.png
    - DATA_OUT/figures/scvi_sweep/scanvi_umap__study_states_age.png
    - DATA_OUT/scib_metrics_obj_a_scanvi.csv

Usage:
    uv run --project /Users/am336941/uv_envs/lgr5_enhanced \
        python analysis_enhanced/scanvi_followup.py [--max-epochs 20] [--force]
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

DATA_OUT = Path("/Users/am336941/PhD/data/Fetal_stem_cells_analysis_enhanced")
REPO = Path(__file__).resolve().parent.parent
FIG_DIR = REPO / "analysis_enhanced" / "figures" / "scvi_sweep"
FIG_DIR.mkdir(parents=True, exist_ok=True)
sc.settings.dpi = 120; sc.settings.dpi_save = 300
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


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-epochs", type=int, default=20)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)

    winner_tag = DATA_OUT / "sweep_scvi_winner.txt"
    if not winner_tag.exists():
        sys.exit("sweep_scvi_winner.txt not found — run 3a1b sweep summary first")
    winner = winner_tag.read_text().strip()
    win_dir = DATA_OUT / "sweep_scvi" / winner
    scvi_model_dir = win_dir / "model"
    integrated_scvi = win_dir / "integrated.h5ad"
    print(f"winner config: {winner}")
    print(f"loading scVI model from {scvi_model_dir}")

    out_scanvi = DATA_OUT / "integrated_obj_a_scanvi.h5ad"
    if out_scanvi.exists() and not args.force:
        print(f"{out_scanvi} already exists — pass --force to redo")
        return 0

    t0 = time.time()
    adata = ad.read_h5ad(integrated_scvi)
    print("loaded", adata.shape, "layers:", list(adata.layers))

    model_dir = DATA_OUT / "model_scanvi_obj_a_scvi"
    accel = device()
    scvi.settings.seed = 0

    if model_dir.exists() and (model_dir / "model.pt").exists() and not args.force:
        # Fast resume: pre-trained scANVI is on disk — just re-attach to the same adata
        print(f"resume — loading existing scANVI model from {model_dir}")
        scvi_model = scvi.model.SCVI.load(str(scvi_model_dir), adata=adata)
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model, labels_key="cell_states", unlabeled_category="unknown",
        )
        # Replace the random-init weights with the saved scANVI checkpoint
        scanvi_model = scvi.model.SCANVI.load(str(model_dir), adata=adata)
        n_epochs_run = None
    else:
        scvi_model = scvi.model.SCVI.load(str(scvi_model_dir), adata=adata)
        print("loaded scVI model; chaining scANVI…")
        scanvi_model = scvi.model.SCANVI.from_scvi_model(
            scvi_model, labels_key="cell_states", unlabeled_category="unknown",
        )
        scanvi_model.train(max_epochs=args.max_epochs, accelerator=accel, n_samples_per_label=100)
        n_epochs_run = len(scanvi_model.history["elbo_train"])
        print(f"scANVI trained {n_epochs_run} epochs in {(time.time()-t0)/60:.1f} min")
        scanvi_model.save(str(model_dir), overwrite=True)

    # Training-loss plot (tolerant of missing validation track)
    hist = scanvi_model.history if hasattr(scanvi_model, "history") else {}
    fig, ax = plt.subplots(figsize=(5, 3.5))
    plotted = False
    for key, lbl in [("elbo_train", "train"), ("elbo_validation", "validation"),
                     ("reconstruction_loss_train", "recon train"),
                     ("classification_loss_train", "classif train")]:
        if key in hist and len(hist[key]):
            s = hist[key].iloc[:, 0]
            ax.plot(s.index, s.values, label=lbl); plotted = True
    if plotted:
        ax.set_xlabel("epoch"); ax.set_ylabel("loss")
        ax.set_title(f"scANVI on top of scVI:{winner}")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG_DIR / "scanvi_training_loss.png")
    plt.close(fig)

    # Latent + UMAP
    adata.obsm["X_scanvi"] = scanvi_model.get_latent_representation(give_mean=True)
    sc.pp.neighbors(adata, use_rep="X_scanvi", n_neighbors=15)
    sc.tl.umap(adata, random_state=0)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sc.pl.umap(adata, color="Study_name",  ax=axes[0], show=False, frameon=False)
    sc.pl.umap(adata, color="cell_states", ax=axes[1], show=False, frameon=False,
               legend_fontsize=6, legend_loc="right margin")
    sc.pl.umap(adata, color="age_group",   ax=axes[2], show=False, frameon=False)
    fig.suptitle(f"scANVI (on scVI:{winner})")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "scanvi_umap__study_states_age.png"); plt.close(fig)

    adata.write_h5ad(out_scanvi, compression="gzip")
    print("saved", out_scanvi)

    # scIB benchmark for scANVI alongside scVI
    print("running scIB on scANVI…")
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
    bio = BioConservation(
        isolated_labels=True, nmi_ari_cluster_labels_kmeans=True,
        silhouette_label=True, clisi_knn=True,
    )
    batch = BatchCorrection(
        graph_connectivity=True, ilisi_knn=True, kbet_per_label=False,
        pcr_comparison=True, bras=True,
    )
    a_lab = adata[adata.obs["cell_states"].astype(str) != "unknown"].copy()
    bm = Benchmarker(
        a_lab, batch_key="Study_name", label_key="cell_states",
        embedding_obsm_keys=["X_scvi", "X_scanvi"] if "X_scvi" in a_lab.obsm else ["X_scanvi"],
        bio_conservation_metrics=bio, batch_correction_metrics=batch, n_jobs=2,
    )
    bm.benchmark()
    res = bm.get_results(min_max_scale=False, clean_names=True)
    res.to_csv(DATA_OUT / "scib_metrics_obj_a_scanvi.csv")
    print(res)
    print(f"done in {(time.time()-t0)/60:.1f} min")

    del adata, a_lab, bm, scanvi_model, scvi_model
    gc.collect()
    if accel == "cuda": torch.cuda.empty_cache()
    elif accel == "mps":
        try: torch.mps.empty_cache()
        except Exception: pass
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
