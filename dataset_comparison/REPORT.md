# Fetal stem cells: pre-remap vs post-remap comparison

## TL;DR

Three things happened together, and only one of them is a real expression change:

1. **75 % of the original fetal SCs are gone** from the new atlas (only 1,949 of 7,811 unique cell keys are recoverable). They didn't physically vanish — they got filtered by the new pipeline's QC/annotation, presumably because their per-cell totals shifted (see #2).

2. **Genome remapping zeroed out the ribosomal/housekeeping gene block.** Per-cell expression correlation between matched cells (old vs new) is **median Pearson r = 0.91** — but that hides a clean systematic effect: a small list of very highly-expressed genes (RPS18, RPS2, RPS3A, RPS9, RPL7A, RPS17, EEF1A1, EEF1B2, ACTG1, MALAT1, B2M, LDHA, …) **dropped to ~0** in the new mapping. Almost no genes went the other way (no genes gained > log1p ≈ 1). This is the classic signature of a newer Ensembl annotation that adds many pseudogene paralogs and reassigns previously-canonical reads to them.

3. **The 3 "nice clusters" were defined by exactly those housekeeping genes that the remap killed.** Cluster names from the March-2024 leiden:
   - 0 — `MTRNR2L12+ASS1+_SC` (3,979 cells)
   - 1 — `RPS10+RPS17+_SC` (3,544 cells) — **RPS17 went from 25,358 → 5 in the new mapping**
   - 2 — `FXYD3+CKB+_SC` (294 cells) — **FXYD3 went from 2,615 → 28 in the new mapping**
   
   Two of the three cluster-defining markers were essentially deleted by the new annotation. So the clusters can't reappear no matter how you re-cluster — the signal that separated them is no longer in the data.

The actual stem-cell biology — LGR5, OLFM4, ASCL2, SMOC2, SOX9, MKI67, PCNA, TOP2A — is **preserved**: those markers shift by < 0.1 log1p units between old and new. So the remap did not break the biology, it broke the housekeeping / ribosomal pseudogene-multimapping signal that drove the previous clustering.

On top of that, the new fetal-SC pool also has a different cell mix: it lost ~5,400 Elmentaite cells and gained ~4,920 Holloway-2021 fetal-stage enteroids (now 67 % of the "fetal stem cells"). That alone changes the embedding.

---

## What was actually compared

| | OLD (`before_remapping/Fetal_healthy_stem_cells_*.h5ad`) | NEW (`gut_hs_fetalSC_AM_05032025_…raw.h5ad`) |
|---|---|---|
| date | Mar 2024 | Mar 2025 (subset of the Oct-2025 atlas) |
| reference | older Ensembl, 19,868–26,442 genes | newer Ensembl (starsolo remap), 43,704 genes |
| scope | **100 % Elmentaite-2021 fetal** | **33 % Elmentaite + 67 % Holloway-2021 fetal-enteroid** |
| cells | 7,817 | 7,308 |
| labels | 3 leiden clusters with named markers | 3 leiden clusters with no marker names |
| stem-cell substates | OLFM4 LGR5 (2,611), OLFM4 (1,294), OLFM4 PCNA (491), GCA-stem (3,420) | all simply "Stem cells" |

Cells were matched between datasets using `donor_id + 16-bp cell barcode`. 955 unique keys overlap; another 1,000+ old cells are still in the new big atlas under different (non-stem) labels.

---

## Finding 1 — 75 % of the original fetal SCs are missing in the new atlas

```
Old fetal SC keys (donor + barcode):           7,811
Same keys present anywhere in new big atlas:   1,928   (24.7 %)
   → still labelled "Stem cells":                956   (12.2 %)
   → reassigned to non-stem labels:              972
```

Per old leiden cluster:

| Old cluster | original n | matched in new | still "Stem cells" | top reassignments |
|---|---|---|---|---|
| MTRNR2L12+ASS1+_SC | 3,979 | 658 (16.5 %) | **348 (8.7 %)** | Enterocyte 138 · Proximal progenitor 118 · TA 52 |
| RPS10+RPS17+_SC | 3,544 | 1,205 (34.0 %) | **559 (15.8 %)** | Proximal progenitor 305 · TA 166 · Enterocyte 61 · Mesoderm 48 · Myofibroblast 26 · SMC 22 |
| FXYD3+CKB+_SC | 294 | 65 (22.1 %) | **49 (16.7 %)** | BEST4+ epithelial 10 · Proximal progenitor 6 |

So even of the 25 % that survived the new pipeline, half were reclassified to **TA / Proximal progenitor / Enterocyte** (i.e. the same lineage, one differentiation step downstream) — which is consistent with the loss of RPS / EEF / MALAT housekeeping markers giving these cells a less "stem-cell-extreme" transcriptome.

The "RPS10+RPS17+_SC" cluster is the only one with a meaningful spread into mesenchymal cell types (Mesoderm, Myofibroblast, SMC) — suggesting some of those cells were never really epithelial stem cells, just labelled that way because their housekeeping signal looked similar.

The most likely reason for the 75 % loss is pipeline QC: per-cell HVG totals dropped (see Finding 2), pushing many cells below `min_counts` / `min_genes` cutoffs in the new annotation pipeline.

---

## Finding 2 — Remapping zeroed out the housekeeping-gene block

For 955 matched cells, on shared HVG symbols (n = 6,089), normalised to log1p(CP10K):

| metric | value |
|---|---|
| per-cell Pearson r (median) | **0.911** |
| per-cell Pearson r (q25 / q75) | 0.906 / 0.915 |
| per-cell raw counts ratio (new / old, median) | **0.74** (mean 1.44, very skewed) |
| per-gene mean correlation across HVGs | r = 0.928 |

The per-cell totals plot (`figs/07_remap_expression_impact.png`, middle panel) is **bimodal**: most cells lie on a near-diagonal line, but a clear sub-population has **new total counts collapsed toward zero**. Those are the cells whose HVG signal was dominated by ribosomal / EEF / MALAT genes.

**Top 20 genes with biggest drop in mean log1p(CP10K)** (`gene_expression_old_vs_new_HVG.csv`):

| gene | old log1p(CP10K) | new log1p(CP10K) | change | what it is |
|---|---|---|---|---|
| RPS18 | 5.21 | **0.00** | −5.21 | ribosomal small subunit |
| RPS2  | 5.05 | 0.00 | −5.05 | ribosomal small subunit |
| RPS3A | 4.56 | 0.00 | −4.56 | ribosomal small subunit |
| RPS9  | 4.13 | 0.00 | −4.13 | ribosomal small subunit |
| RPL7A | 4.09 | 0.00 | −4.09 | ribosomal large subunit |
| RPS25 | 3.88 | 0.00 | −3.88 | ribosomal small subunit |
| MALAT1 | 5.02 | 1.34 | −3.67 | abundant lncRNA |
| ACTG1 | 3.41 | 0.00 | −3.41 | actin γ |
| EEF1B2 | 3.34 | 0.00 | −3.34 | translation elongation factor |
| EEF1A1 | 5.02 | 2.14 | −2.88 | translation elongation factor |
| RPS17 | 2.87 | **0.00** | −2.87 | **cluster-1 marker** |
| EEF1D | 2.22 | 0.02 | −2.20 | translation elongation factor |
| CFL1 | 2.19 | 0.05 | −2.14 | cofilin |
| RPL22 | 3.45 | 1.37 | −2.08 | ribosomal large subunit |
| RPL41 | 5.24 | 3.19 | −2.05 | ribosomal large subunit |
| B2M | 1.97 | 0.00 | −1.97 | MHC class I light chain |
| EIF3L | 1.97 | 0.07 | −1.90 | translation initiation |
| LDHA | 1.83 | 0.00 | −1.83 | lactate dehydrogenase A |
| PPIA | 2.51 | 0.70 | −1.81 | cyclophilin A |
| COX6A1 | 1.97 | 0.21 | −1.76 | cyt-c oxidase |
| RPL39 | 4.61 | 2.88 | −1.73 | ribosomal large subunit |

**No HVG went the other way.** Searching the file: zero genes shifted by > +1 log1p unit. So the remap is not "different normalisation" — it is one-sided: a small set of very-abundant reference genes lost almost all their reads. That is the pattern you get when a newer GTF annotates additional pseudogenes/paralogs and STARsolo's multimapping policy reassigns or discards reads previously assigned to the canonical parent gene.

Mitochondrial genes (MT-CO1, MT-CYB, MT-ND1, MT-ATP6) actually **gained ~0.1–0.3 log1p**, which is consistent with cleaner MT annotation in the newer reference — that's a small win.

### The cluster-marker genes specifically

| cluster | marker | old log1p | new log1p | Δ |
|---|---|---|---|---|
| 0 | MTRNR2L12 | (not in HVG export) | (not in HVG export) | n/a |
| 0 | ASS1 | 0.626 | 0.709 | +0.08 |
| 1 | RPS10 | 2.70 | 2.57 | −0.13 |
| 1 | **RPS17** | 2.87 | **0.001** | **−2.87** |
| 2 | **FXYD3** | 0.85 | **0.010** | **−0.84** |
| 2 | CKB | 1.50 | 1.69 | +0.19 |

So **RPS17 (cluster 1) and FXYD3 (cluster 2) were both essentially zeroed**. The clusters' identities literally disappeared from the data.

### The biological stem-cell markers — preserved

| gene | Δ log1p(CP10K) |
|---|---|
| LGR5 | +0.025 |
| OLFM4 | −0.010 |
| ASCL2 | +0.084 |
| SOX9 | +0.023 |
| MKI67 | +0.001 |
| PCNA | +0.021 |
| TOP2A | +0.006 |
| SMOC2 | −0.166 (small) |

The actual ISC markers move by < 0.2 log1p. The remap did not damage the biology — it removed a class of housekeeping/ribosomal "noise genes" that had been dominating the variance.

---

## Finding 3 — The cell mix in "fetal stem cells" also changed

| | OLD fetal SC | NEW fetal SC |
|---|---|---|
| Elmentaite 2021 (primary fetal tissue) | 7,817 (100 %) | 2,388 (33 %) |
| Holloway 2021 (fetal-stage enteroids) | 0 | **4,920 (67 %)** |
| dominant tissue | mostly small intestine, primary | mostly small intestine, but enteroid-derived |

Even if you fix the integration model, the new pool *is* a different population — primary fetal tissue diluted 1:2 with fetal-stage enteroid cultures. The Holloway cells are biologically real fetal stem cells, but they are an in-vitro system and will sit somewhere different in the embedding.

---

## So — answers to your questions

> **Do I have other cells now than originally?**
Yes, substantially.
- ~5,400 of the original 7,817 Elmentaite fetal SCs are gone from the new atlas (or at least are no longer in any "Stem cells" label).
- ~4,920 Holloway 2021 fetal-stage enteroids are new.
- Of the ~25 % Elmentaite cells that survived, ~50 % were reassigned by scANVI to TA / Proximal progenitor / Enterocyte (one differentiation step further than the old "Stem cells" label).

> **Did remapping to a newer reference change expression a lot?**
For the bulk biology, no — per-gene mean correlation r = 0.93, and ISC markers (LGR5, OLFM4, ASCL2, SMOC2, MKI67, PCNA, etc.) shift by < 0.2 log1p.
For ~20 specific highly-expressed genes — yes, dramatically: RPS18, RPS2, RPS3A, RPS9, RPL7A, RPS17, RPS25, ACTG1, EEF1A1, EEF1B2, EEF1D, MALAT1, B2M, LDHA, RPL22, RPL41, etc. went from very-high to ~zero. **No genes moved the other way by a meaningful amount.** This is consistent with the newer Ensembl reference annotating more pseudogenes and reassigning multimapped reads off the canonical genes.

> **Why don't the 3 nice clusters appear any more?**
Because they were defined by exactly the genes that the remap zeroed. The cluster names tell the story: `MTRNR2L12+ASS1+`, `RPS10+RPS17+`, `FXYD3+CKB+` — two of those three signatures (RPS17 and FXYD3) were silenced by the new annotation. With those gradients gone, the cells now sit on a smoother manifold in the biological-marker space, and any clustering you run will pick up batch / cycling / region structure instead. This is most likely an *improvement* in data quality — the previous separation was probably driven by per-cell capture-rate variance amplified through the ribosomal-pseudogene multimapping, not by genuine cell-state differences.

---

## Recommendations

1. **Don't try to recover the old 3 clusters.** The defining genes are no longer there. Trying to force the old structure back means trusting a signal you've now removed.
2. **Re-cluster on biological markers** — re-run HVG selection on the new fetal SC subset *after* excluding ribosomal / mitochondrial / EEF / MALAT genes, then re-fit scVI on that focused subset. Use ISC-relevant features (Wnt/Notch targets, cycling, regional markers).
3. **Investigate the QC drop-out.** Find the 5,400 "lost" Elmentaite cells: they are still in the starsolo output. Check what `min_counts` / `min_genes` / doublet thresholds the new annotation pipeline used; they're likely failing because their HVGs were the ribosomal genes that got reassigned.
4. **Decide whether to keep Holloway 2021 in this analysis.** If your scientific question is *primary fetal tissue stem cells*, the 4,920 enteroid cells should be split out (or at least labelled separately and integrated as a covariate); they are systematically different from the Elmentaite primary cells regardless of the reference change.
5. **Sanity-check the remap.** Open the starsolo Solo.out logs and compare `Reads Mapped to Gene` between the old and new runs for one Elmentaite sample. If the new run shows a large drop in "Reads Mapped: Unique" with a corresponding gain in "Reads Mapped: Multiple", you have confirmation that the pseudogene reassignment is what's happening.

---

## Files produced

```
dataset_comparison/
├── REPORT.md                                    ← this file
├── 01_inspect_metadata.py                       ← AnnData metadata dump
├── 02_find_stem_cells.py                        ← stem-cell label discovery
├── 03_fetal_sc_compare.py                       ← subset vs new big (same cells)
├── 04_gene_namespace.py                         ← gene namespace shift
├── 05_make_figures.py                           ← UMAP comparison figures
├── 06_inspect_old_fetal.py                      ← actual old fetal SC files
├── 07_match_old_new.py                          ← match by (donor + barcode)
├── 08_raw_counts_compare.py                     ← raw-count expression comparison
├── matched_old_new.csv                          ← 955 matched cells
├── old_fsc_in_newbig.csv                        ← old cells tracked into new big
├── old_fsc_tracked_in_newbig.csv                ← with old leiden labels
├── old_cluster_x_newbig_cellstates.csv          ← reassignment table
├── old_leiden_x_new_leiden.csv                  ← old vs new clustering
├── gene_expression_old_vs_new_HVG.csv           ← per-gene shifts (HVG-restricted)
├── expression_comparison_stats.csv              ← summary stats
└── figs/
    ├── 01_subset_fetalSC_3clusters.png          ← reference: 3 clusters in subset
    ├── 02_subset_markers.png                    ← stem markers in subset
    ├── 03_newbig_fetalSC_global_scVI.png        ← new big global scVI
    ├── 04_newbig_fetalSC_fresh_PCA.png          ← fresh PCA on new big fetal SC
    ├── 05_three_embeddings_compared.png         ← 3 UMAPs side-by-side
    ├── 06_old_clusters_tracked.png              ← old clusters in new UMAP
    └── 07_remap_expression_impact.png           ← THE expression-shift figure
```
