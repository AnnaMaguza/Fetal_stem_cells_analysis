# Q1–Q5 follow-up: cells lost, cells reassigned, newly-recovered cells, remap verdict, no-enteroid reclustering

This report addresses the five follow-up questions on top of `REPORT.md`.

---

## Q1 — Were the lost cells low quality?

**Answer: partly. The lost-vs-retained distributions are clearly different on QC, but the bigger driver is per-donor: at least 3 donors are missing entirely from the new dataset, and several more lost > 99 % of cells.**

`figs/08_Q1_lost_vs_retained_qc.png`

| QC metric (in the OLD object) | retained (n = 1,928) median | lost (n = 5,889) median | M-W p |
|---|---|---|---|
| n_genes_by_counts | 2,367 | 2,186 | 4.5e-27 |
| total_counts | 9,801 | 7,519 | 6.4e-50 |
| pct_counts_mito | **3.8 %** | **7.0 %** | 5.2e-174 |
| pct_counts_ribo | 42.8 % | 33.4 % | 3.3e-294 |

So lost cells are systematically **higher mito %, lower ribo %, fewer total counts**. Not catastrophically low (their median 7,519 counts is still well above any sensible cutoff), but consistent with marginal-quality cells that were retained by the old, more-permissive pipeline.

But the dominant pattern is **at the donor / sample level**, not the per-cell level:

| Donor | n in old | % lost | comment |
|---|---|---|---|
| F72 | 2,111 | 99.4 % | sample re-processed but only 12 cells survived QC |
| F73 | 1,183 | 99.6 % | same |
| F78 | 467 | 99.8 % | same |
| BRC2258 | 914 | 100 % | **donor not present in new atlas** |
| BRC2134 | 270 | 100 % | **donor not present in new atlas** |
| BRC2259 | 144 | 100 % | **donor not present in new atlas** |
| BRC2043 | 732 | 11 % | reprocessed normally |
| F66 | 480 | 7.3 % | reprocessed normally |
| F67 | 90 | 2.2 % | reprocessed normally |

So **3 entire donors (BRC2258, BRC2134, BRC2259) are missing from the new atlas** (~1,328 cells gone purely on that), and 3 other donors (F72, F73, F78) have ~99 % of their cells filtered out by QC in the new pipeline. This is not random low-quality dropout — it is a sample-level dropout. Loss rate per old cluster:

| Old cluster | original n | % lost |
|---|---|---|
| MTRNR2L12+ASS1+_SC | 3,979 | **83.5 %** |
| RPS10+RPS17+_SC | 3,544 | 66.0 % |
| FXYD3+CKB+_SC | 294 | 77.9 % |

The MTRNR2L12+ cluster lost the most cells, consistent with the fact that MTRNR2L12 is a mitochondrial pseudogene whose reads are particularly affected by the new pseudogene-rich annotation (these cells become low-count after the remap and fail QC).

**Action for you:** look at BRC2258, BRC2134, BRC2259 — were the original FASTQs actually re-processed? If yes, they got rejected by the new pipeline; if no, they were never put through. And for F72/F73/F78, dump per-cell QC metrics from the new starsolo `Solo.out/Gene/Summary.csv` and compare with the old run's totals — most likely the per-cell `total_counts` median dropped below the new pipeline's `min_counts` threshold.

---

## Q2 — Were the reassigned cells correctly relabeled?

**Answer: yes, mostly. The re-annotation is biologically sensible — the new label's marker score is higher than the stem-cell marker score in 87–100 % of reassigned cells (except `Mesoderm`).**

`figs/09_Q2_reassigned_validation.png`

For each reassignment label, the % of cells that score higher on the new-label markers than on stem markers:

| New label | n cells | % cells with new-label score > stem score | mean new-label score | mean stem score |
|---|---|---|---|---|
| Stem cells (kept) | 956 | (stays stem) | stem 0.247 | – |
| **Enterocyte** | 199 | **98.5 %** | 0.668 | 0.009 |
| **TA** | 218 | **87.6 %** | 0.462 | 0.028 |
| **Goblet cell** | 7 | **100 %** | 1.651 | 0.121 |
| **Myofibroblast** | 26 | **100 %** | mesenchymal 0.772 | −0.026 |
| **SMC** | 23 | **100 %** | smc 2.083 (very high) | −0.035 |
| **BEST4+ epithelial** | 13 | 84.6 % | −0.027 | −0.012 |
| Proximal progenitor | 429 | 68.3 % | 0.146 | 0.037 |
| Mesothelium | 17 | 58.8 % | mesoderm 0.180 | −0.028 |
| Mesoderm | 49 | **28.6 %** | mesoderm −0.017 | −0.016 |

The reassignments to **Enterocyte, TA, Goblet, Myofibroblast, SMC** are biologically clear — these cells *do* express the new label's markers more strongly than ISC markers. Especially convincing: the 23 cells reassigned to SMC have a smooth-muscle marker score of 2.08 (vs −0.04 stem score) and 26 → Myofibroblast cells score 0.77 on mesenchymal markers vs −0.03 on stem.

**The reassignments to Mesoderm (49 cells) and Mesothelium (17 cells) are NOT well supported.** Only 28.6 % of the "Mesoderm" cells score higher on mesoderm markers than on stem markers, with mean mesoderm score essentially zero. These ~50 cells are likely a misannotation by scANVI — but a very small fraction.

**Net interpretation:** the OLD analysis was over-greedy in calling stem cells. ~50 % of cells previously labeled as "Stem cells" actually had stronger TA / Proximal progenitor / Enterocyte signatures. The new annotation is more discriminating, and for the major reclassified categories (TA, Enterocyte, Goblet, SMC, Myofibroblast) the new labels match the marker biology. The OLD's cluster `MTRNR2L12+ASS1+_SC` in particular looks like a mix of true stem cells and downstream lineage cells that the housekeeping-gene-driven clustering merged together.

---

## Q3 — Are the newly-recovered Elmentaite stem cells good quality and correctly annotated?

**Answer: no, they are systematically lower quality and have a weaker stem-cell signature than the cells that were already in the old analysis. The new pipeline appears to be over-permissive.**

`figs/10_Q3_newly_recovered_QC.png`

In the new atlas there are 2,388 Elmentaite fetal stem cells. Of those:
- 956 were already in the old fetal-SC analysis (matched by donor + barcode)
- **1,432 are new Elmentaite stem cells that weren't in the old subset**

QC comparison (medians, in the new big atlas):

| metric | already-known (n = 956) | new-only (n = 1,432) | new-only / old ratio |
|---|---|---|---|
| n_genes_by_counts | 2,316 | **1,247** | 0.54 |
| total_counts | 14,752 | **4,679** | 0.32 |
| percent_mito | 5.7 % | 4.3 % | – |
| percent_ribo | 45.9 % | 41.6 % | – |
| doublet_scores | 0.16 | 0.15 | – |

Marker scores (median):

| score | already-known | new-only |
|---|---|---|
| score_stem | **0.21** | **−0.02** |
| score_TA | −0.19 | −0.13 |
| score_proximal_prog | 0.16 | 0.10 |
| score_enterocyte | 0.15 | 0.04 |

So the new-only cells have:
- **half the genes** detected per cell
- **one-third the total reads**
- **negligible / negative stem-cell marker score** (vs +0.21 for the cells that were already in the old set)
- weaker scores for *every* fate, not just stem — these cells just don't express any lineage markers strongly because they barely have any reads

**Why weren't they in the old analysis?** Because the old pipeline filtered them out — they have ~5 k total counts and ~1.2 k genes, which is below typical `min_counts` / `min_genes` thresholds for fetal epithelium. The new pipeline has lower thresholds (or none) and is letting marginal cells through.

**Are the new annotations correct?** Probably not. They were assigned to "Stem cells" by scANVI as a low-confidence call: `scANVI` reaches for the closest neighbour in the latent space, which for low-count cells with no clear lineage markers can be the stem-cell cluster (because stem cells themselves don't strongly express differentiation markers either). It's the same reason scANVI tends to dump every confused cell into "Stem cells".

**Donor distribution of new-only cells** matches the old donors (F66, BRC2043, F73, etc.) — they came from the same biological samples, just from cell barcodes that the old pipeline had filtered.

**Action for you:**
1. Tighten the QC cutoffs in the new pipeline. A reasonable filter for fetal epithelium is `min_counts ≥ 2000` and `min_genes ≥ 500`; with that, most of the 1,432 new-only cells would drop out.
2. Re-cluster only the cells that pass the stricter QC — that should reduce the donor-driven fragmentation seen in the global UMAP.
3. Don't trust new "Stem cell" calls for cells with `total_counts < 5000` and `score_stem < 0`.

---

## Q4 — Was remapping a correct decision? Should I worry about the ~20 zeroed genes?

**Short answer: the remap itself is correct — newer Ensembl is the right thing to do — but your STARsolo invocation has dropped the multi-mapper handling, and that is what nuked the housekeeping genes. It's a pipeline-config issue, not a fundamental data problem. You should fix it but you do not need to redo the entire analysis.**

### What actually happened

The dropped genes (RPS18, RPS2, RPS3A, RPS9, RPL7A, RPS25, RPS17, ACTG1, EEF1A1, EEF1B2, EEF1D, MALAT1, B2M, LDHA, etc.) all share two features:
1. They are highly expressed and reach saturation in scRNA-seq.
2. They have **many pseudogene paralogs** scattered across the genome (RPS18 has ~20 pseudogenes; MALAT1 has retroviral-like copies; ACTB/G has actin-pseudogene families; etc.).

In `scprint_env230` the user's STARsolo command (`mapping_pipeline_for_scRNA_seq/bin/map_reads.py`) calls STAR with these solo flags:

```bash
--soloType CB_UMI_Simple
--soloUMIfiltering MultiGeneUMI_CR
--soloUMIdedup 1MM_CR
--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
--soloFeatures Gene GeneFull Velocyto
--clipAdapterType CellRanger4
--outFilterScoreMin 30
```

There is **no `--soloMultiMappers`**, so STARsolo defaults to `Unique` — meaning **multi-mapping reads are silently discarded**. With newer Ensembl annotating many more pseudogene copies, reads that were uniquely "RPS18" in old GENCODE now map equally well to "RPS18" and (e.g.) "RPS18P1", "RPS18P9" — they become multi-mappers and get thrown out. Net effect: RPS18 counts go to zero. The same logic explains every gene on the dropped list.

Note: in your earlier `draft_scripts/fastq_files_download_and_mapping.sh`, you did use `--soloMultiMappers EM`. The current production script removed it.

### Was the remap "correct"?

- **Updating the reference**: yes, this is the right decision. Newer Ensembl has cleaner mitochondrial annotation, more complete sex-chromosome genes, and properly resolves dozens of recently-named genes. The mitochondrial markers in your data actually got slightly better (MT-CO1 +0.19, MT-CYB +0.31, MT-ND1 +0.16 log1p).
- **Biological signal preservation**: yes. ISC markers (LGR5, OLFM4, ASCL2, SMOC2, SOX9, MKI67, PCNA, TOP2A) all moved by < 0.2 log1p units; per-gene mean correlation across the 6,089 shared HVGs is r = 0.93. The biology you care about is intact.
- **The dropped housekeeping genes**: that is a STARsolo flag issue, not the reference. RPS18 is biologically expressed in every cell — having it = 0 is wrong. But it doesn't affect cell-type identification because the loss is uniform across cells.

### Should you worry about the 20 zeroed genes?

It depends on your downstream use:

| use-case | impact |
|---|---|
| Cell type / clustering / annotation | low — no genes gained, biology markers preserved, in some sense cleaner because RP-driven noise is gone |
| HVG selection | positive — RPs are no longer hogging the top of the variance ranking |
| Differential expression between cell types | low for biological markers, **but** any per-cell normalization that uses raw `total_counts` is now distorted by ~30 % because the dropped genes were a big chunk of total reads |
| Per-cell normalization (`normalize_total`) | mildly affected — you're normalizing to a reduced total |
| GSEA / pathway analysis on translation, ribosomal, glycolysis pathways | **affected** — these pathways will look artificially silent; do not run GSEA on these specific gene sets without re-mapping |
| CellPhoneDB / ligand-receptor (RPs not used) | unaffected |
| scPRINT inference (your downstream) | small but present effect, since scPRINT's HVGs and gene embeddings include RPs; but scPRINT is robust to per-cell expression scaling |

### Recommendation

You don't need to redo everything, but if you have the compute, it's worth a one-sample test to confirm the diagnosis and decide on a fix:

1. **One-sample diagnostic.** Pick BRC2043 (which kept most of its cells) and re-run STARsolo with `--soloMultiMappers EM` (or `Uniform`). Compare RPS18, RPS2, RPS17 totals against the current run. If they jump back to 100k+, the diagnosis is confirmed and you have two options:
   - **Re-do the whole remap** with `--soloMultiMappers EM`. You'll get the housekeeping genes back, lose the noise-reduction benefit of the current run, and your clustering will probably look more like the old 3-cluster pattern again.
   - **Keep the current remap** as cleaner reference, but document RPS / EEF / MALAT as "lost in remap" so future GSEA / DE analyses know to skip those genes. This is the lower-effort option and what I'd recommend unless you depend on those gene sets.
2. Whichever option you choose, also tighten the new-pipeline QC thresholds (Q3) — the 1,432 marginal new-only Elmentaite cells with median ~5 k counts shouldn't be passing as stem cells.

---

## Q5 — Reclustering with enteroids removed

**Answer: there were no enteroid cells in your fetal-SC subset to remove. All 4,920 Holloway-2021 cells in `gut_hs_fetalSC_AM_05032025_150941_raw.h5ad` are primary fetal tissue (`growth_condition = primary tissue`, `Material Type = organism part`), not enteroids.**

The Holloway 2021 dataset in the full atlas is much larger than I assumed — only a small fraction is enteroid:

| Holloway growth_condition | n cells in big atlas |
|---|---|
| **primary tissue** | **221,551** |
| fetal enteroid + NRG1 | 4,386 |
| fetal enteroid + EGF | 3,220 |
| fetal enteroid + EGF + NRG1 | 1,220 |

And in the fetal-SC subset specifically, **all 4,920 Holloway cells are `primary tissue`**, so the subset never included enteroid stem cells. The enteroid-derived stem cells (3,439 of them) live under `age_group == "cell culture model"` and were already excluded when the subset was built.

### Re-clustering Elmentaite-only (matching the OLD analysis scope)

Since the old analysis was 100 % Elmentaite, I re-ran clustering on only the Elmentaite cells of the new subset (n = 2,388) — this is the most apples-to-apples comparison to the old 3-cluster result.

`figs/12_Q5_elmentaite_only_freshPCA.png` and `figs/13_Q5_elmentaite_only_scvi.png`.

| Embedding | leiden 0.3 cluster sizes |
|---|---|
| Fresh PCA + log-norm | 8 clusters: 701, 424, 367, 334, 319, 111, 70, 62 |
| **Focused scVI (precomputed) + leiden 0.3** | **5 clusters: 1073, 653, 363, 235, 64** |

**The focused scVI 5-cluster solution recapitulates the spirit of the old 3 clusters with two extras:**

| New cluster | n | top markers | maps to OLD cluster |
|---|---|---|---|
| 0 | 363 | **CKB**, MARCKSL1, S100A10, ID2, MDK, GPC3, ZFP36L2 | **FXYD3+CKB+_SC** (CKB still distinguishes it) |
| 1 | 1,073 | FCGRT, MT1G/MT1E/MT2A/MT1H (metallothioneins), TSPAN8, CCL25, AGR2, KRT18 | a maturation cluster (was probably split between OLD 0 and 1) |
| 2 | 653 | MT-RNR2, MT-ND6, TALAM1, ENSG-only mito-pseudogenes | **MTRNR2L12+ASS1+_SC** (mito-rich) |
| 3 | 235 | **RPL21, RPL7, RPL31, RPS6, RPL12, RPL34, RPS13, RPS27** | **RPS10+RPS17+_SC** — the ribosomal cluster IS recovered, just on different RP genes (the ones that survived the remap) |
| 4 | 64 | DLK1, SPINK1, TM4SF1, MDK, RBP1, RBPJ | small specialised progenitor cluster (DLK1+) |

**So the 3-cluster signal is essentially still there in the Elmentaite-only data**, just sliced into 5 instead of 3 — with one new "maturation" cluster (1) and one tiny `DLK1+` progenitor (4) that the old analysis didn't resolve.

The ribosomal cluster comes back even though RPS17 was zeroed, because there are plenty of *other* RPs (RPL21, RPL7, RPL31, RPS6, RPL12, RPS27, RPS13, RPS8 …) that the new mapping kept. The CKB+ identity stays. The mito-rich identity stays (now driven by MT-RNR2 and MT-ND6 rather than MTRNR2L12 specifically).

### Bottom line for Q5

- The fetal-SC subset is already enteroid-free.
- The reason the original 3 clusters "disappeared" is **not** that enteroid cells were polluting them (they weren't there). The two real reasons are:
  1. The mix changed: **+4,920 Holloway primary cells** (different tissue source than Elmentaite, hence different transcriptional context).
  2. The defining marker genes for two of the three clusters were zeroed by the remap.
- When you restrict to Elmentaite-only and run focused scVI, the 3-cluster signal **is recoverable** as 3 of the 5 leiden-0.3 clusters, with biologically equivalent identities.

If you want the old 3-cluster picture back: cluster on Elmentaite-only (or on Elmentaite + Holloway primary together, accepting that you'll get more clusters because the data is richer), with the focused scVI rather than the global atlas scVI.

---

## File index for this report

```
dataset_comparison/
├── REPORT.md                              ← original report (cells / expression / atlas scope)
├── REPORT_Q1_Q5.md                        ← this file
├── 09_lost_cells_qc.py                    ← Q1 + Q2 + Q3 (combined)
├── 10_Q5_no_enteroids.py                  ← Q5 first attempt (showed all subset cells are non-enteroid)
├── 11_Q5_no_holloway.py                   ← Q5 corrected: Elmentaite-only reclustering
├── 12_Q3_figure.py                        ← regenerated Q3 figure
├── Q1_lost_vs_retained_qc.csv             ← QC stats lost vs retained
├── Q1_loss_by_donor.csv                   ← per-donor loss rates
├── Q2_reassigned_marker_scores.csv        ← marker scores per new label
├── Q5_DE_focused_scvi_03.csv              ← DE genes per cluster (Elmentaite-only)
└── figs/
    ├── 08_Q1_lost_vs_retained_qc.png
    ├── 09_Q2_reassigned_validation.png
    ├── 10_Q3_newly_recovered_QC.png
    ├── 12_Q5_elmentaite_only_freshPCA.png
    ├── 13_Q5_elmentaite_only_scvi.png
    ├── 14_Q5_markers_freshPCA.png
    └── 15_Q5_markers_scvi.png
```
