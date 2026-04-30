"""Find stem-cell-like cells in the old dataset and compare against new+subset."""
import pandas as pd
from pathlib import Path

OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")

old = pd.read_parquet(OUT / "obs_old.parquet")
new = pd.read_parquet(OUT / "obs_new.parquet")
sub = pd.read_parquet(OUT / "obs_subset.parquet")

print("=== OLD: full Cell States distribution ===")
print(old["Cell States"].value_counts(dropna=False).to_string())

print("\n=== OLD: any 'stem' or 'fetal' label anywhere ===")
for col in ["Cell_Type", "Cell States", "Cell States GCA", "Cell States Kong", "Age_group", "Diagnosis"]:
    vals = old[col].astype(str)
    mask = vals.str.contains("stem", case=False, na=False) | vals.str.contains("fetal", case=False, na=False) | vals.str.contains("LGR5", na=False) | vals.str.contains("OLFM4", na=False)
    if mask.any():
        print(f"\n--- column: {col} ---")
        print(old.loc[mask, col].value_counts(dropna=False).head(20).to_string())
        print(f"total cells matched: {mask.sum()}")

print("\n=== NEW: stem-cell labelled cells ===")
for col in ["cell_states", "cellstates_scANVI", "cell_state"]:
    print(f"\n--- column: {col} ---")
    vc = new[col].astype(str).str.lower()
    mask = vc.str.contains("stem", na=False) | vc.str.contains("lgr5", na=False) | vc.str.contains("olfm4", na=False)
    print(new.loc[mask, col].value_counts(dropna=False).head(20).to_string())
    print(f"total: {mask.sum()}")

print("\n=== SUBSET: stem-cell labelled cells ===")
print(sub["cellstates_scANVI"].value_counts(dropna=False).to_string())
print("\nfetal sub by donor + study:")
print(sub.groupby(["Study_name","donor_id"]).size().head(30).to_string())

# Compare barcodes between subset and new big
print("\n=== overlap subset ↔ new big ===")
overlap = sub.index.intersection(new.index)
print(f"subset cells: {len(sub)}, new big cells: {len(new)}, overlap by index: {len(overlap)}")

# barcodes column
print("\n=== compare 'barcode' columns ===")
print("subset barcode head:", sub["barcode"].head(3).tolist())
print("new barcode head:", new["barcode"].head(3).tolist())

# Try to match by barcode + sample_id
def make_key(df):
    return df["barcode"].astype(str) + "|" + df["sample_id"].astype(str)

sub_keys = make_key(sub)
new_keys = make_key(new)
print(f"\nKey overlap (barcode|sample_id): {len(set(sub_keys) & set(new_keys))} / subset {len(sub_keys)}")

# Old dataset uses different barcodes - check
print("\n=== old vs new barcode formats ===")
print("old index head:", old.index[:3].tolist())
print("new index head:", new.index[:3].tolist())
print("subset index head:", sub.index[:3].tolist())

# Try strip prefixes/suffixes - old has 'AAAGCAATCCGTTGTC-1-Human_colon_16S8000511'
# new has 'AAAGAACCACTGAATC-1-HT071-LWRN-EGF_Epithelial'
import re
def strip_barcode(s):
    # extract 16-bp barcode
    m = re.match(r"^([ACGT]{12,20})", s)
    return m.group(1) if m else s

old_bc = old.index.to_series().map(strip_barcode)
new_bc = new.index.to_series().map(strip_barcode)
print(f"\nbarcode-only overlap (no sample disambiguation): old vs new = {len(set(old_bc) & set(new_bc))}")

# Better: use sample_id + barcode
print("\n=== overlap old ↔ new by (Sample_ID/sample_id + barcode) ===")
# old has 'Sample_ID': 'A32-SCL-0-SC-45N-1', new has 'sample_id': 'HT071-LWRN-EGF'
# Extract 16bp barcode from both indexes
old_keys = old_bc + "|" + old["Sample_ID"].astype(str)
new_keys2 = new_bc + "|" + new["sample_id"].astype(str)
print(f"overlap: {len(set(old_keys) & set(new_keys2))}")

# Save filtered tables for downstream
old_stem = old[old["Cell States Kong"].astype(str).str.contains("Stem", na=False) | old["Cell States"].astype(str).str.contains("Stem", na=False)]
print(f"\nOld stem cells: {len(old_stem)}")
print(old_stem["Cell States"].value_counts(dropna=False).head().to_string())
print(old_stem["Cell States Kong"].value_counts(dropna=False).head().to_string())
print(old_stem["Age_group"].value_counts().to_string())
print(old_stem["Study_name"].value_counts().to_string())
old_stem.to_parquet(OUT / "obs_old_stem.parquet")

new_stem = new[new["cellstates_scANVI"].astype(str) == "Stem cells"]
print(f"\nNew stem cells: {len(new_stem)}")
print(new_stem["age_group"].value_counts(dropna=False).to_string())
print(new_stem["Study_name"].value_counts(dropna=False).to_string())
new_stem.to_parquet(OUT / "obs_new_stem.parquet")
