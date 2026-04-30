"""Compare gene namespaces between old and new datasets."""
import pandas as pd
import re
from pathlib import Path

OUT = Path("/Users/am336941/Library/CloudStorage/OneDrive-GSK/Desktop/Fetal_stem_cells_analysis/dataset_comparison")

old_v = pd.read_parquet(OUT / "var_old.parquet")
new_v = pd.read_parquet(OUT / "var_new.parquet")
sub_v = pd.read_parquet(OUT / "var_subset.parquet")

print(f"old genes: {len(old_v)}")
print(f"new genes: {len(new_v)}")
print(f"subset genes: {len(sub_v)}")
print(f"\nold var index head: {list(old_v.index[:5])}")
print(f"new var index head: {list(new_v.index[:5])}")

# Old has gene_id-query as Ensembl
print("\nold var columns:", list(old_v.columns))
print("new var columns:", list(new_v.columns))

# Old: index is gene symbol, gene_id-query is Ensembl
# New: index is mixed (symbol or Ensembl)
ensg_pattern = re.compile(r"^ENSG\d{10,}$")
new_index_str = new_v.index.astype(str)
new_is_ensembl = new_index_str.str.match(ensg_pattern.pattern)
print(f"\nNew genes with Ensembl ID as index: {new_is_ensembl.sum()} / {len(new_v)}")
print(f"New genes with symbol as index: {(~new_is_ensembl).sum()} / {len(new_v)}")

# Build mapping table for old via gene_id-query
old_ensembl = old_v["gene_id-query"].astype(str)
old_symbol = old_v.index.astype(str)

# Check overlap by gene symbol (var index)
sym_overlap = set(old_symbol) & set(new_index_str)
print(f"\nGene symbols overlapping (old index ∩ new index): {len(sym_overlap)}")

# Check Ensembl overlap
ens_overlap = set(old_ensembl) & set(new_index_str[new_is_ensembl])
print(f"Ensembl IDs overlapping (old gene_id-query ∩ new Ensembl-only entries): {len(ens_overlap)}")

# Genes in old not in new (by symbol)
only_old = set(old_symbol) - set(new_index_str)
only_new = set(new_index_str) - set(old_symbol)
print(f"\nIn old only (by symbol): {len(only_old)} (e.g. {list(only_old)[:10]})")
print(f"In new only: {len(only_new)} (e.g. {list(only_new)[:10]})")

# How many of the "only new" are Ensembl IDs (presumably no symbol mapped)?
new_only_arr = pd.Series(list(only_new))
new_only_is_ens = new_only_arr.str.match(ensg_pattern.pattern)
print(f"  ↳ of which are Ensembl-IDs (no symbol): {new_only_is_ens.sum()}")
print(f"  ↳ of which are gene symbols not in old: {(~new_only_is_ens).sum()}")
print(f"  Sample new symbols not in old: {new_only_arr[~new_only_is_ens].sample(min(20, (~new_only_is_ens).sum()), random_state=0).tolist()}")

# Check key intestinal stem cell markers
markers = ["LGR5","OLFM4","ASCL2","SMOC2","PROM1","BMI1","HOPX","SOX9","CDX2","VIL1","MKI67","PCNA","TOP2A",
           "MUC2","TFF3","CHGA","DEFA5","LYZ","ALPI","KRT20","SI","SLC15A1","SLC5A1",
           "NOTCH1","NOTCH2","DLL1","DLL4","HES1","WNT3","RNF43","ZNRF3","ID3","TLR2",
           "ANXA1","ANXA2","S100A6","S100A11"]
print("\n=== Intestinal stem cell marker presence ===")
for m in markers:
    in_old = m in set(old_symbol)
    in_new = m in set(new_index_str)
    print(f"  {m:10s}  old: {in_old}  new: {in_new}")
