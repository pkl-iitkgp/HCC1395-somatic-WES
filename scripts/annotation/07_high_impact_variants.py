#!/usr/bin/env python3

import pandas as pd

f = "annotation/vep/results/HCC1395.clean_variants.tsv"

df = pd.read_csv(f, sep="\t")

for col in ["AF", "DP", "TLOD", "NLOD"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

high = df[df["IMPACT"] == "HIGH"].copy()

print("Total HIGH-impact variants:", len(high))

print("\nHIGH-impact consequence distribution:")
print(high["Consequence"].value_counts().to_string())

print("\nHIGH-impact variants with gene symbols:")

cols = [
    "CHROM", "POS", "REF", "ALT", "SYMBOL",
    "Consequence", "IMPACT", "AF", "DP",
    "TLOD", "NLOD", "HGVSc", "HGVSp",
    "Amino_acids"
]

print(
    high[
        high["SYMBOL"].notna() & (high["SYMBOL"] != "-")
    ][cols].to_string(index=False)
)

out = "annotation/vep/results/HCC1395.high_impact.tsv"
high.to_csv(out, sep="\t", index=False)

print("\nSaved:", out)
