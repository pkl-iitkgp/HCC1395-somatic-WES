#!/usr/bin/env python3

import pandas as pd

f = "annotation/vep/results/HCC1395.clean_variants.tsv"

df = pd.read_csv(f, sep="\t")

df["AF"] = pd.to_numeric(df["AF"], errors="coerce")

print("Total variants:", len(df))
print("Variants with AF:", df["AF"].notna().sum())

print("\nAF summary:")
print(df["AF"].describe())

print("\nAF ranges:")
print("AF < 0.05:", (df["AF"] < 0.05).sum())
print(
    "0.05 <= AF < 0.10:",
    ((df["AF"] >= 0.05) & (df["AF"] < 0.10)).sum()
)
print(
    "0.10 <= AF < 0.25:",
    ((df["AF"] >= 0.10) & (df["AF"] < 0.25)).sum()
)
print(
    "0.25 <= AF < 0.50:",
    ((df["AF"] >= 0.25) & (df["AF"] < 0.50)).sum()
)
print("AF >= 0.50:", (df["AF"] >= 0.50).sum())
