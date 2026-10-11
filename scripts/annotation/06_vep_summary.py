#!/usr/bin/env python3

import pandas as pd

f = "annotation/vep/results/HCC1395.clean_variants.tsv"
out = "annotation/vep/results/HCC1395.VEP.summary.txt"

df = pd.read_csv(f, sep="\t")

lines = []

lines.append("HCC1395 VEP SUMMARY")
lines.append("===================")
lines.append("")
lines.append(f"Total PASS variants: {len(df)}")
lines.append(
    f"Variants with gene symbols: "
    f"{((df['SYMBOL'].notna()) & (df['SYMBOL'] != '-')).sum()}"
)
lines.append(
    f"Missense variants: "
    f"{df['Consequence'].fillna('').str.contains('missense_variant').sum()}"
)
lines.append(
    f"Synonymous variants: "
    f"{df['Consequence'].fillna('').str.contains('synonymous_variant').sum()}"
)
lines.append(
    f"Splice-related variants: "
    f"{df['Consequence'].fillna('').str.contains('splice_').sum()}"
)

lines.append("")
lines.append("IMPACT DISTRIBUTION")
lines.append("-------------------")
lines.append(df["IMPACT"].value_counts().to_string())

lines.append("")
lines.append("TOP CONSEQUENCES")
lines.append("-----------------")
lines.append(df["Consequence"].value_counts().head(20).to_string())

with open(out, "w") as fh:
    fh.write("\n".join(lines) + "\n")

print("\n".join(lines))
print("\nSaved:", out)
