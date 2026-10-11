import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "cancer_WES_project"

input_file = PROJECT / "results/HCC1395.PASS.mutect_evidence.tsv"
output_file = PROJECT / "results/HCC1395.PASS.mutect_evidence_QC.tsv"

df = pd.read_csv(
    input_file,
    sep="\t",
    header=None,
    names=[
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "FILTER",
        "TLOD",
        "NLOD",
        "NALOD",
        "GERMQ",
        "POPAF",
        "CONTQ",
        "ROQ",
        "SEQQ",
        "STRANDQ",
        "MMQ",
        "MBQ",
        "MPOS",
        "PON",
        "STR"
    ]
)

print("\n=== Mutect2 Evidence QC ===")
print(f"Total PASS variants: {len(df)}")

print("\n--- Missing values ---")
missing = df.isna().sum()
print(missing)

print("\n--- TLOD ---")
print(df["TLOD"].describe())

print("\n--- NLOD ---")
print(df["NLOD"].describe())

print("\n--- NALOD ---")
print(df["NALOD"].describe())

print("\n--- GERMQ ---")
print(df["GERMQ"].describe())

print("\n--- POPAF ---")
print(df["POPAF"].describe())

print("\n--- ROQ ---")
print(df["ROQ"].describe())

print("\n--- MPOS ---")
print(df["MPOS"].describe())

print("\n--- PON flag ---")
print(df["PON"].value_counts(dropna=False))

print("\n--- STR flag ---")
print(df["STR"].value_counts(dropna=False))

# Simple descriptive categories
print("\n=== TLOD descriptive groups ===")
print(f"TLOD < 10:  {(df['TLOD'] < 10).sum()}")
print(f"TLOD 10-20: {((df['TLOD'] >= 10) & (df['TLOD'] < 20)).sum()}")
print(f"TLOD 20-50: {((df['TLOD'] >= 20) & (df['TLOD'] < 50)).sum()}")
print(f"TLOD >=50:  {(df['TLOD'] >= 50).sum()}")

print("\n=== NLOD descriptive groups ===")
print(f"NLOD < 2:   {(df['NLOD'] < 2).sum()}")
print(f"NLOD 2-5:   {((df['NLOD'] >= 2) & (df['NLOD'] < 5)).sum()}")
print(f"NLOD >= 5:  {(df['NLOD'] >= 5).sum()}")

print("\n=== GERMQ descriptive groups ===")
print(f"GERMQ < 20: {(df['GERMQ'] < 20).sum()}")
print(f"GERMQ 20-50: {((df['GERMQ'] >= 20) & (df['GERMQ'] < 50)).sum()}")
print(f"GERMQ >=50: {(df['GERMQ'] >= 50).sum()}")

# Save cleaned evidence table
df.to_csv(output_file, sep="\t", index=False)

print("\nSaved:")
print(output_file)
