import pandas as pd
import numpy as np

input_file = "results/HCC1395.PASS.integrated_QC.tsv"
output_file = "results/HCC1395.PASS.evidence_profile.tsv"

df = pd.read_csv(input_file, sep="\t")

# --------------------------------------------------
# Convert numeric columns
# --------------------------------------------------

numeric_cols = [
    "DP",
    "AF",
    "ALT_READS",
    "NORMAL_DP",
    "NORMAL_AF",
    "NORMAL_ALT_READS",
    "TLOD",
    "NLOD",
    "NALOD",
    "GERMQ",
    "POPAF",
    "CONTQ",
    "ROQ",
    "SEQQ",
    "MMQ",
    "MPOS",
]

for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

# --------------------------------------------------
# Convert literal "." to missing values
# --------------------------------------------------

for col in ["CONTQ", "ROQ", "SEQQ", "STRANDQ", "PON", "STR"]:
    if col in df.columns:
        df[col] = df[col].replace(".", np.nan)

# --------------------------------------------------
# Basic information
# --------------------------------------------------

print("\n=== MUTECT2 EVIDENCE PROFILE ===")

print(f"Total PASS variants: {len(df)}")

# --------------------------------------------------
# Tumor evidence
# --------------------------------------------------

print("\n--- Tumor evidence ---")

print("Tumor DP:")
print(df["DP"].describe())

print("\nTumor AF:")
print(df["AF"].describe())

print("\nTumor ALT reads:")
print(df["ALT_READS"].describe())

# --------------------------------------------------
# Normal evidence
# --------------------------------------------------

print("\n--- Normal evidence ---")

print("Normal DP:")
print(df["NORMAL_DP"].describe())

print("\nNormal AF:")
print(df["NORMAL_AF"].describe())

print("\nNormal ALT reads:")
print(df["NORMAL_ALT_READS"].describe())

# --------------------------------------------------
# Mutect2 likelihood evidence
# --------------------------------------------------

print("\n--- Mutect2 likelihood evidence ---")

for col in ["TLOD", "NLOD", "NALOD"]:
    print(f"\n{col}:")
    print(df[col].describe())

# --------------------------------------------------
# Germline / population evidence
# --------------------------------------------------

print("\n--- Germline / population evidence ---")

for col in ["GERMQ", "POPAF"]:
    print(f"\n{col}:")
    print(df[col].describe())

# --------------------------------------------------
# Mapping / base / position evidence
# --------------------------------------------------

print("\n--- Read-level evidence ---")

for col in ["ROQ", "MMQ", "MPOS"]:
    print(f"\n{col}:")
    print(df[col].describe())

# --------------------------------------------------
# Missing-value assessment
# --------------------------------------------------

print("\n--- Missing values ---")

for col in [
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
    "MPOS",
    "PON",
    "STR",
]:
    if col in df.columns:
        missing = df[col].isna().sum()
        print(f"{col}: {missing} missing ({missing/len(df)*100:.2f}%)")

# --------------------------------------------------
# Existing FILTER status
# --------------------------------------------------

print("\n--- FILTER status ---")

print(df["FILTER"].value_counts(dropna=False))

# --------------------------------------------------
# Save cleaned evidence table
# --------------------------------------------------

df.to_csv(output_file, sep="\t", index=False)

print("\nSaved:")
print(output_file)
