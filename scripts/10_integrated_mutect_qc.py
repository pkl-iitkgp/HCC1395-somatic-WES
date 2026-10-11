import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "cancer_WES_project"

tumor_file = PROJECT / "results/HCC1395.PASS.variant_QC.tsv"
normal_file = PROJECT / "results/HCC1395.PASS.normal_evidence.tsv"
mutect_file = PROJECT / "results/HCC1395.PASS.mutect_evidence_QC.tsv"

output_file = PROJECT / "results/HCC1395.PASS.integrated_QC.tsv"


# -----------------------------
# Read tumor QC
# -----------------------------
tumor = pd.read_csv(tumor_file, sep="\t")

print(f"Tumor QC variants: {len(tumor)}")


# -----------------------------
# Read normal evidence
# -----------------------------
normal = pd.read_csv(
    normal_file,
    sep="\t",
    header=None,
    names=[
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "Normal_DP",
        "Normal_AD",
        "Normal_AF",
        "Normal_GQ"
    ]
)

print(f"Normal evidence variants: {len(normal)}")


# -----------------------------
# Split normal AD
# -----------------------------
normal[["Normal_REF_reads", "Normal_ALT_reads"]] = (
    normal["Normal_AD"]
    .astype(str)
    .str.split(",", expand=True)
)

normal["Normal_REF_reads"] = pd.to_numeric(
    normal["Normal_REF_reads"], errors="coerce"
)

normal["Normal_ALT_reads"] = pd.to_numeric(
    normal["Normal_ALT_reads"], errors="coerce"
)


# -----------------------------
# Read Mutect evidence
# -----------------------------
mutect = pd.read_csv(mutect_file, sep="\t")

print(f"Mutect evidence variants: {len(mutect)}")


# -----------------------------
# Merge normal evidence
# -----------------------------
keys = ["CHROM", "POS", "REF", "ALT"]

df = tumor.merge(
    normal,
    on=keys,
    how="left"
)

# -----------------------------
# Merge Mutect2 evidence
# -----------------------------
df = df.merge(
    mutect,
    on=keys,
    how="left",
    suffixes=("", "_mutect")
)


# -----------------------------
# Calculate descriptive fields
# -----------------------------

df["AF_difference"] = df["AF"] - df["Normal_AF"]

df["Normal_has_ALT"] = df["Normal_ALT_reads"] > 0

df["Normal_AF_ge_0.05"] = df["Normal_AF"] >= 0.05

df["Tumor_ALT_reads"] = df["ALT_READS"]


# -----------------------------
# Summary
# -----------------------------

print("\n=== Integrated Mutect2 QC ===")

print(f"Total integrated variants: {len(df)}")

print("\n--- Normal ALT support ---")
print(
    f"Normal ALT reads > 0: "
    f"{df['Normal_has_ALT'].sum()}"
)

print(
    f"Normal ALT reads = 0: "
    f"{(~df['Normal_has_ALT']).sum()}"
)

print("\n--- Normal AF ---")
print(
    f"Normal AF >= 0.05: "
    f"{df['Normal_AF_ge_0.05'].sum()}"
)

print(
    f"Normal AF < 0.05: "
    f"{(~df['Normal_AF_ge_0.05']).sum()}"
)

print("\n--- Tumor vs normal AF difference ---")
print(df["AF_difference"].describe())

print("\n--- Tumor ALT reads ---")
print(df["Tumor_ALT_reads"].describe())

print("\n--- Normal ALT reads ---")
print(df["Normal_ALT_reads"].describe())

print("\n--- Tumor AF ---")
print(df["AF"].describe())

print("\n--- Normal AF ---")
print(df["Normal_AF"].describe())


# -----------------------------
# Save integrated table
# -----------------------------

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("\nSaved:")
print(output_file)
