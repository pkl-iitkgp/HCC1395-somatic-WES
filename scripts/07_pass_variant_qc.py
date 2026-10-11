import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "cancer_WES_project"

tumor_file = PROJECT / "results/HCC1395.PASS.tumor.tsv"
normal_file = PROJECT / "results/HCC1395.PASS.normal.tsv"

tumor = pd.read_csv(tumor_file, sep="\t")
normal = pd.read_csv(normal_file, sep="\t")

print("Tumor variants:", len(tumor))
print("Normal variants:", len(normal))

# --------------------------------------------------
# Convert numeric columns
# --------------------------------------------------

for df in [tumor, normal]:
    df["DP"] = pd.to_numeric(df["DP"], errors="coerce")
    df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
    df["GQ"] = pd.to_numeric(df["GQ"], errors="coerce")

# --------------------------------------------------
# Split AD into reference and alternate reads
# --------------------------------------------------

def split_ad(df):
    ad = df["AD"].astype(str).str.split(",", expand=True)

    df["REF_READS"] = pd.to_numeric(ad[0], errors="coerce")
    df["ALT_READS"] = pd.to_numeric(ad[1], errors="coerce")

    return df

tumor = split_ad(tumor)
normal = split_ad(normal)

# --------------------------------------------------
# Merge tumor and normal
# --------------------------------------------------

key = ["CHROM", "POS", "REF", "ALT"]

normal_subset = normal[
    key + ["DP", "AD", "AF", "GQ", "REF_READS", "ALT_READS"]
].copy()

normal_subset.columns = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "NORMAL_DP",
    "NORMAL_AD",
    "NORMAL_AF",
    "NORMAL_GQ",
    "NORMAL_REF_READS",
    "NORMAL_ALT_READS",
]

qc = tumor.merge(
    normal_subset,
    on=key,
    how="left"
)

# --------------------------------------------------
# AF bins
# --------------------------------------------------

qc["AF_BIN"] = pd.cut(
    qc["AF"],
    bins=[-float("inf"), 0.05, 0.10, 0.20, 0.50, 1.0],
    labels=[
        "<0.05",
        "0.05-0.10",
        "0.10-0.20",
        "0.20-0.50",
        ">0.50"
    ]
)

# --------------------------------------------------
# Basic summary
# --------------------------------------------------

print("\n=== BASIC SUMMARY ===")
print("Total PASS variants:", len(qc))
print("SNPs/MNPs/indels were calculated separately with bcftools.")

print("\n=== TUMOR DEPTH ===")
print(qc["DP"].describe())

print("\n=== TUMOR AF ===")
print(qc["AF"].describe())

print("\n=== TUMOR ALT READS ===")
print(qc["ALT_READS"].describe())

print("\n=== AF BINS ===")
print(qc["AF_BIN"].value_counts().sort_index())

print("\n=== LOW DEPTH ===")
print("DP < 10:", (qc["DP"] < 10).sum())
print("DP < 20:", (qc["DP"] < 20).sum())
print("DP < 30:", (qc["DP"] < 30).sum())

print("\n=== LOW AF ===")
print("AF < 0.05:", (qc["AF"] < 0.05).sum())
print("AF < 0.10:", (qc["AF"] < 0.10).sum())

print("\n=== NORMAL ALT SUPPORT ===")
print("Normal ALT reads > 0:",
      (qc["NORMAL_ALT_READS"] > 0).sum())

print("Normal AF >= 0.05:",
      (qc["NORMAL_AF"] >= 0.05).sum())

# --------------------------------------------------
# Save complete QC table
# --------------------------------------------------

output = PROJECT / "results/HCC1395.PASS.variant_QC.tsv"

qc.to_csv(output, sep="\t", index=False)

print("\nSaved:")
print(output)
