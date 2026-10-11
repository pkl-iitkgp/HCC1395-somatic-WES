import pandas as pd
import numpy as np

input_file = "results/HCC1395.PASS.evidence_profile.tsv"
output_file = "results/HCC1395.PASS.variant_evidence_ranking.tsv"

df = pd.read_csv(input_file, sep="\t")

# --------------------------------------------------
# Numeric conversion
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
    "ROQ",
    "MPOS",
]

for col in numeric_cols:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")


# --------------------------------------------------
# Helper: percentile rank
# --------------------------------------------------

def percentile_rank(series):
    return series.rank(pct=True)


# --------------------------------------------------
# Evidence percentiles
#
# These are descriptive rankings only.
# They are NOT clinical thresholds.
# --------------------------------------------------

df["Tumor_DP_percentile"] = percentile_rank(df["DP"])
df["Tumor_AF_percentile"] = percentile_rank(df["AF"])
df["Tumor_ALT_percentile"] = percentile_rank(df["ALT_READS"])
df["TLOD_percentile"] = percentile_rank(df["TLOD"])

# For normal evidence, lower is more consistent
# with absence of the alternate allele.
df["Normal_AF_low_percentile"] = 1 - percentile_rank(
    df["NORMAL_AF"]
)

df["Normal_ALT_low_percentile"] = 1 - percentile_rank(
    df["NORMAL_ALT_READS"]
)

# --------------------------------------------------
# Exploratory evidence index
#
# Equal weighting of descriptive components.
# This is NOT a validated somatic mutation score.
# --------------------------------------------------

score_columns = [
    "Tumor_DP_percentile",
    "Tumor_AF_percentile",
    "Tumor_ALT_percentile",
    "TLOD_percentile",
    "Normal_AF_low_percentile",
    "Normal_ALT_low_percentile",
]

df["Exploratory_evidence_index"] = (
    df[score_columns].mean(axis=1)
)


# --------------------------------------------------
# Evidence descriptors
# --------------------------------------------------

df["Normal_evidence"] = np.select(
    [
        (df["NORMAL_ALT_READS"] == 0) &
        (df["NORMAL_AF"] < 0.05),

        (df["NORMAL_ALT_READS"] > 0) |
        (df["NORMAL_AF"] >= 0.05)
    ],
    [
        "Low_normal_ALT_support",
        "Normal_ALT_or_AF_detected"
    ],
    default="Other"
)


# --------------------------------------------------
# Tumor evidence descriptors
#
# Use quartiles rather than arbitrary biological
# thresholds.
# --------------------------------------------------

dp_q25 = df["DP"].quantile(0.25)
dp_q75 = df["DP"].quantile(0.75)

af_q25 = df["AF"].quantile(0.25)
af_q75 = df["AF"].quantile(0.75)

alt_q25 = df["ALT_READS"].quantile(0.25)
alt_q75 = df["ALT_READS"].quantile(0.75)

df["Tumor_depth_group"] = np.select(
    [
        df["DP"] <= dp_q25,
        df["DP"] >= dp_q75
    ],
    [
        "Lower_quartile",
        "Upper_quartile"
    ],
    default="Middle_50_percent"
)

df["Tumor_AF_group"] = np.select(
    [
        df["AF"] <= af_q25,
        df["AF"] >= af_q75
    ],
    [
        "Lower_quartile",
        "Upper_quartile"
    ],
    default="Middle_50_percent"
)

df["Tumor_ALT_reads_group"] = np.select(
    [
        df["ALT_READS"] <= alt_q25,
        df["ALT_READS"] >= alt_q75
    ],
    [
        "Lower_quartile",
        "Upper_quartile"
    ],
    default="Middle_50_percent"
)


# --------------------------------------------------
# Normal evidence flags
# --------------------------------------------------

df["Normal_ALT_flag"] = np.where(
    df["NORMAL_ALT_READS"] > 0,
    "ALT_detected",
    "No_ALT_detected"
)

df["Normal_AF_flag"] = np.where(
    df["NORMAL_AF"] >= 0.05,
    "AF_ge_0.05",
    "AF_lt_0.05"
)


# --------------------------------------------------
# Sort by exploratory evidence index
# --------------------------------------------------

df = df.sort_values(
    "Exploratory_evidence_index",
    ascending=False
)


# --------------------------------------------------
# Add rank
# --------------------------------------------------

df.insert(
    0,
    "Evidence_rank",
    range(1, len(df) + 1)
)


# --------------------------------------------------
# Select useful columns first
# --------------------------------------------------

priority_columns = [
    "Evidence_rank",
    "CHROM",
    "POS",
    "REF",
    "ALT",
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
    "ROQ",
    "MPOS",
    "AF_difference",
    "Normal_evidence",
    "Tumor_depth_group",
    "Tumor_AF_group",
    "Tumor_ALT_reads_group",
    "Exploratory_evidence_index",
]

remaining_columns = [
    c for c in df.columns
    if c not in priority_columns
]

df = df[priority_columns + remaining_columns]


# --------------------------------------------------
# Save
# --------------------------------------------------

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("=== VARIANT EVIDENCE RANKING ===")
print(f"Total variants: {len(df)}")

print("\nTop 20 variants by exploratory evidence index:")

print(
    df[
        [
            "Evidence_rank",
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "DP",
            "AF",
            "ALT_READS",
            "NORMAL_AF",
            "NORMAL_ALT_READS",
            "TLOD",
            "GERMQ",
            "Exploratory_evidence_index"
        ]
    ].head(20).to_string(index=False)
)

print("\n--- Normal evidence groups ---")
print(df["Normal_evidence"].value_counts())

print("\n--- Tumor depth groups ---")
print(df["Tumor_depth_group"].value_counts())

print("\n--- Tumor AF groups ---")
print(df["Tumor_AF_group"].value_counts())

print("\n--- Tumor ALT-read groups ---")
print(df["Tumor_ALT_reads_group"].value_counts())

print("\nSaved:")
print(output_file)
