import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

PROJECT = Path.home() / "cancer_WES_project"

input_file = PROJECT / "results/HCC1395.PASS.variant_QC.tsv"
figure_dir = PROJECT / "results/figures"

figure_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(input_file, sep="\t")

# --------------------------------------------------
# Plot 1: Tumor allele fraction
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.hist(df["AF"], bins=50)

plt.xlabel("Tumor allele fraction (AF)")
plt.ylabel("Number of PASS variants")
plt.title("HCC1395 PASS Variant Tumor Allele Fraction")

plt.tight_layout()

plt.savefig(
    figure_dir / "HCC1395_PASS_AF_distribution.png",
    dpi=300
)

plt.close()

# --------------------------------------------------
# Plot 2: Tumor depth
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.hist(df["DP"], bins=50)

plt.xlabel("Tumor read depth (DP)")
plt.ylabel("Number of PASS variants")
plt.title("HCC1395 PASS Variant Tumor Depth")

plt.tight_layout()

plt.savefig(
    figure_dir / "HCC1395_PASS_DP_distribution.png",
    dpi=300
)

plt.close()

# --------------------------------------------------
# Plot 3: AF vs DP
# --------------------------------------------------

plt.figure(figsize=(8, 5))

plt.scatter(
    df["DP"],
    df["AF"],
    s=10,
    alpha=0.5
)

plt.xlabel("Tumor read depth (DP)")
plt.ylabel("Tumor allele fraction (AF)")
plt.title("HCC1395 PASS Variants: AF vs Depth")

plt.tight_layout()

plt.savefig(
    figure_dir / "HCC1395_PASS_AF_vs_DP.png",
    dpi=300
)

plt.close()

print("Plots created:")
print(figure_dir / "HCC1395_PASS_AF_distribution.png")
print(figure_dir / "HCC1395_PASS_DP_distribution.png")
print(figure_dir / "HCC1395_PASS_AF_vs_DP.png")
