# HCC1395 Somatic Whole-Exome Sequencing (WES) Analysis

## Overview

This project documents an educational and computational benchmarking workflow
for analyzing somatic variants in the HCC1395 tumor/normal cell-line pair.

The analysis covers variant-level quality control, matched-normal evidence,
Mutect2 evidence summaries, VEP functional annotation, exploratory evidence
prioritization, and manual IGV review of selected candidates.

**This is a research and learning project, not a clinical diagnostic
pipeline.** Candidate variants and annotations must not be interpreted as
proof of pathogenicity, cancer-driver status, or clinical actionability.

## Study material

- Tumor sample: HCC1395_T
- Matched normal sample: HCC1395_N
- Tumor SRA run: SRR37849527
- Normal SRA run: SRR37849526
- Reference assembly: GRCh38
- Somatic variant caller: GATK Mutect2
- Functional annotation: Ensembl Variant Effect Predictor (VEP)

The input reads, reference files, alignment files, and annotation caches are
not included in this repository because of their size. Obtain the appropriate
source data and reference resources before attempting to reproduce the
workflow.

## Workflow

The overall workflow is:

1. Acquire tumor and matched-normal sequencing reads.
2. Perform raw-read quality control.
3. Trim reads and assess post-trimming quality.
4. Align reads to the GRCh38 reference.
5. Perform alignment processing and recalibration as configured in the
   original analysis.
6. Call somatic variants with Mutect2 and apply filtering.
7. Summarize PASS-variant read evidence in tumor and normal samples.
8. Integrate variant-level QC with Mutect2 evidence annotations.
9. Annotate the PASS VCF using VEP with the GRCh38 offline cache.
10. Generate consequence summaries and high-impact candidate tables.
11. Rank variants using an exploratory evidence index.
12. Review selected candidate loci in IGV.

### Reproducibility status

The repository contains the downstream QC, evidence-analysis, VEP-processing,
and candidate-review files listed below. The complete upstream acquisition,
trimming, alignment, recalibration, and Mutect2 execution scripts are not
currently included. Consequently, this repository does **not yet provide a
fully automated, end-to-end reproduction** of every workflow stage.

The VEP scripts also need their local paths and software environment adapted
before they can run on another machine.

## Analysis scripts

### Variant quality control and evidence analysis

- `scripts/07_pass_variant_qc.py` — summarizes tumor variant QC and merges
  available matched-normal evidence.
- `scripts/08_pass_variant_plots.py` — creates allele-fraction, depth, and
  allele-fraction-versus-depth plots.
- `scripts/09_mutect_evidence_qc.py` — summarizes Mutect2 evidence fields.
- `scripts/10_integrated_mutect_qc.py` — integrates tumor, normal, and Mutect2
  evidence tables.
- `scripts/11_mutect_evidence_profile.py` — profiles evidence distributions
  and missing values.
- `scripts/12_variant_evidence_ranking.py` — calculates exploratory
  percentile-based evidence rankings.

### VEP annotation and consequence analysis

- `scripts/annotation/02_vep_annotate_pass.sh` — annotates the PASS VCF using
  offline VEP.
- `scripts/annotation/03_vep_table_extraction.sh` — extracts selected VEP
  annotation fields into a tabular output.
- `scripts/annotation/04_make_clean_variant_table.py` — creates a clean
  variant table using a predefined consequence-priority rule.
- `scripts/annotation/05_vaf_quality_analysis.py` — performs additional
  allele-fraction and quality analysis.
- `scripts/annotation/06_vep_summary.py` — summarizes annotated consequences.
- `scripts/annotation/07_high_impact_variants.py` — extracts annotations
  marked HIGH impact.
- `scripts/annotation/08_high_impact_detailed.py` — produces a detailed
  high-impact annotation table.

The exploratory test scripts and backup copies are not part of the primary
documented workflow.

## Main output directories

### `annotation/vep/results/`

Contains selected VEP-derived tables, including:

- `HCC1395.clean_variants.tsv`
- `HCC1395.VEP.summary.txt`
- `HCC1395.high_impact.tsv`
- `HCC1395.high_impact_detailed.tsv`
- `HCC1395.high_impact_ranked.tsv`
- `HCC1395.high_impact_with_evidence_ranking.tsv`
- `HCC1395.IGV_candidates.tsv`

The VEP VCF and the extracted annotation table are also retained as selected
analysis outputs. The full VEP cache is not included.

### `results/`

Contains selected variant QC and evidence tables, including:

- `HCC1395.PASS.stats.txt`
- `HCC1395.PASS.variant_QC.tsv`
- `HCC1395.PASS.mutect_evidence_QC.tsv`
- `HCC1395.PASS.integrated_QC.tsv`
- `HCC1395.PASS.evidence_profile.tsv`
- `HCC1395.PASS.variant_evidence_ranking.tsv`

The `results/figures/` directory contains the three PASS-variant QC plots.

The `results/igv_review/` directory contains five candidate review notes and
`IGV_candidate_review_summary.md`.

### `qc/`

Raw and trimmed-read quality reports were generated during the original
analysis. Selected MultiQC reports may be added to the public repository if
they provide useful summary information without unnecessary duplication.

## Interpreting the results

### PASS status

PASS indicates that a variant passed the filtering criteria applied in this
analysis. It does not independently establish that the variant is a true
somatic mutation or biologically important.

### VEP consequences

VEP reports predicted consequences based on transcript annotations and the
selected reference assembly and cache. A HIGH-impact annotation is a
computational consequence category, not a direct demonstration of functional
effect or pathogenicity.

Some scripts select one transcript annotation using a predefined consequence
priority. This simplifies reporting but may omit relevant information from
other transcripts.

### Exploratory evidence index

The evidence-ranking script combines percentile ranks for tumor depth,
tumor allele fraction, tumor alternate-read count, TLOD, and inverse-ranked
normal allele-fraction and alternate-read evidence.

This index is an exploratory heuristic with equal-weighted components. It
is **not a validated somatic mutation score, statistical significance test,
clinical classification, or driver-gene prediction model**. Rankings depend
on the variants in the analyzed set and the evidence available for each
variant.

### IGV review

The IGV review notes record visual observations at selected candidate loci.
Manual inspection can help identify read-level support and potential
technical concerns, but it cannot by itself establish pathogenicity or
causality. Repetitive sequence, mapping ambiguity, strand/orientation
patterns, and other artifacts should be considered.

## Software and resources

The analysis uses or may depend on the following tools and libraries:

- GATK / Mutect2
- Ensembl VEP
- GRCh38 reference genome and matching VEP cache
- FastQC
- fastp
- MultiQC
- IGV
- Python
- pandas
- NumPy
- Matplotlib
- bcftools and samtools for relevant variant/alignment operations

Exact versions and complete command logs should be recorded before claiming
full reproducibility. The repository does not currently provide a complete
locked environment specification.

## Reproducing downstream analyses

The downstream scripts expect files to exist at the paths defined in each
script. Before running them on another system:

1. Obtain the original input data and the appropriate GRCh38 reference.
2. Recreate the upstream processing and variant-calling steps.
3. Install compatible versions of the required software.
4. Adapt project paths and VEP installation/cache settings.
5. Check that input table column names and VCF sample names match the scripts.
6. Run each stage in dependency order and compare outputs with the documented
   results.

Do not assume that a script is portable or that an output is reproducible
until these checks have been completed.

## Limitations

- This is a single tumor/normal cell-line analysis, not a cohort study.
- Cell-line findings may not generalize to patient tumors.
- Variant annotations and read-level evidence require independent validation
  before biological conclusions can be made.
- Exploratory rankings are not validated clinical or biological scores.
- The full upstream workflow and environment specification remain incomplete
  in this repository.
- ANNOVAR and Funcotator workflows were not completed and are excluded from
  the documented analysis.

## License and data provenance

A license has not yet been selected. Confirm the applicable terms for the
source data, software, and any redistributed resources before choosing a
repository license. Large sequencing and reference resources are intentionally
not redistributed here.
