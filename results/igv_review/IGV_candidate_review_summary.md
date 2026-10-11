# HCC1395 Tumor–Normal WES: IGV Candidate Review Summary

## Purpose

This report consolidates manual IGV inspection and VCF evidence for five
technically prioritized, HIGH-impact candidate variants from the HCC1395
tumor–normal whole-exome sequencing workflow.

This is an educational and benchmarking analysis, not a clinical report.
Read-level support and tumor–normal differences do not by themselves establish
pathogenicity, functional effects, or cancer-driver status.

## Candidate summary

| Gene | GRCh38 locus | Variant type / predicted consequence | Tumor AD (ref, alt) | Tumor AF | Normal AD (ref, alt) | IGV observations | Important caveats |
|---|---|---|---|---:|---|---|---|
| PLG | chr6:160711032 | 13-bp deletion; splice-acceptor / polypyrimidine-tract annotations | 5, 368 | 0.989 | 452, 0 | Multiple high-MAPQ tumor reads show a 13D CIGAR pattern; inspected normal reads align without the deletion. | Confirm transcript-specific splice effects and local alignment behavior. |
| C8orf76 | chr8:123220221 | 7-bp insertion; predicted frameshift | 0, 195 | 0.995 | 78, 0 | Insertion GGGCTAC observed in tumor reads, including an insertion CIGAR and a soft-clipped read; inspected normal reads show reference A. | STR/repeat context may complicate alignment; validate the indel and its representation. |
| TRMT6 | chr20:5943958 | C>A SNV; predicted stop-gained | 74, 128 | 0.622 | 200, 0 | Tumor reads show alternate A with high mapping and base qualities; inspected normal read shows reference C. | The manually inspected alternate reads were all F2R1, although VCF orientation counts show support in both orientation groups. Review broader read evidence. |
| MUC16 | chr19:8969814 | 17-bp deletion; predicted frameshift | 0, 96 | 0.991 | 157, 0 | Multiple high-MAPQ tumor reads show 17D CIGAR patterns without clipping; the normal track was reported to show reference sequence across the region. | MUC16 is large and repetitive; inspect local sequence complexity and confirm the indel representation. |
| TRIM36 | chr5:115130597 | 7-bp deletion; predicted frameshift | 0, 94 | 0.990 | 117, 0 | Tumor reads show 7D CIGAR patterns in F1R2 and F2R1 examples; an inspected normal read spans the region without the deletion. | Confirm transcript context and the variant with independent evidence where possible. |

**AD** = allele depths for reference and alternate alleles, respectively.
**AF** = reported tumor alternate-allele fraction. Values are reproduced from
the reviewed VCF evidence and may reflect caller-specific genotype fields.

## Read-level review highlights

### PLG
- Multiple inspected tumor reads support the same 13-base deletion.
- Inspected reads have high mapping qualities and no obvious clipping at the
  deletion boundaries.
- Inspected normal reads did not show the deletion.
- The variant is strongly tumor-enriched in the reviewed call evidence.

### C8orf76
- The insertion sequence is GGGCTAC.
- Tumor evidence includes reads with an explicit 7I CIGAR and a soft-clipped
  read at the locus.
- Inspected normal reads show the reference allele.
- Repeat context makes additional validation particularly important.

### TRMT6
- The alternate A is observed in inspected tumor reads at good mapping and
  base quality.
- The inspected normal read carries the reference C.
- The predicted stop-gained consequence is an annotation, not direct evidence
  of a protein-level effect.

### MUC16
- The corrected net deletion length is **17 bases**.
- Three inspected tumor reads show a 17D CIGAR pattern with high mapping
  quality and no clipping.
- The inspected normal reads near the locus carry reference bases; the normal
  AD field reports zero alternate reads.
- Because MUC16 is repetitive, the call should be interpreted cautiously.

### TRIM36
- Two inspected tumor reads support the same 7-base deletion and represent
  both F1R2 and F2R1 orientations.
- The inspected normal read spans the region without a deletion.
- The frameshift consequence is predicted and requires biological
  interpretation before drawing functional conclusions.

## Interpretation and limitations

1. These five variants were selected for manual review from a technically
   prioritized HIGH-impact set. The selection is not a clinical classification.
2. High mapping quality, consistent indel CIGAR patterns, good base qualities,
   and a tumor–normal difference support technical plausibility, but do not
   eliminate all mapping, repeat-context, or calling artifacts.
3. Allele fractions near 1.0 do not, by themselves, establish clonality,
   loss of heterozygosity, or a driver role. Purity, copy number, local
   mapping, and assay characteristics can affect allele fractions.
4. Predicted HIGH-impact consequences do not prove that a transcript or protein
   is functionally disrupted.
5. Further work could include independent variant validation, broader read
   review, transcript-aware annotation, population and clinical database
   checks, and relevant functional or literature evidence.

## Individual review files

- `PLG_chr6_160711032_review.txt`
- `C8orf76_chr8_123220221_review.txt`
- `TRMT6_chr20_5943958_review.txt`
- `MUC16_chr19_8969814_review.txt`
- `TRIM36_chr5_115130597_review.txt`

## Project status

Manual IGV review of these five candidates has been documented. This report
summarizes the observations available during review; it does not replace the
underlying VCF, annotation table, BAM files, or individual review notes.
