#!/usr/bin/env python3

import csv
import re
import pandas as pd

input_vcf = "annotation/vep/HCC1395.PASS.VEP.vcf"
output_tsv = "annotation/vep/results/HCC1395.clean_variants.tsv"

priority = [
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "splice_region_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "synonymous_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
]

def get_csq_fields(header_line):
    m = re.search(r'Format: ([^"]+)', header_line)
    if not m:
        raise RuntimeError("Could not find CSQ Format in VEP header")
    return m.group(1).split("|")

def choose_annotation(annotations):
    for consequence in priority:
        for ann in annotations:
            consequences = ann.get("Consequence", "").split("&")
            if consequence in consequences:
                return ann

    return annotations[0] if annotations else {}

rows = []

with open(input_vcf) as fh:
    csq_fields = None

    for line in fh:
        line = line.rstrip()

        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = get_csq_fields(line)
            continue

        if line.startswith("#CHROM"):
            header = line.split("\t")
            sample_names = header[9:]
            try:
                tumor_idx = sample_names.index("HCC1395_T") + 9
            except ValueError:
                raise RuntimeError("HCC1395_T sample not found")
            continue

        if line.startswith("#"):
            continue

        fields = line.split("\t")

        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        info = fields[7]
        fmt = fields[8]
        tumor = fields[tumor_idx]

        info_dict = {}
        for item in info.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_dict[k] = v

        fmt_keys = fmt.split(":")
        tumor_values = tumor.split(":")
        tumor_dict = dict(zip(fmt_keys, tumor_values))

        csq_string = info_dict.get("CSQ", "")
        annotations = []

        for ann_string in csq_string.split(","):
            vals = ann_string.split("|")
            ann = dict(zip(csq_fields, vals))
            annotations.append(ann)

        ann = choose_annotation(annotations)

        af = tumor_dict.get("AF", "")
        dp = tumor_dict.get("DP", "")

        rows.append({
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
            "FILTER": fields[6],
            "DP": dp,
            "AF": af,
            "TLOD": info_dict.get("TLOD", ""),
            "NLOD": info_dict.get("NLOD", ""),
            "SYMBOL": ann.get("SYMBOL", ""),
            "Consequence": ann.get("Consequence", ""),
            "IMPACT": ann.get("IMPACT", ""),
            "Feature": ann.get("Feature", ""),
            "BIOTYPE": ann.get("BIOTYPE", ""),
            "EXON": ann.get("EXON", ""),
            "HGVSc": ann.get("HGVSc", ""),
            "HGVSp": ann.get("HGVSp", ""),
            "Amino_acids": ann.get("Amino_acids", ""),
        })

columns = [
    "CHROM", "POS", "REF", "ALT", "FILTER", "DP", "AF",
    "TLOD", "NLOD", "SYMBOL", "Consequence", "IMPACT",
    "Feature", "BIOTYPE", "EXON", "HGVSc", "HGVSp",
    "Amino_acids"
]

df = pd.DataFrame(rows, columns=columns)

df.to_csv(output_tsv, sep="\t", index=False)

print("Saved:", output_tsv)
print("Total variants:", len(df))
