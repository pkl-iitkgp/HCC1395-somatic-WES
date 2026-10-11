#!/usr/bin/env python3

import re
import pandas as pd

vcf = "annotation/vep/HCC1395.PASS.VEP.vcf"
out = "annotation/vep/results/HCC1395.high_impact_detailed.tsv"

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
]

def parse_csq_header(line):
    m = re.search(r'Format: ([^"]+)', line)
    if not m:
        raise RuntimeError("CSQ format not found")
    return m.group(1).split("|")

def best_annotation(annotations):
    # Prefer protein-coding annotations with a gene symbol
    candidates = [
        a for a in annotations
        if a.get("SYMBOL") not in ("", "-")
        and a.get("BIOTYPE") == "protein_coding"
    ]

    if not candidates:
        candidates = [
            a for a in annotations
            if a.get("SYMBOL") not in ("", "-")
        ]

    if not candidates:
        candidates = annotations

    # Prioritize severe consequences
    for consequence in priority:
        for a in candidates:
            cons = a.get("Consequence", "").split("&")
            if consequence in cons:
                return a

    return candidates[0] if candidates else {}

rows = []
csq_fields = None
tumor_idx = None

with open(vcf) as fh:

    for line in fh:
        line = line.rstrip()

        if line.startswith("##INFO=<ID=CSQ"):
            csq_fields = parse_csq_header(line)
            continue

        if line.startswith("#CHROM"):
            header = line.split("\t")
            samples = header[9:]

            if "HCC1395_T" not in samples:
                raise RuntimeError("HCC1395_T not found")

            tumor_idx = 9 + samples.index("HCC1395_T")
            continue

        if line.startswith("#"):
            continue

        f = line.split("\t")

        chrom = f[0]
        pos = f[1]
        ref = f[3]
        alt = f[4]
        filt = f[6]
        info = f[7]
        fmt = f[8]
        tumor = f[tumor_idx]

        info_dict = {}

        for item in info.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                info_dict[k] = v

        fmt_keys = fmt.split(":")
        tumor_values = tumor.split(":")
        tumor_dict = dict(zip(fmt_keys, tumor_values))

        csq = info_dict.get("CSQ", "")

        annotations = []

        for record in csq.split(","):
            values = record.split("|")
            annotations.append(dict(zip(csq_fields, values)))

        high_annotations = [
            a for a in annotations
            if a.get("IMPACT") == "HIGH"
        ]

        if not high_annotations:
            continue

        ann = best_annotation(high_annotations)

        rows.append({
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt,
            "FILTER": filt,
            "SYMBOL": ann.get("SYMBOL", ""),
            "Gene": ann.get("Gene", ""),
            "Feature": ann.get("Feature", ""),
            "BIOTYPE": ann.get("BIOTYPE", ""),
            "Consequence": ann.get("Consequence", ""),
            "IMPACT": ann.get("IMPACT", ""),
            "EXON": ann.get("EXON", ""),
            "HGVSc": ann.get("HGVSc", ""),
            "HGVSp": ann.get("HGVSp", ""),
            "Amino_acids": ann.get("Amino_acids", ""),
            "AF": tumor_dict.get("AF", ""),
            "DP": tumor_dict.get("DP", ""),
            "TLOD": info_dict.get("TLOD", ""),
            "NLOD": info_dict.get("NLOD", ""),
        })

df = pd.DataFrame(rows)

df.to_csv(out, sep="\t", index=False)

print("Saved:", out)
print("HIGH-impact variants:", len(df))
print()
print("Variants with HGVSc:", (df["HGVSc"] != "").sum())
print("Variants with HGVSp:", (df["HGVSp"] != "").sum())
print("Variants with gene:", (df["SYMBOL"] != "").sum())
print()
print(df[
    ["CHROM","POS","REF","ALT","SYMBOL",
     "Consequence","AF","DP","HGVSc","HGVSp"]
].to_string(index=False))
