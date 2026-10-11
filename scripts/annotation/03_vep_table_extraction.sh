#!/usr/bin/env bash

cd ~/cancer_WES_project

mkdir -p annotation/vep/results

cd ~/vep/ensembl-vep

PERL5LIB=~/vep/ensembl-vep/modules:$PERL5LIB \
perl vep \
  -i ~/cancer_WES_project/annotation/vep/HCC1395.PASS.VEP.vcf \
  -o ~/cancer_WES_project/annotation/vep/results/HCC1395.VEP.table.txt \
  --tab \
  --cache \
  --offline \
  --assembly GRCh38 \
  --dir_cache ~/cancer_WES_project/annotation/vep/cache \
  --force_overwrite \
  --fields "Uploaded_variation,Location,Allele,Gene,SYMBOL,Feature,Consequence,IMPACT,BIOTYPE,EXON,HGVSc,HGVSp"
