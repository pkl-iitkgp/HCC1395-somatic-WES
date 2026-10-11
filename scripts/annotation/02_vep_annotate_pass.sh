#!/usr/bin/env bash
set -euo pipefail

VEP_DIR="$HOME/vep/ensembl-vep"
VEP_PERL="$HOME/miniconda3/envs/te_project/bin/perl"
PROJECT="$HOME/cancer_WES_project"

export PERL5LIB="$VEP_DIR/modules:$HOME/miniconda3/envs/te_project/lib/perl5/site_perl/5.22.2:$HOME/miniconda3/envs/te_project/lib/perl5/site_perl/5.22.2/x86_64-linux-thread-multi${PERL5LIB:+:$PERL5LIB}"

cd "$VEP_DIR"

"$VEP_PERL" vep \
  -i "$PROJECT/variants/HCC1395.PASS.vcf.gz" \
  -o "$PROJECT/annotation/vep/HCC1395.PASS.VEP.vcf" \
  --vcf \
  --cache \
  --offline \
  --assembly GRCh38 \
  --dir_cache "$PROJECT/annotation/vep/cache" \
  --fasta "$PROJECT/reference/Homo_sapiens_assembly38.fasta" \
  --hgvs \
  --force_overwrite
