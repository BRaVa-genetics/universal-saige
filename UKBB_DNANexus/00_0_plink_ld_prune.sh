#!/bin/bash
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

set -e # Stop job if any command fails


# [OPTION] pop
# Population to use 
#   Options: 
#   - 'eur': Genetically European
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop=$1

# [OPTION] chrom
# Chromosome to run
readonly chrom=$2


## INPUT
readonly bfile="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c${chrom}_b0_v2"
readonly samples_w_superpop="/mnt/project/saige_pipeline/data/00_set_up/ukb_wes_450k.qced.sample_list_w_superpops.tsv"

## OUTPUT
readonly out="ukb_array.wes_450k_qc_pass_${pop}.pruned.chr${chrom}"

if [[ "${pop}" == "allpop" ]]; then
  plink \
    --bfile "${bfile}" \
    --indep-pairwise 50 5 0.05 \
    --out "${out}"

  # Extract set of pruned variants and export to bfile
  plink \
    --bfile "${bfile}" \
    --extract "${out}.prune.in" \
    --make-bed \
    --out "${out}"

elif [[ "${pop}" == "eur" ]]; then
  # Subset to a genetic ancestry
  # LD prune on European subset of samples passing WES 450k QC
  plink \
    --bfile "${bfile}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --indep-pairwise 50 5 0.05 \
    --out "${out}"

  # Extract set of pruned variants and export to bfile
  plink \
    --bfile "${bfile}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --extract "${out}.prune.in" \
    --make-bed \
    --out "${out}"
fi