#!/bin/bash
#
# Create PLINK dataset with 2000 randomly selected markers, using the hard-called genotypes.
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

## OPTIONS

# Population to use 
#   Options: 
#   - 'eur': Genetically European
#   - 'allpop': Individuals from all populations who pass QC (i.e. no population filter)
readonly pop=$1 

# Chromosome
readonly chrom=$2


## INPUT
readonly array_bfile="/mnt/project/Bulk/Genotype Results/Genotype calls/ukb22418_c${chrom}_b0_v2"
readonly samples_w_superpop="/mnt/project/saige_pipeline/data/00_set_up/ukb_wes_450k.qced.sample_list_w_superpops.tsv"


## OUTPUT
readonly out="ukb_array.wes_450k_qc_pass_${pop}.for_vr.chr${chrom}"


#1. Calculate allele counts for each marker in the large PLINK file with hard called genotypes
if [[ "${pop}" == "allpop" ]]; then 
  plink2 \
    --bfile "${array_bfile}" \
    --freq counts \
    --out "${out}"

elif [[ "${pop}" == "eur" ]]; then 
  plink2 \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --bfile "${array_bfile}" \
    --freq counts \
    --out "${out}"
fi

#2. Randomly extract IDs for markers falling in the two MAC categories:
# * 1,000 markers with 10 <= MAC < 20
# * 1,000 markers with MAC >= 20
cat <(
  tail -n +2 "${out}.acount" \
  | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' \
  | shuf -n 1000 ) \
<( \
  tail -n +2 "${out}.acount" \
  | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' \
  | shuf -n 1000 \
  ) > "${out}.markerid.list"


#3. Extract markers from the large PLINK file
if [[ "${pop}" == "allpop" ]]; then 
  plink2 \
    --bfile "${array_bfile}" \
    --extract "${out}.markerid.list" \
    --make-bed \
    --out "${out}"

elif [[ "${pop}" == "eur" ]]; then 
  # Make sure to still subset to Europeans
  plink2 \
    --bfile "${array_bfile}" \
    --keep <( awk '{ if ($2=="EUR") print $1,$1 }' "${samples_w_superpop}" ) \
    --extract "${out}.markerid.list" \
    --make-bed \
    --out "${out}"
fi
