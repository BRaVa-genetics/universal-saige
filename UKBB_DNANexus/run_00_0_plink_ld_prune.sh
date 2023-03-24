#!/bin/bash

set -u # throws error if variables are undefined

export project=`dx pwd`
WD="/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_gwas"

readonly script="00_0_plink_ld_prune.sh"

script_local="${WD}/bash/${script}"
script_dnax="/saige_pipeline/scripts/${script}"

source "/Users/nbaya/gms/lindgren/ukb_wes/ukb_wes_450k_qc/bash/dnax_utils.sh"
upload_file "${script_local}" "${script_dnax}"

# [OPTION] pop
# Population to use
# - 'eur': Subset to genetic European
# - 'allpop': Do not subset to a single ancestry (i.e. use all individuals who pass QC)
pop="allpop"


run_job() {

	chrom=$1

	if [ "${chrom}" -eq 23 ]; then
		chrom="X";
	fi

	dx run swiss-army-knife \
		-iin="/saige_pipeline/scripts/${script}" \
		-icmd="bash ${script} ${pop} ${chrom}" \
		--name="plink_ld_prune_${pop}_c${chrom}" \
		--instance-type="mem1_ssd1_v2_x2" \
		--destination="/saige_pipeline/data/00_set_up" \
		--priority="low" \
		--brief \
		-y
}

max_tasks=8
i=0
(
for chrom in {1..23}; do 
   ((i=i%max_tasks)); ((i++==0)) && wait
   run_job "${chrom}" & 
done
)

# run_job "$1"