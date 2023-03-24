#!/bin/bash
# 02_saige_variant_test 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    mkdir plink_files
    mkdir -p out/{output,log}_file
    
    dx-mount-all-inputs --except bim

    dx download "$bim" -o "bfile.bim"
    sed -i 's/^chr//g' "${WD}/bfile.bim"

    # Download saige-1.1.6.3.tar.gz docker image from ukbb_meta/docker/
    dx download file-GK53YGjJg8JX4yqg925zY7x5
    docker load --input saige-1.1.6.3.tar.gz

    ## Run script
    docker run \
      -e HOME=${WD}  \
      -e output_prefix="${output_prefix}"  \
      -v ${WD}/:$HOME/ \
      wzhou88/saige:1.1.6.3 \
      step2_SPAtests.R  \
        --bedFile ${HOME}/in/bed/* \
        --bimFile ${HOME}/bfile.bim \
        --famFile ${HOME}/in/fam/* \
        --AlleleOrder=ref-first \
        --minMAF=0  \
        --minMAC=20  \
        --GMMATmodelFile ${HOME}/in/model_file/*  \
        --varianceRatioFile ${HOME}/in/variance_ratios/*  \
        --sparseGRMFile ${HOME}/in/sparse_grm/*  \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm_samples/* \
        --LOCO=FALSE  \
        --is_Firth_beta=TRUE  \
        --pCutoffforFirth=0.1 \
        --is_output_moreDetails=TRUE  \
        --is_fastTest=TRUE  \
        --SAIGEOutputFile=${HOME}/${output_prefix}.tsv 2>&1 | tee ${HOME}/${output_prefix}.log

    gzip "${output_prefix}.tsv" 
    mv "${output_prefix}.tsv.gz" "out/output_file/" 
    mv "${output_prefix}.log" "out/log_file/"

    dx-upload-all-outputs
}