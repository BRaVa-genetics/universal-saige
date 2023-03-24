#!/bin/bash
# 01_saige_fit_null 0.0.1

main() {

    # For debugging
    set -exo pipefail

    ## Set up directories
    WD=$( pwd )
    # mkdir -p "plink_files"
    mkdir -p out/{model_file,variance_ratios}

    dx-mount-all-inputs --except pheno_file
    
    # Download phenotype file
    # dx download "$pheno_file" -o "pheno_file.tsv.gz"
    # zcat "pheno_file.tsv.gz" > pheno_file
    dx download "$pheno_file" -o "pheno_file"

    # Download SAIGE docker image from ukbb_meta/docker/
    dx download file-GK53YGjJg8JX4yqg925zY7x5
    docker load --input saige-1.1.6.3.tar.gz

    # Get number of threads
    n_threads=$(( $(nproc --all) - 1 ))

    # Get inverse-normalize flag if trait_type=="quantitative"
    if [ ${trait_type} == "quantitative" ]; then
      trait_flags="--traitType=${trait_type}   --invNormalize=TRUE"
    else
      trait_flags="--traitType=${trait_type}"
    fi

    ## Run script
    docker run \
      -e HOME=${WD}  \
      -e pheno_col="${pheno_col}" \
      -e trait_type="${trait_type}" \
      -e inv_normalize_flag=${inv_normalize_flag} \
      -e output_prefix="${output_prefix}"  \
      -e n_threads="${n_threads}" \
      -v ${WD}/:$HOME/ \
      wzhou88/saige:1.1.6.3 step1_fitNULLGLMM.R  \
        --bedFile ${HOME}/in/plink_for_vr_bed/* \
        --bimFile ${HOME}/in/plink_for_vr_bim/* \
        --famFile ${HOME}/in/plink_for_vr_fam/* \
        --sparseGRMFile ${HOME}/in/sparse_grm/* \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm_samples/*  \
        --useSparseGRMtoFitNULL=TRUE  \
        --phenoFile ${HOME}/pheno_file \
        --skipVarianceRatioEstimation FALSE \
        --phenoCol "${pheno_col}" \
        --covarColList "${covar_col_list}" \
        --qCovarColList="${qcovar_col_list}"  \
        --sampleIDColinphenoFile="IID" \
        ${trait_flags} \
        --outputPrefix="${HOME}/${output_prefix}" \
        --IsOverwriteVarianceRatioFile=TRUE \
        --nThreads=${n_threads}

    mv *.rda out/model_file/
    mv *.varianceRatio.txt out/variance_ratios/

    dx-upload-all-outputs
}