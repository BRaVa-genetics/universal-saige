# Import functions from other scripts:
. check.sh
. setup.sh

run_container_cmd(CMD, is_singularity){
    if (is_singularity){
        # check if saige.sif exists:
        if (!file.exists("saige.sif"))
            singularity pull saige-1.1.6.3.sif docker://wzhou88/saige:1.1.6.3

        singularity run saige-1.1.6.3 $CMD

    } else if (!is_singularity){
        docker pull wzhou88/saige:1.1.6.3
        docker run -v $PWD:/data wzhou88/saige:1.1.6.3 $CMD
    }
}

preprocess_step1(){

}

preprocess_step2(){

}

step1(sparseGRMFile, sparseGRMSampleIDFile, phenoFile, phenoCol, sampleIDColinphenoFile, covarColList, qCovarColList, bedFile, bimFile, famFile, outputPrefix, traitType, SampleIDIncludeFile, is_singularity){
    CMD = 
    """
    step1_fitNULLGLMM.R --sparseGRMFile $sparseGRMFile \
                        --sparseGRMSampleIDFile $sparseGRMSampleIDFile \
                        --useSparseGRMtoFitNULL TRUE \
                        --phenoFile $phenoFile \
                        --phenoCol="$phenoCol" \
                        --sampleIDColinphenoFile=$sampleIDColinphenoFile \
                        --covarColList=$covarColList \
                        --qCovarColList="$qCovarColList" \
                        --bedFile $bedFile \
                        --bimFile $bimFile \
                        --famFile $famFile \
                        --outputPrefix $outputPrefix \
                        --nThreads=$(nproc) \
                        --traitType="$trait_type" \
                        --invNormalize=TRUE \
                        --skipVarianceRatioEstimation=FALSE \
                        --SampleIDIncludeFile $SampleIDIncludeFile \
                        --includeNonautoMarkersforVarRatio=TRUE \
                        --isCateVarianceRatio=TRUE
    """

    run_container_cmd(CMD, is_singularity)
}

step2(bedFile, bimFile, famFile, GMMATmodelFile, varianceRatioFile, sparseGRMFile, sparseGRMSampleIDFile, SAIGEOutputFile, chrom, subSampleFile, is_singularity){
    CMD = 
    """
    step2_SPAtests.R --bedFile $bedFile \
                    --bimFile $bimFile \
                    --famFile $famFile \
                    --GMMATmodelFile $GMMATmodelFile \
                    --varianceRatioFile $varianceRatioFile \
                    --sparseGRMFile $sparseGRMFile \
                    --sparseGRMSampleIDFile $sparseGRMSampleIDFile \
                    --SAIGEOutputFile $SAIGEOutputFile \
                    --LOCO FALSE \
                    --is_Firth_beta TRUE \
                    --pCutoffforFirth=0.05 \
                    --chrom "$chrom" \
                    --subSampleFile $subSampleFile \
                    --minMAF=0 \
                    --minMAC=0.5 \
                    --is_fastTest TRUE \
    """

    run_container_cmd(CMD, is_singularity)
}

phenoCol, GMMATmodelFile, SAIGEOutputFile, subSampleFile, chrom, bimFile, sparseGRMFile, outputPrefix, qCovarColList, sparseGRMSampleIDFile, sampleIDColinphenoFile, covarColList, varianceRatioFile, SampleIDIncludeFile, phenoFile, bedFile, famFile, is_singularity){

main(phenoCol, subSampleFile, covarColList, qCovarColList, chrom, # Phenotype+covariate variables
     GMMATmodelFile, SAIGEOutputFile, sparseGRMFile, outputPrefix, # Where to write files
     exomeFileBim, exomeFileBed, exomeFileFam, vcfFile, # Exome data
     genotypeFileBim, genotypeFileBed, genotypeFileFam)  # Genotype data
{    
    # check if either exome plink files or VCF file are provided but not both:
    check_exome_data(exomeFileBim, exomeFileBed, exomeFileFam vcfFile)
    # check if either genotype plink files or exome plink files are provided for sparse GRM:
    check_sparse_grm_data(genotypeFileBim, genotypeFileBed, genotypeFileFam, exomeFileBim, exomeFileBed, exomeFileFam)
    # check if singularity or docker is installed:
    check_container_env(is_singularity)

    preprocess_step1()
    step1(sparseGRMFile, sparseGRMSampleIDFile, phenoFile, phenoCol, sampleIDColinphenoFile, covarColList, qCovarColList, bedFile, bimFile, famFile, outputPrefix, traitType, SampleIDIncludeFile, is_singularity)
    preprocess_step2()
    step2(bedFile, bimFile, famFile, GMMATmodelFile, varianceRatioFile, sparseGRMFile, sparseGRMSampleIDFile, SAIGEOutputFile, chrom, subSampleFile, is_singularity)
}

main(phenoCol, GMMATmodelFile, SAIGEOutputFile, subSampleFile, chrom, bimFile, sparseGRMFile, outputPrefix, qCovarColList, traitType, sparseGRMSampleIDFile, sampleIDColinphenoFile, covarColList, varianceRatioFile, SampleIDIncludeFile, phenoFile, bedFile, famFile, is_singularity)