run_container_cmd(CMD, is_singularity){
    if (is_singularity){
        singularity pull saige.sif docker://wzhou88/saige:1.1.6.3
        singularity run saige.sif $CMD

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
                        --includeNonautoMarkersforVarRatio=TRUE
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
                    --minMAC=20 \
                    --is_fastTest TRUE \
    """

    run_container_cmd(CMD, is_singularity)
}

main(){
    preprocess_step1()
    step1()
    preprocess_step2()
    step2()
}

main()