check_exome_data(bimFile, bedFile, famFile, vcfFile){
    if [[ $bimFile == "" && $bedFile == "" && $famFile == "" && $vcfFile == "" ]];
    then
        echo "Error: either plink files or VCF file is required"
        exit(1)
    if [[ $bimFile != "" && $bedFile != "" && $famFile != "" && $vcfFile != "" ]];
    then
        echo "Error: either plink files or VCF file is required"
        exit(1)
    fi
}

check_container_env(is_singularity){
    # check if singularity or docker is installed:
    if [ is_singularity ]
    then
        if [ system("which singularity") != 0 ]
        then
            echo "Error: singularity is not installed"
            exit(1)
        fi
    elif [ !is_singularity ]
        if [system("which docker") != 0]
        then
            echo "Error: docker is not installed"
            exit(1)
        fi
    }
}

check_sparse_grm_data(genotypeFileBim, genotypeFileBed, genotypeFileFam, exomeFileBim, exomeFileBed, exomeFileFam){
    if [[ $genotypeFileBim == "" && $genotypeFileBed == "" && $genotypeFileFam == "" && $exomeFileBim == "" && $exomeFileBed == "" && $exomeFileFam == "" ]];
        echo "Error: either genotype plink files or exome plink files are required for sparse GRM"
        exit(1)
    }
}