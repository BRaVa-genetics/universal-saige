#!/bin/bash

# check_exome_data {
#   if [[ $bimFile == "" && $bedFile == "" && $famFile == "" && $vcfFile == "" ]]; then
#     echo "Error: either plink files or VCF file is required"
#     exit(1)
#   if [[ $bimFile != "" && $bedFile != "" && $famFile != "" && $vcfFile != "" ]]; then
#     echo "Error: either plink files or VCF file is required"
#     exit(1)
#   fi
# }

function check_container_env() {
  # Check if singularity or docker is installed
  SINGULARITY="$1"
  if [[ ${SINGULARITY} = true ]]; then
    if [[ $(which singularity) == "" ]]; then
      echo "Error: singularity is not installed"
      exit 1
    fi
  else
    if [[ $(which docker) == "" ]]; then
        echo "Error: docker is not installed"
        exit 1
    fi
  fi
  echo "singularity or docker found"
}

# check_sparse_grm_data {
#     if [[ $genotypeFileBim == "" && $genotypeFileBed == "" && $genotypeFileFam == "" && $exomeFileBim == "" && $exomeFileBed == "" && $exomeFileFam == "" ]];
#         echo "Error: either genotype plink files or exome plink files are required for sparse GRM
#         exit(1)
#     }
# }