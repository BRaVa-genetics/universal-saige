#!/bin/bash

source ./setup.sh
source ./check.sh

POSITIONAL_ARGS=()

SINGULARITY=false
SAMPLEIDCOL="IID"
OUT="out"
TRAITTYPE=""
PLINK=""
SPARSEGRM=""
SPARSEGRMID=""
PHENOFILE=""
PHENOCOL=""
COVARCOLLIST=""
CATEGCOVARCOLLIST=""
WD=$(pwd)

run_container () {
  if [[ ${SINGULARITY} = true ]]; then
    singularity exec \
      --env HOME=${WD} \
      --bind ${WD}/:$HOME/ \
      "saige-${saige_version}.sif" $cmd
  else
    docker run \
      -e HOME=${WD} \
      -v ${WD}/:$HOME/ \
      "wzhou88/saige:${saige_version}" $cmd
  fi
}

generate_GRM(){
  echo "Generating GRM..."
  ls -l
  pwd
  unameOut="$(uname -s)"
  case "${unameOut}" in
      Linux*)     machine=Linux;;
      Darwin*)    machine=Mac;;
      *)          machine="UNKNOWN:${unameOut}"
  esac
  echo ${machine}

  # download plink binary to resources:
  if [[ $machine == "Mac" ]]; then
    echo "Downloading OSX version of plink"
    wget -nc https://s3.amazonaws.com/plink1-assets/plink_mac_20230116.zip --no-check-certificate -P resources/
    unzip -o resources/plink_mac_20230116.zip -d resources/
  elif [[ $machine == "Linux" ]]; then
    echo "Downloading linux version of plink"
    wget -nc https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip --no-check-certificate -P resources/
    unzip -o resources/plink_linux_x86_64_20230116.zip -d resources/
  else
    echo "Operating system not compatible with the code"
  fi

  echo $SAMPLEIDS

  if [[ $GENETIC_DATA_FORMAT == "vcf" || $GENETIC_DATA_FORMAT == "bgen" ]]; then
    # get list of files with format in dir:
    FILES=$(ls ${GENETIC_DATA_DIR}/*${GENETIC_DATA_FORMAT}*)
    # iterate files (bash):
    if [[ $GENETIC_DATA_FORMAT == "vcf" ]]; then
        for FILE in $FILES
        do
            echo "Processing $f file..."
            # take action on each file. $f store current file name
            ./resources/plink \
                --vcf ${GENETIC_DATA_DIR}/${FILE} \
                --keep-fam ${SAMPLEIDS} \
                --indep-pairwise 50 5 0.05 \
                --out /tmp/${FILE}
        done
        for FILE in $FILES
        do
            ./resources/plink \
                --vcf "${GENOTYPE_PLINK}" \
                --keep-fam ${SAMPLEIDS} \
                --extract "${OUT}.prune.in" \
                --make-bed \
                --out "${OUT}"
            
        done
    elif [[ $GENETIC_DATA_FORMAT == "bgen" ]]; then
        for FILE in $FILES
        do
            echo "Processing $f file..."
            # take action on each file. $f store current file name
            ./resources/plink \
                --bgen ${GENETIC_DATA_DIR}/${FILE} \
                --keep-fam ${SAMPLEIDS} \
                --indep-pairwise 50 5 0.05 \
                --out /tmp/${FILE}
        done
        for FILE in $FILES
        do
            ./resources/plink \
                --bgen "${GENOTYPE_PLINK}" \
                --keep-fam ${SAMPLEIDS} \
                --extract /tmp/${FILE}.prune.in \
                --make-bed \
                --out "${OUT}"
        done
    fi
  elif [[ $GENETIC_DATA_FORMAT == "plink" ]]; then
    FILES=$(ls ${GENETIC_DATA_DIR}/*.bed)

    for FILE in $FILES
    do
        echo "Processing $f file..."
        # take action on each file. $f store current file name
        ./resources/plink \
            --plink ${GENETIC_DATA_DIR}/${FILE%.bed} \
            --keep-fam ${SAMPLEIDS} \
            --indep-pairwise 50 5 0.05 \
            --out /tmp/${FILE}
    done
    for FILE in $FILES
    do
        ./resources/plink \
            --plink ${GENETIC_DATA_DIR}/${FILE%.bed} \
            --keep-fam ${SAMPLEIDS} \
            --extract /tmp/${FILE}.prune.in \
            --make-bed \
            --out "${OUT}"
    done

  fi

  # Extract set of pruned variants and export to bfile
  
  cmd="createSparseGRM.R \
    --plinkFile="${HOME}/${OUT}" \
    --nThreads=$(nproc) \
    --outputPrefix="${HOME}/${OUT}" \
    --numRandomMarkerforSparseKin=5000 \
    --relatednessCutoff=0.05"

  run_container
 
  echo "GRM generated!"
 
  SPARSEGRM="${OUT}.sparseGRM.mtx"
  SPARSEGRMID="${OUT}.sparseGRM.mtx.sampleIDs.txt"

}

generate_plink_for_vr(){
  echo "Generating plink file for vr..."

  wget -nc https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20230325.zip -P resources/
  unzip -o resources/plink2_linux_x86_64_20230325.zip -d resources/

  #1. Calculate allele counts for each marker in the large PLINK file with hard called genotypes

  ./resources/plink2 \
    --keep <(awk '{ print $1,$1 }' ${SAMPLEIDS})  \
    --bfile "${GENOTYPE_PLINK}" \
    --freq counts \
    --out "${OUT}"

  #2. Randomly extract IDs for markers falling in the two MAC categories:
  # * 1,000 markers with 10 <= MAC < 20
  # * 1,000 markers with MAC >= 20

  cat <(
    tail -n +2 "${OUT}.acount" \
    | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' \
    | shuf -n 1000 ) \
  <( \
    tail -n +2 "${OUT}.acount" \
    | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' \
    | shuf -n 1000 \
    ) > "${OUT}.markerid.list"

  # Make sure to still subset to Europeans
  ./resources/plink2 \
    --bfile "${GENOTYPE_PLINK}" \
    --keep <(awk '{ print $1,$1 }' ${SAMPLEIDS}) \
    --extract "${OUT}.markerid.list" \
    --make-bed \
    --out "${OUT}"

}

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      shift # past argument
      shift # past value
      ;;
    -p|--geneticDataDirectory)
      GENETIC_DATA_DIR="$2"
      shift # past argument
      shift # past value
      ;;
    --geneticDataFormat)
      GENETIC_DATA_FORMAT="$2"
      shift # past argument
      shift # past value
      ;;
    --generate_GRM)
      shift # past argument
      shift # past value
      ;;
    --generate_plink_for_vr)
      shift # past argument
      shift # past value
      ;;
    --sampleIDs)
      SAMPLEIDS="$2" 
      shift
      shift
      ;; 
    -h|--help)
      echo "usage: 00_step0_VR_and_GRM.sh
  required:
    --geneticDataDirectory: directory containing the genetic data (genotype/WES/WGS data in the format plink/vcf/bgen)
    --geneticDataFormat: format of the genetic data {plink,vcf,bgen}.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 1 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    --sampleIDCol (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.
    --generate_GRM (default: false): generate GRM for the genetic data.
    --generate_plink_for_vr (default: false): generate plink file for vr.
      "
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

# Checks

# check if genetic data format in vcf, bgen or plink:
if [[ ${GENETIC_DATA_FORMAT} != "vcf" ]] && [[ ${GENETIC_DATA_FORMAT} != "bgen" ]] && [[ ${GENETIC_DATA_FORMAT} != "plink" ]]; then
  echo "geneticDataFormat must be in {vcf,bgen,plink}"
  exit 1
fi

# check if genetic data directory exists and if files with the correct extension are present:
if [[ ! -d ${GENETIC_DATA_DIR} ]]; then
  echo "geneticDataDirectory does not exist"
  exit 1
fi

if [[ ${SAMPLEIDS} != "" ]]; then
  SAMPLEIDS=${HOME}/$SAMPLEIDS
fi

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix to ${PHENOCOL}. Check that this will not overwrite existing files."
  OUT="${PHENOCOL}"
fi

echo "OUT               = ${OUT}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "PLINK             = ${PLINK_WES}.{bim/bed/fam}"
echo "SAMPLEIDS         = ${SAMPLEIDS}"

check_container_env $SINGULARITY
 
if [[ ${SINGULARITY} = true && ! $( test -f "saige-${saige_version}.sif" ) ]]; then
  singularity pull "saige-${saige_version}.sif" "docker://wzhou88/saige:${saige_version}"
elif [[ ${SINGULARITY} = false ]]; then
  docker pull wzhou88/saige:${saige_version}
fi

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

# Get number of threads
n_threads=$(( $(nproc --all) - 1 ))

if [[ ${GENERATE_GRM} = true ]]; then
  generate_GRM
elif [[ ${GENERATE_PLINK_FOR_VR} = true ]]; then
  generate_plink_for_vr
fi