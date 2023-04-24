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

subset_variants(){
    echo "Subsetting genetic data for GRM / VR"

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

    # get list of files with format in dir:
    if [[ $GENETIC_DATA_FORMAT == "vcf" ]]; then

        FILES=$(ls ${GENETIC_DATA_DIR}/*vcf)
        echo "files found: ${FILES}"

        for file in ${FILES}; do
            ./resources/plink --vcf "${file}" --make-bed --out "/tmp/${file%.vcf}.plink"
        done

        ls /tmp/*.plink.bed | sed 's/\.bed$//g' > /tmp/merge_list.txt

        ./resources/plink --merge-list /tmp/merge_list.txt --make-bed --out /tmp/merged

    elif [[ $GENETIC_DATA_FORMAT == "plink" ]]; then
        FILES=$(ls ${GENETIC_DATA_DIR}/*bed)

        for file in ${FILES}; do
            echo "${file%.*}" >> /tmp/plink_prefixes.txt
        done

        # Remove duplicate prefixes
        sort -u /tmp/plink_prefixes.txt > /tmp/unique_plink_prefixes.txt

        # Merge the PLINK files
        ./resources/plink --merge-list /tmp/unique_plink_prefixes.txt --make-bed --out /tmp/merged

    fi
}

generate_GRM(){
    echo "LD pruning file for GRM generation"

    ./resources/plink \
        --bfile "/tmp/merged" \
        --indep-pairwise 50 5 0.05 \
        --out "/tmp/${out}"

    # Extract set of pruned variants and export to bfile
    ./resources/plink \
        --bfile "/tmp/merged" \
        --extract "/tmp/${out}.prune.in" \
        --make-bed \
        --out "${HOME}/out/plink_for_grm"

    cmd="createSparseGRM.R \
        --plinkFile="${HOME}/out/plink_for_grm" \
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
    # get count of variants in merged plink file:

    ./resources/plink \
        --bfile "/tmp/merged" \
        --freq counts \
        --out "/tmp/${out}"

    variants_lessthan_20_MAC=1000
    variants_greaterthan_20_MAC=1000

    cat <(
        tail -n +2 "/tmp/${out}.frq.counts" \
        | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' \
        | shuf -n $variants_lessthan_20_MAC ) \
    <( \
        tail -n +2 "/tmp/${out}.frq.counts" \
        | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' \
        | shuf -n $variants_greaterthan_20_MAC \
        ) > "/tmp/${out}.markerid.list"

    actual_variants_lessthan_20_MAC=$(awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10)' "/tmp/${out}.frq.counts" | wc -l)
    actual_variants_greaterthan_20_MAC=$(awk '$5 >= 20 && ($6-$5)>= 20' "/tmp/${out}.frq.counts" | wc -l)

    if [[ $actual_variants_lessthan_20_MAC -ne $variants_lessthan_20_MAC ]]; then
        echo "Error: ${actual_variants_lessthan_20_MAC} variants (MAC<20) found - less than the required ${variants_lessthan_20_MAC} variants."
        exit 1
    elif [[ $actual_variants_greaterthan_20_MAC -ne $variants_greaterthan_20_MAC ]]; then
        echo "Error: ${actual_variants_greaterthan_20_MAC} variants (MAC>20) found - less than the required ${variants_greaterthan_20_MAC} variants."
        exit 1
    fi

    # Extract markers from the large PLINK file
    ./resources/plink \
        --bfile "/tmp/merged" \
        --extract "/tmp/${out}.markerid.list" \
        --make-bed \
        --out "./in/variants_subset"
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
    --geneticDataType)
      GENETIC_DATA_TYPE="$2"
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
                --geneticDataFormat: format of the genetic data {plink,vcf}.
                --geneticDataType: type of the genetic data {WES,WGS,genotype}.
            optional:
                -o,--outputPrefix:  output prefix of the SAIGE step 0 output.
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

# check if genetic data format in vcf or plink:
if [[ ${GENETIC_DATA_FORMAT} != "vcf" ]] && [[ ${GENETIC_DATA_FORMAT} != "plink" ]]; then
  echo "geneticDataFormat must be in {vcf,plink}"
  exit 1
fi

if [[ ${GENETIC_DATA_TYPE} != "WES" ]] && [[ ${GENETIC_DATA_TYPE} != "WGS" ]] && [[ ${GENETIC_DATA_TYPE} != "genotype" ]]; then
  echo "geneticDataType must be in {WES,WGS,genotype}"
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

subset_variants

if [[ ${GENERATE_GRM} = true ]]; then
  generate_GRM
elif [[ ${GENERATE_PLINK_FOR_VR} = true ]]; then
  generate_plink_for_vr
fi