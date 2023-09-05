#!/bin/bash

source ./run_container.sh

POSITIONAL_ARGS=()

SINGULARITY=false
generate_grm=false
generate_plink_for_vr=false

WD=$(pwd)

subset_variants(){
    echo "Subsetting genetic data for GRM / VR"

    # get list of files with format in dir:
    if [[ $GENETIC_DATA_FORMAT == "vcf" ]]; then

        FILES=$(ls ${GENETIC_DATA_DIR}/*vcf.gz)
        echo "files found: ${FILES}"

        for file in ${FILES}; do
            file_basename=$(basename "${file}")
            ./resources/plink --vcf "${file}" --make-bed --out "/tmp/${file_basename%.*}.plink"
        done

        ls /tmp/*.plink.bed | sed 's/\.bed$//g' > /tmp/merge_list.txt          

    elif [[ $GENETIC_DATA_FORMAT == "plink" ]]; then
        FILES=$(ls ${GENETIC_DATA_DIR}/*bed)

        for file in ${FILES}; do
            echo "${file%.*}" >> /tmp/plink_prefixes.txt
        done

        # Remove duplicate prefixes
        sort -u /tmp/plink_prefixes.txt > /tmp/merge_list.txt
    fi

    if [[ -n "$SAMPLEIDS" ]]; then
      ./resources/plink --merge-list /tmp/merge_list.txt --make-bed --out /tmp/merged --keep <(awk '{print $1, $1}' "$SAMPLEIDS")
    else
      ./resources/plink --merge-list /tmp/merge_list.txt --make-bed --out /tmp/merged 
    fi

}

generate_GRM(){
    echo "LD pruning file for GRM generation"

    numRandomMarkerforSparseKin=5000

    ./resources/plink \
        --bfile "/tmp/merged" \
        --indep-pairwise 50 5 0.05 \
        --out "/tmp/merged"

    # Extract set of pruned variants and export to bfile
    ./resources/plink \
        --bfile "/tmp/merged" \
        --extract "/tmp/merged.prune.in" \
        --make-bed \
        --out "${HOME}/${OUT}.plink_for_grm"

    cmd="createSparseGRM.R \
        --plinkFile="${HOME}/${OUT}.plink_for_grm" \
        --nThreads=$(nproc) \
        --outputPrefix="${HOME}/${OUT}" \
        --numRandomMarkerforSparseKin=$numRandomMarkerforSparseKin \
        --relatednessCutoff=0.05"

    variant_count=$(wc -l < "${HOME}/${OUT}.plink_for_grm.bim")
    if [[ $variant_count -ge $numRandomMarkerforSparseKin ]]; then
      run_container
    else
      echo "Error: ${variant_count} variants found in ${OUT}.plink_for_grm, which is less than the required ${numRandomMarkerforSparseKin} variants."
      exit 1
    fi
    
    echo "GRM generated!"
}

generate_plink_for_vr(){
    # get count of variants in merged plink file:

    ./resources/plink \
        --bfile "/tmp/merged" \
        --freq counts \
        --out "/tmp/merged"

    variants_lessthan_20_MAC=1000
    variants_greaterthan_20_MAC=1000

    cat <(
        tail -n +2 "/tmp/merged.frq.counts" \
        | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' \
        | shuf -n $variants_lessthan_20_MAC ) \
    <( \
        tail -n +2 "/tmp/merged.frq.counts" \
        | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' \
        | shuf -n $variants_greaterthan_20_MAC \
        ) > "/tmp/merged.markerid.list"

    actual_variants_lessthan_20_MAC=$(awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10)' "/tmp/merged.frq.counts" | wc -l)
    actual_variants_greaterthan_20_MAC=$(awk '$5 >= 20 && ($6-$5)>= 20' "/tmp/merged.frq.counts" | wc -l)

    if [[ $actual_variants_lessthan_20_MAC -gt $variants_lessthan_20_MAC ]]; then
        echo "Error: ${actual_variants_lessthan_20_MAC} variants (MAC<20) found - less than the required ${variants_lessthan_20_MAC} variants."
        exit 1
    elif [[ $actual_variants_greaterthan_20_MAC -gt $variants_greaterthan_20_MAC ]]; then
        echo "Error: ${actual_variants_greaterthan_20_MAC} variants (MAC>20) found - less than the required ${variants_greaterthan_20_MAC} variants."
        exit 1
    fi

    # Extract markers from the large PLINK file
    ./resources/plink \
        --bfile "/tmp/merged" \
        --extract "/tmp/merged.markerid.list" \
        --make-bed \
        --out "${OUT}.plink_for_var_ratio"
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
      generate_grm=true
      shift # past argument
      ;;
    --generate_plink_for_vr)
      generate_plink_for_vr=true
      shift # past argument
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
                -o,--outputPrefix: output prefix of the SAIGE step 0 output.
                -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
                --generate_GRM (default: false): generate GRM for the genetic data.
                --generate_plink_for_vr (default: false): generate plink file for vr.
                --sampleIDs: path to a file containing sampleIDs (as a single column) to be used to define the GRM.
                Note that if nothing is passed, then all of the samples in the plink/vcf files will be used.
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

# check if either generate_GRM or generate_plink_for_vr:
if [[ ${generate_grm} = false ]] && [[ ${generate_plink_for_vr} = false ]]; then
  echo "Error: either generate_GRM or generate_plink_for_vr must be set to true"
  exit 1
fi

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

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix. Check that this will not overwrite existing files."
  OUT="${PHENOCOL}"
fi

echo "OUT               = ${OUT}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "PLINK             = ${PLINK_WES}.{bim/bed/fam}"
echo "SAMPLEIDS         = ${SAMPLEIDS}"

check_container_env $SINGULARITY

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

subset_variants

if [[ ${generate_grm} = true ]]; then
  echo "generating GRM"
  generate_GRM
fi

if [[ ${generate_plink_for_vr} = true ]]; then
  echo "generating plink for vr"
  generate_plink_for_vr
fi
