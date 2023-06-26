#!/bin/bash

source ./run_container.sh
source ./check_pheno.sh

POSITIONAL_ARGS=()

SINGULARITY=false
OUT="out"
TESTTYPE=""
PLINK=""
VCF=""
MODELFILE=""
VARIANCERATIO=""
SPARSEGRM=""
SPARSEGRMID=""
GROUPFILE=""

while [[ $# -gt 0 ]]; do
  case $1 in
    --ancestry)
      ANC="$2"
      shift # past argument
      shift # past value
      ;;
    --sex)
      SEX="$2"
      shift # past argument
      shift # past value
      ;;
    --dataset)
      DATASET="$2"
      shift # past argument
      shift # past value
      ;;
    --lastName)
      LAST_NAME="$2"
      shift # past argument
      shift # past value
      ;;
    --freezeNumber)
      FREEZE_NUMBER="$2"
      shift # past argument
      shift # past value
      ;;  
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    --chr)
      CHR="$2"
      shift # past argument
      shift # past value
      ;;
    --testType)
      TESTTYPE="$2"
      if ! ( [[ ${TESTTYPE} = "variant" ]] || [[ ${TESTTYPE} = "group" ]] ); then
        echo "Test type is not in {variant,group}"
        exit 1
      fi
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      shift # past argument
      shift # past value
      ;;
    -p|--plink)
      PLINK="$2"
      shift # past argument
      shift # past value
      ;;
    --vcf)
      VCF="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--modelFile)
      MODELFILE="$2"
      shift # past argument
      shift # past value
      ;;
    -v|--varianceRatio)
      VARIANCERATIO="$2"
      shift # past argument
      shift # past value
      ;;
    -g|--groupFile)
      GROUPFILE="$2"
      shift # past argument
      shift # past value
      ;; 
    --annotations)
      ANNOTATIONS="$2"
      shift # past argument
      shift # past value
      ;;
    --subSampleFile)
      SUBSAMPLES="$2"
      shift # past argument
      shift # past value
      ;;
    --sparseGRM)
      SPARSEGRM="$2"
      shift # past argument
      shift # past value
      ;;
    --sparseGRMID)
      SPARSEGRMID="$2"
      shift # past argument
      shift # past value
      ;;
    --phenotype)
      PHENOCOL="$2"
      shift # past argument
      shift # past value
      ;;
    --phenoFile)
      PHENOFILE="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "usage: 02_step2_SPAtests_variant_and_gene.sh
  required:
    --testType: type of test {variant,group}.
    -p,--plink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --vcf vcf exome file. If a plink exome file is not available then this vcf file will be used. These must be present in the working directory at ./in/vcf/
    --modelFile: filename of the model file output from step 1. This must be in relation to the working directory.
    --varianceRatio: filename of the varianceRatio file output from step 1. This must be in relation to the working directory.
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
	--chr: chromosome to test.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 2 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -g,--groupFile: required if group test is selected. Filename of the annotation file used for group tests. This must be in relation to the working directory.
    --annotations: required if group test is selected. comma seperated list of annotations to test found in groupfile.
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

if is_valid_r_var "$PHENOCOL"; then
    echo "The variable name '$PHENOCOL' is not valid for an R variable."
fi

# Checks
if [[ ${TESTTYPE} == "" ]]; then
  echo "Test type not set"
  exit 1
fi

if [[ ${PLINK} == "" ]] && [[ ${VCF} == "" ]]; then
  echo "plink files plink.{bim,bed,fam} and vcf not set"
  exit 1
fi

if [[ ${SPARSEGRM} == "" ]]; then
  echo "sparse GRM .mtx file not set"
  exit 1
fi

if [[ ${SPARSEGRMID} == "" ]]; then
  echo "sparse GRM ID file not set"
  exit 1
fi

if [[ ${MODELFILE} == "" ]]; then
  echo "model file not set"
  exit 1
fi

if [[ ${VARIANCERATIO} == "" ]]; then
  echo "variance ration file not set"
  exit 1
fi

if [[ $GROUPFILE == "" ]] && [[ ${TESTTYPE} == "group" ]]; then
  echo "attempting to run group tests without an annotation file"
  exit 1
fi

if [[ $ANNOTATIONS == "" ]] && [[ ${TESTTYPE} == "group" ]]; then
  echo "attempting to run group tests without selected annotations"
  exit 1
fi

if [[ $SUBSAMPLES != "" ]]; then
  SUBSAMPLES="${HOME}/${SUBSAMPLES}"
fi

# Get output file name

# Get column numbers for 'sex' and the specified phenotype
sex_col_num=$(head -n 1 $pheno_file | tr ' ' '\n' | grep -n -w 'sex' | cut -d: -f1)
pheno_col_num=$(head -n 1 $pheno_file | tr ' ' '\n' | grep -n -w $PHENOCOL | cut -d: -f1)

# Subtract 1 from column numbers because awk counts from 1 but grep from 0
sex_col_num=$((sex_col_num - 1))
pheno_col_num=$((pheno_col_num - 1))

if [[ $trait_type = "continuous" ]]; then
  # Count valid phenotype values for males (M), females (F), and both sexes
  if [[ $SEX == "BOTH" ]]; then
    N=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '(($sex_col=="M") || ($sex_col=="F")) && ($pheno_col!="NA") {count++} END {print count}' $pheno_file)
  elif [[ $SEX == "F" ]]; then
    N=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '($sex_col=="F") && ($pheno_col!="NA") {count++} END {print count}' $pheno_file)
  elif [[ $SEX == "M" ]]; then
    N=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '($sex_col=="M") && ($pheno_col!="NA") {count++} END {print count}' $pheno_file)
  fi

  OUT="${DATASET}.${LAST_NAME}.chr${CHR}_${phenoCol}.${FREEZE_NUMBER}.${SEX}.${ANC}.${N}.SAIGE.$(date '+%Y%m%d')"

elif [[ $trait_type = "binary" ]]; then
  if [[ $SEX == "BOTH" ]]; then
    N_case=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '(($sex_col=="M") || ($sex_col=="F")) && ($pheno_col=="1") {count++} END {print count}' $pheno_file)
    N_control=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '(($sex_col=="M") || ($sex_col=="F")) && ($pheno_col=="0") {count++} END {print count}' $pheno_file)
  elif [[ $SEX == "F" ]]; then
    N_case=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '($sex_col=="F") && ($pheno_col=="1") {count++} END {print count}' $pheno_file)
    N_control=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '($sex_col=="F") && ($pheno_col=="0") {count++} END {print count}' $pheno_file)
  elif [[ $SEX == "M" ]]; then
    N_case=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '($sex_col=="M") && ($pheno_col=="1") {count++} END {print count}' $pheno_file)
    N_control=$(awk -v sex_col=$sex_col_num -v pheno_col=$pheno_col_num '($sex_col=="M") && ($pheno_col=="0") {count++} END {print count}' $pheno_file)
  fi

  OUT="${DATASET}.${LAST_NAME}.chr${CHR}_${phenoCol}.${FREEZE_NUMBER}.${SEX}.${ANC}.${N_case}.${N_control}.SAIGE.$(date '+%Y%m%d')"
fi

echo "OUT               = ${OUT}"
echo "TESTTYPE          = ${TESTTYPE}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "PLINK             = ${PLINK}.{bim/bed/fam}"
echo "MODELFILE         = ${MODELFILE}"
echo "VARIANCERATIO     = ${VARIANCERATIO}"
echo "GROUPFILE         = ${GROUPFILE}"
echo "ANNOTATIONS"      = ${ANNOTATIONS}
echo "SPARSEGRM         = ${SPARSEGRM}"
echo "SPARSEGRMID       = ${SPARSEGRMID}"

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

# Get number of threads
n_threads=$(( $(nproc --all) - 1 ))

## Set up directories
WD=$( pwd )


if [[ "$TESTTYPE" = "variant" ]]; then
  echo "variant testing"
  min_mac="20"
  GROUPFILE=""
else
  echo "gene testing"
  min_mac="0.5"
  GROUPFILE="${HOME}/${GROUPFILE}"
fi

if [[ ${PLINK} != "" ]]; then
  PLINK="${HOME}/${PLINK}"
  BED=${PLINK}".bed"
  BIM=${PLINK}".bim"
  FAM=${PLINK}".fam"
  VCF=""
elif [[ ${VCF} != "" ]]; then 
  BED=""
  BIM=""
  FAM="" 
  VCF="${VCF}"
else
  echo "No plink or vcf found!"
  exit 1
fi

cmd="step2_SPAtests.R \
        --bedFile=$BED \
        --bimFile=$BIM \
        --famFile=$FAM \
        --groupFile=$GROUPFILE \
        --annotation_in_groupTest=$ANNOTATIONS \
        --vcfFile=${VCF} \
        --vcfField="DS" \
        --chrom="$CHR" \
        --minMAF=0 \
        --minMAC=${min_mac} \
        --GMMATmodelFile=${HOME}/${MODELFILE} \
        --varianceRatioFile=${HOME}/${VARIANCERATIO} \
        --sparseGRMFile=${HOME}/${SPARSEGRM} \
        --sparseGRMSampleIDFile=${HOME}/${SPARSEGRMID} \
        --subSampleFile=${SUBSAMPLES} \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.1 \
        --is_output_moreDetails=TRUE \
        --is_fastTest=TRUE \
        --is_output_markerList_in_groupTest=TRUE \
        --is_single_in_groupTest=TRUE \
        --maxMAF_in_groupTest=0.0001,0.001,0.01 \
        --SAIGEOutputFile=${HOME}/${OUT}.txt
    "

run_container
