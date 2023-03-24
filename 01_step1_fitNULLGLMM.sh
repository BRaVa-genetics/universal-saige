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

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      SINGULARITY=true
      shift # past argument
      shift # past value
      ;;
    -t|--traitType)
      TRAITTYPE="$2"
      if ! ( [[ ${TRAITTYPE} = "quantitative" ]] || [[ ${TRAITTYPE} = "binary" ]] ); then
        echo "Trait type is not in {quantitative,binary}"
        exit 1
      fi
      shift # past argument
      shift # past value
      ;;
    --plink)
      PLINK="$2"
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
    --phenoFile)
      PHENOFILE="$2"
      shift # past argument
      shift # past value
      ;;
    --phenoCol)
      PHENOCOL="$2"
      shift # past argument
      shift # past value
      ;;
    -c|--covarColList)
      COVARCOLLIST="$2"
      shift # past argument
      shift # past value
      ;;
    --categCovarColList)
      CATEGCOVARCOLLIST="$2"
      shift # past argument
      shift # past value
      ;;
    -i|--sampleIDCol)
      SAMPLEIDCOL="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "usage: 01_step1_fitNULLGLMM.sh
  required:
    -t,--traitType: type of the trait {quantitative,binary}.
    -p,--plink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
    --phenoFile: filename of the phenotype file. This must be present in the working directory at ./in/pheno_files/
    --phenoCol: the column names of the phenotype to be analysed in the file specified in --phenoFile.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 1 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -c,--covarColList: comma separated column names (e.g. age,pc1,pc2) of continuous covariates to include as fixed effects in the file specified in --phenoFile.
    --categCovarColList: comma separated column names of categorical variables to include as fixed effects in the file specified in --phenoFile.
    --sampleIDCol (default: IID): column containing the sample IDs in teh phenotype file, which must match the sample IDs in the plink files.
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
if [[ ${TRAITTYPE} == "" ]]; then
  echo "traitType not set"
  exit 1
fi

if [[ ${PLINK} == "" ]]; then
  echo "plink files plink.{bim,bed,fam} not set"
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

if [[ ${PHENOFILE} == "" ]]; then
  echo "phenoFile not set"
  exit 1
fi

if [[ ${PHENOCOL} == "" ]]; then
  echo "phenoCol not set"
  exit 1
fi

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix to ${PHENOCOL}. Check that this will not overwrite existing files."
  OUT="${PHENOCOL}"
fi

if [[ $COVARCOLLIST = "" ]]; then
  echo "Warning: no continuous fixed effeect covariates included."
fi

if [[ $CATEGCOVARCOLLIST = "" ]]; then
  echo "Warning: no categorical fixed effeect covariates included."
fi

echo "OUT               = ${OUT}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "TRAITTYPE         = ${TRAITTYPE}"
echo "PLINK             = ${PLINK}.{bim/bed/fam}"
echo "SPARSEGRM         = ${SPARSEGRM}"
echo "SPARSEGRMID       = ${SPARSEGRMID}"
echo "PHENOFILE         = ${PHENOFILE}"
echo "PHENOCOL          = ${PHENOCOL}"
echo "COVARCOLLIST      = ${COVARCOLLIST}"
echo "CATEGCOVARCOLLIST = ${CATEGCOVARCOLLIST}"
echo "SAMPLEIDCOL       = ${SAMPLEIDCOL}"

check_container_env $SINGULARITY

echo $SINGULARITY
if [[ ${SINGULARITY} = true ]]; then
  # check if saige sif exists:
  if !(test -f "saige-${saige_version}.sif"); then
      singularity pull "saige-${saige_version}.sif" "docker://wzhou88/saige:${saige_version}"
  fi
else
  docker pull wzhou88/saige:${saige_version}
fi

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

# Get number of threads
n_threads=$(( $(nproc --all) - 1 ))

# Get inverse-normalize flag if trait_type=="quantitative"
if [[ ${TRAITTYPE} == "quantitative" ]]; then
  trait_flags="--traitType=${TRAITTYPE} --invNormalize=TRUE"
else
  trait_flags="--traitType=${TRAITTYPE}"
fi

if [[ ${SINGULARITY} = true ]]; then
  singularity exec \
     --env HOME=${WD} \
     --bind ${WD}/:$HOME/ \
     "saige-${saige_version}.sif" file ${HOME}/in/pheno_files/${PHENOFILE} 
  singularity exec \
    --env HOME=${WD} \
    --bind ${WD}/:$HOME/ \
    "saige-${saige_version}.sif" step1_fitNULLGLMM.R \
      --bedFile ${HOME}/in/plink_vr/${PLINK}.bed \
      --bimFile ${HOME}/in/plink_vr/${PLINK}.bim \
      --famFile ${HOME}/in/plink_vr/${PLINK}.fam \
      --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
      --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
      --useSparseGRMtoFitNULL=TRUE \
      --phenoFile ${HOME}/in/pheno_files/${PHENOFILE} \
      --skipVarianceRatioEstimation FALSE \
      --phenoCol ""${PHENOCOL}"" \
      --covarColList ""${COVARCOLLIST}"" \
      --qCovarColList=""${CATEGCOVARCOLLIST}"" \
      --sampleIDColinphenoFile=${IDCOL} \
      ${trait_flags} \
      --outputPrefix="${HOME}/${OUT}" \
      --IsOverwriteVarianceRatioFile=TRUE \
      --nThreads=${n_threads} \
      --isCateVarianceRatio=TRUE
else
  docker run \
    -e HOME=${WD} \
    -v ${WD}/:$HOME/ \
    "wzhou88/saige:${saige_version}" step1_fitNULLGLMM.R \
      --bedFile ${HOME}/in/plink_vr/${PLINK}.bed \
      --bimFile ${HOME}/in/plink_vr/${PLINK}.bim \
      --famFile ${HOME}/in/plink_vr/${PLINK}.fam \
      --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
      --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
      --useSparseGRMtoFitNULL=TRUE \
      --phenoFile ${HOME}/in/pheno_files/${PHENOFILE} \
      --skipVarianceRatioEstimation FALSE \
      --phenoCol "${PHENOCOL}" \
      --covarColList "${COVARCOLLIST}" \
      --qCovarColList="${CATEGCOVARCOLLIST}" \
      --sampleIDColinphenoFile=${IDCOL} \
      ${trait_flags} \
      --outputPrefix="${HOME}/${OUT}" \
      --IsOverwriteVarianceRatioFile=TRUE \
      --nThreads=${n_threads} \
      --isCateVarianceRatio=TRUE
fi
