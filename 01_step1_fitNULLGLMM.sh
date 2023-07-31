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
  ./resources/plink \
    --bfile "${GENOTYPE_PLINK}" \
    --keep-fam ${SAMPLEIDS} \
    --indep-pairwise 50 5 0.05 \
    --out "${OUT}"

  # Extract set of pruned variants and export to bfile
  ./resources/plink \
    --bfile "${GENOTYPE_PLINK}" \
    --keep-fam ${SAMPLEIDS} \
    --extract "${OUT}.prune.in" \
    --make-bed \
    --out "${OUT}"
  
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
      SINGULARITY=true
      shift # past argument
      shift # past value
      ;;
    -t|--traitType)
      TRAITTYPE="$2"
      if ! ( [[ ${TRAITTYPE} == "quantitative" ]] || [[ ${TRAITTYPE} == "binary" ]] ); then
        echo "Trait type is not in {quantitative,binary}"
        exit 1
      fi
      shift # past argument
      shift # past value
      ;;
    -p|--genotypePlink)
      GENOTYPE_PLINK="$2"
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
    --sampleIDs)
      SAMPLEIDS="$2" 
      shift
      shift
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
    --genotypePlink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
    --phenoFile: filename of the phenotype file. This must be present in the working directory at ./in/pheno_files/
    --phenoCol: the column names of the phenotype to be analysed in the file specified in --phenoFile.
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 1 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -c,--covarColList: comma separated column names (e.g. age,pc1,pc2) of continuous covariates to include as fixed effects in the file specified in --phenoFile.
    --categCovarColList: comma separated column names of categorical variables to include as fixed effects in the file specified in --phenoFile.
    --sampleIDCol (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.
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

if [[ ${SPARSEGRM} == "" || ${SPARSEGRMID} == "" ]]; then
  echo "Sparse GRM .mtx file not set. Generating sparse GRM from genotype or exome sequence data."
  echo "Will attempt to generate GRM"
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
  echo "Warning: no continuous fixed effect covariates included."
fi

if [[ $CATEGCOVARCOLLIST = "" ]]; then
  echo "Warning: no categorical fixed effect covariates included."
fi

echo "OUT               = ${OUT}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "TRAITTYPE         = ${TRAITTYPE}"
echo "PLINK             = ${PLINK_WES}.{bim/bed/fam}"
echo "SPARSEGRM         = ${SPARSEGRM}"
echo "SPARSEGRMID       = ${SPARSEGRMID}"
echo "PHENOFILE         = ${PHENOFILE}"
echo "PHENOCOL          = ${PHENOCOL}"
echo "COVARCOLLIST      = ${COVARCOLLIST}"
echo "CATEGCOVARCOLLIST = ${CATEGCOVARCOLLIST}"
echo "SAMPLEIDS         = ${SAMPLEIDS}"
echo "SAMPLEIDCOL       = ${SAMPLEIDCOL}"

check_container_env $SINGULARITY

if [[ ${SPARSEGRM} == "" || ${SPARSEGRMID} == "" ]]; then
  if [[ ${GENOTYPE_PLINK} == "" ]]; then
    echo "Genotype plink files plink.{bim,bed,fam} not set - cannot generate GRM!"
    exit 1
  fi
  generate_GRM
fi

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

# Get inverse-normalize flag if trait_type=="quantitative"
if [[ ${TRAITTYPE} == "quantitative" ]]; then
  echo "Quantitative trait passed to SAIGE, perform IRNT"
  INVNORMALISE=TRUE
  TOL="0.00001"
else
  echo "Binary trait passed to SAIGE"
  INVNORMALISE=FALSE
  TOL="0.2" # SAIGE DEFAULT
fi

cmd="""step1_fitNULLGLMM.R \
      --plinkFile "${HOME}/${OUT}" \
	    --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
      --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
      --useSparseGRMtoFitNULL=TRUE \
      --phenoFile ${HOME}/in/pheno_files/${PHENOFILE} \
      --skipVarianceRatioEstimation FALSE \
      --traitType=${TRAITTYPE} \
      --invNormalize=${INVNORMALISE} \
      --phenoCol ""${PHENOCOL}"" \
      --covarColList ""${COVARCOLLIST}"" \
      --qCovarColList=""${CATEGCOVARCOLLIST}"" \
      --sampleIDColinphenoFile=${SAMPLEIDCOL} \
      --outputPrefix="${HOME}/${OUT}" \
      --IsOverwriteVarianceRatioFile=TRUE \
      --nThreads=${n_threads} \
      --isCateVarianceRatio=TRUE \
      --tol ${TOL}"""

run_container
