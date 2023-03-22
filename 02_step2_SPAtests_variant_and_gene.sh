# Import functions from other scripts:
. check.sh
. setup.sh

#!/bin/bash
# 01_saige_fit_null
source ./setup.sh
. check.sh

POSITIONAL_ARGS=()

SINGULARITY=false
OUT="out"
TESTTYPE=""
PLINK=""
MODELFILE=""
VARIANCERATIO=""
SPARSEGRM=""
SPARSEGRMID=""
GROUPFILE=""

variant_or_group
  echo "Running variant based tests for all variants in with MAC > 20"
        "saige-${saige_version}.sif" step2_SPAtests.R \
        --GMMATmodelFile ${HOME}/in/model_file/${MODELFILE} \
        --varianceRatioFile ${HOME}/in/variance_ratios/${VARIANCERATIO} \
        --groupFile ${HOME}/in/${GROUPFILE} \


while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    --testType)
      TESTTYPE="$2"
      if ! ( [[ ${TRAITTYPE} = "variant" ]] || [[ ${TRAITTYPE} = "group" ]] ); then
        echo "Test type is not in {variant,group}"
        exit 1
      fi
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      SINGULARITY=true
      shift # past argument
      shift # past value
      ;;
    -p|--plink)
      PLINK="$2"
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
    -h|--help)
      echo "usage: 02_step2_SPAtests_variant_and_gene.sh
  required:
    --testType: type of test {variant,group}.
    -p,--plink: plink filename prefix of bim/bed/fam files. These must be present in the working directory at ./in/plink_for_vr_bed/
    --modelFile: filename of the model file output from step 1. This must be in relation to the working directory.
    --varianceRatio: filename of the varianceRatio file output from step 1. This must be in relation to the working directory.
    --sparseGRM: filename of the sparseGRM .mtx file. This must be present in the working directory at ./in/sparse_grm/
    --sparseGRMID: filename of the sparseGRM ID file. This must be present in the working directory at ./in/sparse_grm/
  optional:
    -o,--outputPrefix:  output prefix of the SAIGE step 2 output.
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -g,--groupFile: required if group test is selected. Filename of the annotation file used for group tests. This must be in relation to the working directory.
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
if [[ ${TESTTYPE} == "" ]]; then
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

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix to 'out'. Check that this will not overwrite existing files."
fi

echo "OUT               = ${OUT}"
echo "TESTTYPE          = ${TESTTYPE}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "PLINK             = ${PLINK}.{bim/bed/fam}"
echo "MODELFILE         = ${MODELFILE}"
echo "VARIANCERATIO     = ${VARIANCERATIO}"
echo "GROUPFILE         = ${GROUPFILE}"
echo "SPARSEGRM         = ${SPARSEGRM}"
echo "SPARSEGRMID       = ${SPARSEGRMID}"

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

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

if [[ "$test_type" = "variant" ]]; then
  echo "Running variant based tests for all variants in with MAC > 20"
  if [[ ${SINGULARITY} = true ]]; then
    singularity exec \
      --env HOME=${WD} \
      --bind ${WD}/:$HOME/ \
      "saige-${saige_version}.sif" step2_SPAtests.R \
        --bedFile ${HOME}/in/plink_exome/${PLINK}.bed \
        --bimFile ${HOME}/in/plink_exome/${PLINK}.bim \
        --famFile ${HOME}/in/plink_exome/${PLINK}.fam \
        --minMAF=0 \
        --minMAC=20 \
        --GMMATmodelFile ${HOME}/${MODELFILE} \
        --varianceRatioFile ${HOME}/${VARIANCERATIO} \
        --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.1 \
        --is_output_moreDetails=TRUE \
        --is_fastTest=TRUE \
        --SAIGEOutputFile=${HOME}/${OUT}_variant.tsv
  else
    # Check --AlleleOrder=ref-first
    docker run \
      -e HOME=${WD} \
      -v ${WD}/:$HOME/ \
      "wzhou88/saige:${saige_version}" step2_SPAtests.R \
        --bedFile ${HOME}/in/plink_exome/${PLINK}.bed \
        --bimFile ${HOME}/in/plink_exome/${PLINK}.bim \
        --famFile ${HOME}/in/plink_exome/${PLINK}.fam \
        --minMAF=0 \
        --minMAC=20 \
        --GMMATmodelFile ${HOME}/${MODELFILE} \
        --varianceRatioFile ${HOME}/${VARIANCERATIO} \
        --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.1 \
        --is_output_moreDetails=TRUE \
        --is_fastTest=TRUE \
        --SAIGEOutputFile=${HOME}/${OUT}_variant.tsv
  fi
else
    echo '''Running gene based tests and variant based tests for all variants present in the annotations.
    This includes the collapsed variants in the set-based tests'''
    if [[ ${SINGULARITY} = true ]]; then
      singularity exec \
        --env HOME=${WD} \
        --bind ${WD}/:$HOME/ \
        "saige-${saige_version}.sif" step2_SPAtests.R \
        --bedFile ${HOME}/in/plink_exome/${PLINK}.bed \
        --bimFile ${HOME}/in/plink_exome/${PLINK}.bim \
        --famFile ${HOME}/in/plink_exome/${PLINK}.fam \
        --minMAF=0 \
        --minMAC=0.5 \
        --GMMATmodelFile ${HOME}/${MODELFILE} \
        --varianceRatioFile ${HOME}/${VARIANCERATIO} \
        --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
        --LOCO=FALSE \
        --groupFile ${HOME}/in/${GROUPFILE} \
        --annotation_in_groupTest="pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous" \
        --maxMAF_in_groupTest=0.0001,0.001,0.01 \
        --is_output_markerList_in_groupTest=TRUE \
        --is_single_in_groupTest=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.1 \
        --is_fastTest=TRUE \
        --SAIGEOutputFile=${HOME}/${OUT}_variant_and_group.tsv
fi
    # Check --AlleleOrder=ref-first
    docker run \
      -e HOME=${WD} \
      -v ${WD}/:$HOME/ \
      "wzhou88/saige:${saige_version}" step2_SPAtests.R \
        --bedFile ${HOME}/in/plink_exome/${PLINK}.bed \
        --bimFile ${HOME}/in/plink_exome/${PLINK}.bim \
        --famFile ${HOME}/in/plink_exome/${PLINK}.fam \
        --minMAF=0 \
        --minMAC=0.5 \
        --GMMATmodelFile ${HOME}/${MODELFILE} \
        --varianceRatioFile ${HOME}/${VARIANCERATIO} \
        --sparseGRMFile ${HOME}/in/sparse_grm/${SPARSEGRM} \
        --sparseGRMSampleIDFile ${HOME}/in/sparse_grm/${SPARSEGRMID} \
        --LOCO=FALSE \
        --groupFile ${HOME}/in/${GROUPFILE} \
        --annotation_in_groupTest="pLoF,damaging_missense,other_missense,synonymous,pLoF:damaging_missense,pLoF:damaging_missense:other_missense:synonymous" \
        --maxMAF_in_groupTest=0.0001,0.001,0.01 \
        --is_output_markerList_in_groupTest=TRUE \
        --is_single_in_groupTest=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.1 \
        --is_fastTest=TRUE \
        --SAIGEOutputFile=${HOME}/${OUT}_variant_and_group.tsv
fi
