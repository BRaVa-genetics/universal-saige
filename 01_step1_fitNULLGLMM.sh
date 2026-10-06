#!/bin/bash

source ./run_container.sh
source ./check_pheno.sh

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
HOME=$WD

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      echo "out: $OUT"
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      SINGULARITY="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--traitType)
      TRAITTYPE="$2"
      if ! ( [[ ${TRAITTYPE} == "quantitative" ]] || [[ ${TRAITTYPE} == "binary" ]] ); then
        echo "Trait type is not in {quantitative,binary}"
        exit 1
      fi
      echo "trait type: ${TRAITTYPE}"
      shift # past argument
      shift # past value
      ;;
    -p|--genotypePlink)
      GENOTYPE_PLINK="$2"
      echo "genotype plink: ${GENOTYPE_PLINK}"
      shift # past argument
      shift # past value
      ;;
    --genotypePgen)
      GENOTYPE_PGEN="$2"
      echo "genotype pgen: ${GENOTYPE_PGEN}"
      shift # past argument
      shift # past value
      ;;
    --dryRun)
      DRYRUN=true
      shift # past argument
      ;;
    --sparseGRM)
      SPARSEGRM="$2"
      echo "sparse GRM: ${SPARSEGRM}"
      shift # past argument
      shift # past value
      ;;
    --sparseGRMID)
      SPARSEGRMID="$2"
      echo "sparse GRM ID: ${SPARSEGRMID}"
      shift # past argument
      shift # past value
      ;;
    --phenoFile)
      PHENOFILE="$2"
      echo "phenotype file: ${PHENOFILE}"
      shift # past argument
      shift # past value
      ;;
    --phenoCol)
      PHENOCOL="$2"
      echo "phenotype column: ${PHENOCOL}"
      shift # past argument
      shift # past value
      ;;
    -c|--covarColList)
      COVARCOLLIST="$2"
      echo "covariate column list: ${COVARCOLLIST}"
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
      echo "sample ID column: ${SAMPLEIDCOL}"
      shift # past argument
      shift # past value
      ;;
    --sex)
      SEX="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      echo "usage: 01_step1_fitNULLGLMM.sh
  required:
    -t,--traitType: type of the trait {quantitative,binary}.
    --genotypePlink: plink 1 filename prefix of bim/bed/fam files (the variance-ratio markers from step 0). This must be relative to, and contained within, the working directory from which the docker/singularity was launched.
    --genotypePgen: plink 2 filename prefix of pgen/pvar/psam files, the alternative to --genotypePlink.
    --sparseGRM: filename of the sparseGRM .mtx file. This must be relative to, and contained within, the current working directory.
    --sparseGRMID: filename of the sparseGRM ID file. This must be relative to, and contained within, the current working directory.
    --phenoFile: filename of the phenotype file. This must be relative to, and contained within, the current working directory.
    --phenoCol: the column names of the phenotype to be analysed in the file specified in --phenoFile.
  optional:
    -o,--outputPrefix: output prefix of the SAIGE step 1 output. This must be relative to, and contained within, the current working directory.
    -s,--isSingularity (default: false): is singularity (or apptainer) available? If not, it is assumed that docker is available.
    -c,--covarColList: comma separated column names (e.g. age,pc1,pc2) of continuous covariates to include as fixed effects in the file specified in --phenoFile.
    --categCovarColList: comma separated column names of categorical variables to include as fixed effects in the file specified in --phenoFile.
    --sampleIDCol (default: IID): column containing the sample IDs in the phenotype file, which must match the sample IDs in the plink files.
    --sex ('M' or 'F'): for a sex-specific trait. Every sample with a non-missing phenotype must share one value of the
      'sex' column (and equal --sex, when that column is coded M/F). Do not also pass sex as a covariate.
    --dryRun: print the SAIGE command instead of running it.
  fit gates (the slim build): a binary trait fitted on fewer than 100 cases, a categorical covariate level with fewer
    than 10 cases or controls, a separated covariate model or a fit that did not converge is REFUSED, with the gate
    named in the log. SAIGE_FIT_GATES=0 in the environment downgrades the refusal to a warning (not recommended).
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

if [[ ${SAMPLEIDS} != "" ]]; then
  SAMPLEIDS=${HOME}/$SAMPLEIDS
fi

if [[ ${SPARSEGRM} == "" || ${SPARSEGRMID} == "" ]]; then
  echo "Sparse GRM .mtx file not set. Generate a GRM in step 0."
fi

if [[ ${PHENOFILE} == "" ]]; then
  echo "phenoFile not set"
  exit 1
fi

if [[ ${PHENOCOL} == "" ]]; then
  echo "phenoCol not set"
  exit 1
fi

# A sex-specific trait: every sample with a phenotype must carry the same value
# in the 'sex' column. The coding of that column varies by cohort, so only an
# M/F coding can also be compared with --sex itself.
if [[ ${SEX:-} != "" ]]; then
  if [[ ${SEX} != "M" && ${SEX} != "F" ]]; then
    echo "--sex must be M or F"
    exit 1
  fi
  fs=' '   # awk's whitespace splitting, unless the file is tab-delimited (empty fields keep their place)
  [[ $(head -n 1 "${PHENOFILE}") == *$'\t'* ]] && fs='\t'
  header=$(head -n 1 "${PHENOFILE}" | awk -F"${fs}" '{for (i = 1; i <= NF; i++) print $i}')
  sex_col_num=$(grep -n -x 'sex' <<< "${header}" | cut -d: -f1)
  pheno_col_num=$(grep -n -x -F "${PHENOCOL}" <<< "${header}" | cut -d: -f1)
  if [[ ${sex_col_num} == "" || ${pheno_col_num} == "" ]]; then
    echo "--sex ${SEX}: ${PHENOFILE} needs a 'sex' column and a '${PHENOCOL}' column"
    exit 1
  fi
  sex_values=$(awk -F"${fs}" -v s="${sex_col_num}" -v p="${pheno_col_num}" \
    'NR > 1 && $p != "NA" && $p != "" {print $s}' "${PHENOFILE}" | sort | uniq -c)
  if (( $(grep -c . <<< "${sex_values}") > 1 )); then
    echo "--sex ${SEX}: samples with a non-missing ${PHENOCOL} have more than one value of 'sex' (count, value):"
    echo "${sex_values}"
    exit 1
  fi
  observed_sex=$(awk '{print $2}' <<< "${sex_values}")
  if [[ ( ${observed_sex} == "M" || ${observed_sex} == "F" ) && ${observed_sex} != "${SEX}" ]]; then
    echo "--sex ${SEX}: samples with a non-missing ${PHENOCOL} are all sex ${observed_sex}"
    exit 1
  fi
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


if is_valid_r_var "$PHENOCOL"; then
    echo "The variable name '$PHENOCOL' is valid for an R variable."
fi

if [[ "$PHENOCOL" =~ .*"-".* || "$PHENOCOL" =~ .*",".* || "$PHENOCOL" =~ .*"=".* ]]; then
  echo "Phenotype name cannot contain \"-\" or \",\" or \"=\""
  exit 1
fi

GENOTYPE_PLINK="${GENOTYPE_PLINK:-}"; GENOTYPE_PGEN="${GENOTYPE_PGEN:-}"
if [[ ${GENOTYPE_PLINK} == "" && ${GENOTYPE_PGEN} == "" ]]; then
  echo "Genotype file not specified (--genotypePlink or --genotypePgen). Note, this is the file used to determine the variance ratio."
  echo "You can automatically generate it in step 0, the collection of plink files end with .plink_for_var_ratio.{bim,bed,fam}."
  exit 1
fi
if [[ ${GENOTYPE_PLINK} != "" && ${GENOTYPE_PGEN} != "" ]]; then
  echo "Pass ONE of --genotypePlink and --genotypePgen"
  exit 1
fi

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

# Get number of threads
n_threads=$(( $(ncpu) - 1 )); (( n_threads < 1 )) && n_threads=1

# Get inverse-normalize flag if trait_type=="quantitative"
# --useSparseGRMforVarRatio=TRUE: step 2 is given the sparse GRM, and a model
# whose variance-ratio file carries no `sparse` rows is REFUSED there (the
# All of Us models carry them). Without it every step-2 run fails on load.
# --tol is 0.02 (SAIGE's default) for BOTH trait types: the All of Us
# production runs used it for both, and the earlier 1e-5 for quantitative
# traits bought nothing measurable at a real cost in fit time.
TOL="0.02"
if [[ ${TRAITTYPE} == "quantitative" ]]; then
  echo "Quantitative trait passed to SAIGE, perform IRNT"
  INVNORMALISE=TRUE
else
  echo "Binary trait passed to SAIGE"
  INVNORMALISE=FALSE
fi

if [[ ${GENOTYPE_PGEN} != "" ]]; then
  GENO_ARGS="--pgenFile ${HOME}/${GENOTYPE_PGEN}.pgen --pvarFile ${HOME}/${GENOTYPE_PGEN}.pvar --psamFile ${HOME}/${GENOTYPE_PGEN}.psam"
else
  GENO_ARGS="--plinkFile ${HOME}/${GENOTYPE_PLINK}"
fi

cmd="""step1_fitNULLGLMM.R \
      ${GENO_ARGS} \
      --relatednessCutoff 0.05 \
      --sparseGRMFile ${HOME}/${SPARSEGRM} \
      --sparseGRMSampleIDFile ${HOME}/${SPARSEGRMID} \
      --useSparseGRMtoFitNULL=TRUE \
      --useSparseGRMforVarRatio=TRUE \
      --phenoFile ${HOME}/${PHENOFILE} \
      --skipVarianceRatioEstimation FALSE \
      --traitType=${TRAITTYPE} \
      --invNormalize=${INVNORMALISE} \
      --phenoCol ""${PHENOCOL}"" \
      --covarColList=""${COVARCOLLIST}"" \
      --qCovarColList=""${CATEGCOVARCOLLIST}"" \
      --sampleIDColinphenoFile=${SAMPLEIDCOL} \
      --outputPrefix="${HOME}/${OUT}" \
      --IsOverwriteVarianceRatioFile=TRUE \
      --nThreads=${n_threads} \
      --isCateVarianceRatio=TRUE \
      --tol ${TOL} \
      --SampleIDIncludeFile=${SAMPLEIDS} \
      --isCovariateOffset TRUE"""

# A REFUSED fit (a gate, or any error after the output files were opened)
# leaves 0-byte <prefix>.rda and .varianceRatio.txt behind, which a later step
# 2 would load blindly and fail on uselessly. Remove them, so a refused trait
# leaves nothing that looks like a model; the log has the reason.
set +e
run_container; rc=$?
set -e
if (( rc != 0 )); then
  for f in "${OUT}.rda" "${OUT}.varianceRatio.txt"; do
    [[ -f ${f} && ! -s ${f} ]] && rm -f "${f}" && echo "removed empty ${f} (the fit was refused or failed; see the log)" >&2
  done
  exit ${rc}
fi
