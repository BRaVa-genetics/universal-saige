#!/bin/bash
# Effective sample size (Nglmm) of each phenotype's null GLMM on the sparse
# GRM, one row per phenotype in a csv. extractNglmm.R is on the slim image's
# PATH; it prints a line "Nglmm <value>" and writes no files.

source ./run_container.sh

POSITIONAL_ARGS=()

binary_phenos=""
cont_phenos=""
PHENO_FILE=""
COVAR_LIST=""
SPARSE_GRM_FILE=""
SPARSE_GRM_ID_FILE=""
OUT_FILE="neff.csv"
SINGULARITY=false

WD=$(pwd)
HOME=$WD

while [[ $# -gt 0 ]]; do
  case $1 in
    --binaryPhenos)
      binary_phenos="$2"
      shift # past argument
      shift # past value
      ;;
    --contPhenos)
      cont_phenos="$2"
      shift # past argument
      shift # past value
      ;;
    --phenoFile)
      PHENO_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    --covarList)
      COVAR_LIST="$2"
      shift # past argument
      shift # past value
      ;;
    --sparseGRM)
      SPARSE_GRM_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    --sparseGRMID)
      SPARSE_GRM_ID_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--outputFile)
      OUT_FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      SINGULARITY="$2"
      shift # past argument
      shift # past value
      ;;
    --dryRun)
      DRYRUN=true
      shift # past argument
      ;;
    -h|--help)
      echo "usage: 03_estimate_nGlmm.sh
  required:
    --binaryPhenos: space separated list of binary phenotypes.
    --contPhenos: space separated list of continuous phenotypes.
    --phenoFile: filename of the phenotype file.
    --sparseGRM: filename of the sparseGRM .mtx file.
    --sparseGRMID: filename of the sparseGRM ID file.
    Each is relative to the current working directory, or an absolute path under it, under /mnt/project (the UKB RAP
    project mount, bound read-only when it exists) or under a directory in SAIGE_EXTRA_MOUNTS (colon-separated).
  optional:
    --covarList: comma separated column names of covariates in --phenoFile.
    -o,--outputFile (default: neff.csv): the csv written, one 'pheno,nglmm' row per phenotype. Each phenotype's log is
      <outputFile without .csv>.<pheno>.log.
    -s,--isSingularity (default: false): is singularity (or apptainer) available? If not, it is assumed that docker is available.
    --dryRun: print the commands instead of running them.
  A phenotype whose fit fails (or is refused by the fit gates) gets no row, and the script exits non-zero.
      "
      exit 0
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

if [[ ${PHENO_FILE} == "" || ${SPARSE_GRM_FILE} == "" || ${SPARSE_GRM_ID_FILE} == "" ]]; then
  echo "--phenoFile, --sparseGRM and --sparseGRMID are required"
  exit 1
fi
if [[ ${binary_phenos} == "" && ${cont_phenos} == "" ]]; then
  echo "no phenotypes: pass --binaryPhenos and/or --contPhenos"
  exit 1
fi

check_container_env $SINGULARITY

# relative paths are inside the working directory; an absolute path is used as
# given, and must be under the working directory, /mnt/project or a
# SAIGE_EXTRA_MOUNTS directory (see run_container.sh)
container_path () { [[ $1 == /* ]] && echo "$1" || echo "${HOME}/$1"; }

COVAR_ARGS=""
[[ ${COVAR_LIST} != "" ]] && COVAR_ARGS="--covarColList ${COVAR_LIST}"

echo "pheno,nglmm" > "${OUT_FILE}"
failed=()
estimate () {   # $1 = phenotype, $2 = trait type
  local log="${OUT_FILE%.csv}.$1.log" rc nglmm
  cmd="extractNglmm.R \
    --phenoFile $(container_path "${PHENO_FILE}") \
    --phenoCol $1 \
    ${COVAR_ARGS} \
    --traitType $2 \
    --sparseGRMFile $(container_path "${SPARSE_GRM_FILE}") \
    --sparseGRMSampleIDFile $(container_path "${SPARSE_GRM_ID_FILE}") \
    --useSparseGRMtoFitNULL TRUE"
  echo "Estimating Nglmm for $1 ($2); log: ${log}"
  run_container > "${log}" 2>&1; rc=$?
  if [[ ${DRYRUN:-false} = true ]]; then
    cat "${log}"; (( rc == 0 )) || failed+=("$1"); return
  fi
  nglmm=$(awk '$1 == "Nglmm" {v = $2} END {print v}' "${log}")
  if (( rc != 0 )) || [[ ${nglmm} == "" ]]; then
    echo "  FAILED (exit ${rc}); the end of ${log}:" >&2
    tail -n 5 "${log}" >&2
    failed+=("$1")
  else
    echo "$1,${nglmm}" >> "${OUT_FILE}"
    echo "  Nglmm = ${nglmm}"
  fi
}

for pheno in ${cont_phenos}; do estimate "${pheno}" quantitative; done
for pheno in ${binary_phenos}; do estimate "${pheno}" binary; done

echo "Contents of ${OUT_FILE}:"
cat "${OUT_FILE}"
if (( ${#failed[@]} > 0 )); then
  echo "No Nglmm for: ${failed[*]}" >&2
  exit 1
fi
