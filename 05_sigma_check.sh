#!/bin/bash
# Is the Sigma step 2 will solve against positive definite? One null model and
# the GRM it was fitted on, at the cutoff steps 1 and 2 use. Run it on EVERY
# null model before step 2 (saige-slim LEDGER #183, 2026-10-28).
#
# Why: on the sparse-GRM band step 2 takes the variance of a test from a
# conjugate-gradient solve against Sigma = tau1 I + tau2 K, where K is the GRM
# with every entry below --relatednessCutoff deleted. Deleting entries does not
# keep a matrix positive definite, and the solve checks only its residual, so
# on an indefinite Sigma a burden test can get a negative variance (printed as
# p = 1, BETA NA) or a tiny one (an inflated Z printed as a valid p: 1e-73 on
# a random-score control in All of Us). The step-2 guard repairs such tests by
# default from image 66b80710 and names them in <prefix>.txt.sigmaGuard.txt,
# but a GRM whose Sigma is indefinite is still the wrong GRM: in All of Us the
# admixed amr matrix kept ~90 % ancestry background at 0.125 and was rebuilt
# from KING-robust kinship (see the README). aou_sigma_pd_check.R is on the
# image's PATH: it cuts the thinned GRM into families, eigendecomposes each
# block (up to --maxDenseComponent samples), hands the whole sparse Sigma to a
# Cholesky once, and writes a small key/value report. Nothing N x N is formed.

source ./run_container.sh

POSITIONAL_ARGS=()

MODEL_FILE=""
SPARSE_GRM_FILE=""
SPARSE_GRM_ID_FILE=""
OUT_FILE=""
LABEL=""
RELCUTOFF="0.05"
MAX_DENSE="2000"
PLANT=""
SINGULARITY=false

WD=$(pwd)
HOME=$WD

while [[ $# -gt 0 ]]; do
  case $1 in
    --modelFile)
      MODEL_FILE="$2"
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
    --relatednessCutoff)
      RELCUTOFF="$2"
      shift # past argument
      shift # past value
      ;;
    --label)
      LABEL="$2"
      shift # past argument
      shift # past value
      ;;
    --maxDenseComponent)
      MAX_DENSE="$2"
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
    --plantIndefinite)
      PLANT="--plantIndefinite"
      shift # past argument
      ;;
    --dryRun)
      DRYRUN=true
      shift # past argument
      ;;
    -h|--help)
      echo "usage: 05_sigma_check.sh
  required:
    --modelFile: the step-1 model (<out_step1>.rda).
    --sparseGRM: filename of the sparseGRM .mtx file the model was fitted on.
    --sparseGRMID: filename of the sparseGRM ID file.
    Each is relative to the current working directory, or an absolute path under it or under a directory listed,
    colon-separated, in SAIGE_EXTRA_MOUNTS (bound read-only).
  optional:
    --relatednessCutoff (default: 0.05): the cutoff steps 1 and 2 use (the GRM is thinned to entries above it). Pass
      the value steps 0-2 used; the check describes the Sigma step 2 will build, so it must be the same.
    -o,--outputFile (default: <modelFile without .rda>.sigma.tsv): the report, one 'label key value' row per line
      (tau1, tau2, pairs kept, components evaluated / indefinite, the smallest eigenvalue of Sigma and of the diagonal it
      sits on, and the verdict). The log is <outputFile>.log.
    --label (default: the model file's name): the report's label column.
    --maxDenseComponent (default: 2000): the largest family block eigendecomposed; larger blocks are judged by the
      Cholesky of the whole Sigma alone.
    -s,--isSingularity (default: false): is singularity (or apptainer) available? If not, it is assumed that docker is available.
    --dryRun: print the command instead of running it.
    --plantIndefinite: the positive control -- a planted indefinite block, which must be reported as NOT positive definite.
  Exit 0 when Sigma is positive definite (a warning when it has near-null directions), 2 when it is NOT, with the
  rows of the report on stdout. A NOT verdict means the GRM, thinned at this cutoff, is not a covariance matrix: rebuild
  it (raise the cutoff only if the dropped pairs are background, not relatives; in an admixed cohort build the GRM
  from KING-robust kinship instead -- see the README), refit step 1 on it, and run this check again. Do not run
  step 2 on a model whose Sigma failed this check.
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

if [[ ${MODEL_FILE} == "" || ${SPARSE_GRM_FILE} == "" || ${SPARSE_GRM_ID_FILE} == "" ]]; then
  echo "--modelFile, --sparseGRM and --sparseGRMID are required"
  exit 1
fi
for f in "${MODEL_FILE}" "${SPARSE_GRM_FILE}" "${SPARSE_GRM_ID_FILE}"; do
  [[ -s ${f} ]] || { echo "not found or empty: ${f}" >&2; exit 1; }
done
[[ ${OUT_FILE} == "" ]] && OUT_FILE="${MODEL_FILE%.rda}.sigma.tsv"
[[ ${LABEL} == "" ]] && LABEL=$(basename "${MODEL_FILE%.rda}")

check_relcutoff "${RELCUTOFF}"
grm_density_check "${SPARSE_GRM_FILE}" "${RELCUTOFF}" || true   # its density, and a cutoff that differs from the GRM's
check_container_env $SINGULARITY

# relative paths are inside the working directory; an absolute path is used as
# given, and must be under the working directory or a SAIGE_EXTRA_MOUNTS
# directory (see run_container.sh)
container_path () { [[ $1 == /* ]] && echo "$1" || echo "${HOME}/$1"; }

log="${OUT_FILE}.log"
cmd="aou_sigma_pd_check.R \
  --modelFile=$(container_path "${MODEL_FILE}") \
  --sparseGRMFile=$(container_path "${SPARSE_GRM_FILE}") \
  --sparseGRMSampleIDFile=$(container_path "${SPARSE_GRM_ID_FILE}") \
  --relatednessCutoff=${RELCUTOFF} \
  --maxDenseComponent=${MAX_DENSE} \
  --label=${LABEL} \
  --reportFile=$(container_path "${OUT_FILE}") ${PLANT}"
echo "Sigma check for ${MODEL_FILE} on ${SPARSE_GRM_FILE} at --relatednessCutoff ${RELCUTOFF}; log: ${log}"
if [[ ${DRYRUN:-false} = true ]]; then
  run_container; exit $?
fi
run_container > "${log}" 2>&1; rc=$?
if (( rc != 0 )) || [[ ! -s ${OUT_FILE} ]]; then
  echo "  FAILED (exit ${rc}); the end of ${log}:" >&2
  tail -n 10 "${log}" >&2
  exit 1
fi
echo "Contents of ${OUT_FILE}:"
cat "${OUT_FILE}"
verdict=$(awk -F'\t' '$2 == "verdict" {v = $3} END {print v}' "${OUT_FILE}")
case "${verdict}" in
  "positive definite")
    echo "Sigma is positive definite: step 2 can run on this model."
    ;;
  "positive definite with near-null directions")
    echo "WARNING: Sigma is positive definite but has near-null directions (a family block with an eigenvalue" >&2
    echo "  near zero). Step 2 can run; its guard (<prefix>.txt.sigmaGuard.txt) will name any test that lands on" >&2
    echo "  such a direction." >&2
    ;;
  "tau2 = 0: Sigma is diagonal, positive definite by construction")
    echo "Sigma is diagonal (no random effect was fitted): positive definite by construction."
    ;;
  *)
    {
      echo "################################################################################"
      echo "REFUSED: Sigma is NOT positive definite for ${MODEL_FILE} at --relatednessCutoff ${RELCUTOFF}"
      echo "  (verdict: '${verdict}'). The GRM thinned at this cutoff is not a covariance matrix, and"
      echo "  step 2's sparse-GRM variances on it are the thinning's artefact: negative (p 1, BETA NA)"
      echo "  or tiny (an inflated Z printed as a valid p). Rebuild the GRM -- in an admixed cohort from"
      echo "  KING-robust kinship (README, step 0) -- refit step 1 on it, and run this check again."
      echo "  Do NOT run step 2 on this model."
      echo "################################################################################"
    } >&2
    exit 2
    ;;
esac
