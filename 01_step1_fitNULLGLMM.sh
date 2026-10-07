#!/bin/bash

source ./run_container.sh
source ./check_pheno.sh

POSITIONAL_ARGS=()

SINGULARITY=false
SAMPLEIDCOL="IID"
RELCUTOFF="0.05"
FORCE_DENSE_GRM=false
FEMALE_CODE="0"
MALE_CODE="1"
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
      # an empty value (e.g. an unset variable) would let --sampleIDs swallow the next flag
      [[ ${SAMPLEIDS} == "" || ${SAMPLEIDS} == -* ]] && { echo "--sampleIDs needs a file; to use every sample, leave the flag out"; exit 1; }
      shift
      shift
      ;;
    -i|--sampleIDCol)
      SAMPLEIDCOL="$2"
      echo "sample ID column: ${SAMPLEIDCOL}"
      shift # past argument
      shift # past value
      ;;
    --relatednessCutoff)
      RELCUTOFF="$2"
      shift # past argument
      shift # past value
      ;;
    --forceDenseGRM)
      FORCE_DENSE_GRM=true
      shift # past argument
      ;;
    --sex)
      SEX="$2"
      shift # past argument
      shift # past value
      ;;
    --femaleCode)
      FEMALE_CODE="$2"
      shift # past argument
      shift # past value
      ;;
    --maleCode)
      MALE_CODE="$2"
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
    --forceDenseGRM: fit even when the sparse GRM is dense (more than 100 relatives per sample), which is otherwise
      refused; not recommended.
    --relatednessCutoff (default: 0.05): the GRM is thinned to entries at or above it. MUST equal step 0's and step 2's
      (All of Us: 0.05, and 0.125 for amr); nothing in SAIGE checks it.
    --sex ('M' or 'F'): for a sex-specific trait. SAIGE fits only the samples whose 'sex' column holds that sex's code
      (--FemaleOnly/--MaleOnly, as All of Us did) and drops the rest. Leave every sex term out of the covariates.
    --femaleCode, --maleCode (default 0, 1: BRaVa's phenotype coding): the values of the 'sex' column for each sex, when it
      is numeric; a column of M and F is read as such. SAIGE's own default, and the All of Us file, use 1 for female.
    --dryRun: print the SAIGE command instead of running it.
  fit gates (the slim build): a binary trait fitted on fewer than 100 cases, a categorical covariate level with fewer
    than 10 cases or controls, a separated covariate model or a fit that did not converge is REFUSED, with the gate
    named in the log. SAIGE_FIT_GATES=0 in the environment downgrades the refusal to a warning (not recommended).
  dense GRM: a sparse GRM with more than 100 relatives per sample on average is REFUSED before the fit (All of Us amr
    at 0.05: ~644, a fit that never finished). Raise --relatednessCutoff in steps 0, 1 and 2; --forceDenseGRM
    fits anyway (not recommended).
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
check_relcutoff "${RELCUTOFF}"
# A dense GRM is refused before the fit is paid for (see grm_density_check)
dense=0; grm_density_check "${SPARSEGRM}" "${RELCUTOFF}" || dense=$?
if (( dense == 3 )); then
  if [[ ${FORCE_DENSE_GRM} = true ]]; then
    echo "--forceDenseGRM: fitting on the dense GRM anyway" >&2
  else
    echo "REFUSED: the sparse GRM is too dense to fit (above). Rerun step 0 with a higher --relatednessCutoff and" >&2
    echo "pass the same value here and to step 2, or pass --forceDenseGRM to fit anyway (not recommended)." >&2
    exit 1
  fi
fi
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

# A sex-specific trait, as All of Us fitted them: SAIGE keeps the samples whose
# 'sex' equals the code (--FemaleOnly/--MaleOnly) and appends _FemaleOnly or
# _MaleOnly to the output prefix, renamed back after the fit. The code is the
# phenotype file's: F/M when the column holds F and M, else --femaleCode and
# --maleCode, default 0 and 1 as BRaVa's phenotype curation writes them
# (extract_BRaVa_phenotypes.r: female = sex 0). SAIGE's defaults and the All of
# Us file use 1 for female, so the counts printed here are the check that the
# code names the sex meant; no sample of that code with a phenotype is refused.
SEX_ARGS=""; SEX_SUFFIX=""
if [[ ${SEX:-} != "" ]]; then
  if [[ ${SEX} != "M" && ${SEX} != "F" ]]; then
    echo "--sex must be M or F"
    exit 1
  fi
  if [[ ,${COVARCOLLIST},${CATEGCOVARCOLLIST}, == *,sex,* ]]; then
    echo "--sex ${SEX}: 'sex' is also a covariate; leave every sex term (sex, age_sex, age2_sex) out of a sex-specific fit"
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
  if awk -F"${fs}" -v s="${sex_col_num}" 'NR > 1 && $s != "" && $s != "NA" && $s != "M" && $s != "F" { exit 1 }' "${PHENOFILE}"; then
    female="F"; male="M"
  else
    female="${FEMALE_CODE}"; male="${MALE_CODE}"
  fi
  if [[ ${SEX} == "F" ]]; then code="${female}"; other="${male}"; else code="${male}"; other="${female}"; fi
  # one pass: the sex column's values; samples with a phenotype kept and dropped;
  # for a binary trait, the cases (1) kept and dropped
  read -r n_kept n_dropped cases_kept cases_dropped sex_values <<< "$(awk -F"${fs}" -v s="${sex_col_num}" -v p="${pheno_col_num}" -v c="${code}" '
    NR > 1 && $s != "" && $s != "NA" && !($s in seen) { seen[$s]; vals = vals (vals == "" ? "" : ",") $s }
    NR > 1 && $p != "NA" && $p != "" { if ($s == c) { k++; ck += ($p == 1) } else { d++; cd += ($p == 1) } }
    END { print k + 0, d + 0, ck + 0, cd + 0, vals }' "${PHENOFILE}")"
  case_note=""; [[ ${TRAITTYPE} == "binary" ]] && case_note=" (${cases_kept} cases kept, ${cases_dropped} dropped)"
  echo "--sex ${SEX}: the 'sex' column holds ${sex_values}; ${SEX} is coded ${code}"
  # refusals: the code names nobody, or nobody with a phenotype, or (binary) no case
  if [[ ,${sex_values}, != *,${code},* ]]; then
    echo "REFUSED: --sex ${SEX} is coded ${code}, which is not a value of the 'sex' column (${sex_values})." >&2
    echo "Set --femaleCode/--maleCode to the file's coding (default 0 female, 1 male: BRaVa's; SAIGE and All of Us use 1 for female)." >&2
    exit 1
  fi
  if (( n_kept == 0 )); then
    echo "REFUSED: no sample with sex == ${code} has a non-missing ${PHENOCOL} (${n_dropped} of the other sex do)." >&2
    echo "Check --femaleCode/--maleCode (default 0 female, 1 male: BRaVa's coding; SAIGE and All of Us use 1 for female)." >&2
    exit 1
  fi
  if [[ ${TRAITTYPE} == "binary" ]] && (( cases_kept == 0 && cases_dropped > 0 )); then
    echo "REFUSED: all ${cases_dropped} cases of ${PHENOCOL} have sex == ${other}, the sex being dropped: the codes look flipped." >&2
    echo "Check --femaleCode/--maleCode (default 0 female, 1 male: BRaVa's coding; SAIGE and All of Us use 1 for female)." >&2
    exit 1
  fi
  # loud warnings: what flipped or mismatched codes look like when they still leave samples to fit
  warn=()
  # Genetic sex, where the genotype file has it, is a reference whose coding the
  # format fixes (PLINK: 1 male, 2 female), so it can check the phenotype file's
  # codes. The phenotype file still decides who is fitted (SAIGE filters on it).
  # A .psam may have no SEX column, and a .fam may carry 0 (unknown) for all, as
  # data converted from VCF usually does: then the codes cannot be checked.
  if [[ ${GENOTYPE_PGEN} != "" ]]; then genofile="${GENOTYPE_PGEN}.psam"; else genofile="${GENOTYPE_PLINK}.fam"; fi
  id_col_num=$(grep -n -x -F "${SAMPLEIDCOL}" <<< "${header}" | cut -d: -f1)
  if [[ ! -r ${genofile} || ${id_col_num} == "" ]]; then
    echo "--sex ${SEX}: genetic sex not checked (${genofile} not readable, or no ${SAMPLEIDCOL} column in ${PHENOFILE})"
  else
    # genetic sex per IID: 1 male, 2 female, anything else unknown
    read -r g_males g_females n_compared n_discordant <<< "$(awk -v psam="${GENOTYPE_PGEN}" -v gf="${genofile}" \
        -v idc="${id_col_num}" -v s="${sex_col_num}" -v fc="${female}" -v mc="${male}" '
      FILENAME == gf {
        if (psam != "") {                       # .psam: a #FID/#IID header naming IID and SEX
          if (FNR == 1) { for (i = 1; i <= NF; i++) { h = $i; sub(/^#/, "", h); if (h == "IID") ic = i; if (toupper(h) == "SEX") sc = i }; next }
          id = $ic; g = sc ? $sc : ""
        } else { id = $2; g = $5 }              # .fam: FID IID father mother sex phenotype
        g = (g == "1" || g == "M" || g == "m") ? 1 : (g == "2" || g == "F" || g == "f") ? 2 : 0
        if (g) { gsex[id] = g; if (g == 1) nm++; else nf++ }
        next
      }
      FNR > 1 && ($idc in gsex) && ($s == fc || $s == mc) { n++; if (($s == fc ? 2 : 1) != gsex[$idc]) d++ }
      END { print nm + 0, nf + 0, n + 0, d + 0 }' "${genofile}" FS="${fs}" "${PHENOFILE}")"   # whitespace for the genotype file, the phenotype file's own after
    if (( g_males == 0 || g_females == 0 || n_compared == 0 )); then
      warn+=("The genotype file (${genofile}) has no usable genetic sex (${g_males} male, ${g_females} female, ${n_compared} matched), so the codes cannot be checked against it. The phenotype file's 'sex' column decides who is fitted.")
    else
      pct=$(awk -v d="${n_discordant}" -v n="${n_compared}" 'BEGIN { printf "%.1f", 100 * d / n }')
      echo "--sex ${SEX}: genetic sex (${genofile}) agrees with the phenotype file's for $(( n_compared - n_discordant )) of ${n_compared} samples (${pct}% discordant)"
      if awk -v d="${n_discordant}" -v n="${n_compared}" 'BEGIN { exit !(d > n / 2) }'; then
        echo "REFUSED: under --femaleCode ${female} / --maleCode ${male}, ${pct}% of samples disagree with their genetic sex" >&2
        echo "(${genofile}, where 1 is male and 2 female): the codes are flipped. Set --femaleCode/--maleCode to the phenotype file's coding." >&2
        exit 1
      elif awk -v d="${n_discordant}" -v n="${n_compared}" 'BEGIN { exit !(d > n / 100) }'; then
        warn+=("${n_discordant} of ${n_compared} samples (${pct}%) have a 'sex' that disagrees with their genetic sex in ${genofile}: more than sample QC usually leaves. The phenotype file's 'sex' decides who is fitted.")
      fi
    fi
  fi
  if [[ ,${sex_values}, != ",${code},${other}," && ,${sex_values}, != ",${other},${code}," ]]; then
    warn+=("The 'sex' column holds ${sex_values}, not just ${female} (female) and ${male} (male): its coding may not be the one assumed.")
  fi
  # (quantitative only: a binary trait's other sex may be coded 0, controls, and its cases are the sharper test)
  if [[ ${TRAITTYPE} != "binary" ]] && (( n_dropped > n_kept )); then
    warn+=("More samples with a ${PHENOCOL} are of the sex being dropped (${n_dropped}) than kept (${n_kept}): are the codes flipped?")
  fi
  if [[ ${TRAITTYPE} == "binary" ]] && (( cases_dropped > 0 )); then
    warn+=("${cases_dropped} cases of ${PHENOCOL} have sex == ${other}, the sex being dropped: flipped codes, or a trait that is not sex-specific?")
  fi
  if (( ${#warn[@]} > 0 )); then
    {
      echo "################################################################################"
      echo "WARNING: CHECK THE SEX CODING for --sex ${SEX} (coded ${code}; female ${female}, male ${male})."
      for w in "${warn[@]}"; do echo "  - ${w}"; done
      echo "  --femaleCode/--maleCode default to 0 and 1, BRaVa's coding; SAIGE and the All of"
      echo "  Us file use 1 for female. Check the counts below before using this model."
      echo "################################################################################"
    } >&2
  fi
  echo "--sex ${SEX}: fitting the ${n_kept} samples with sex == ${code} and a non-missing ${PHENOCOL}; ${n_dropped} of the other sex are dropped${case_note}"
  if [[ ${SEX} == "F" ]]; then
    SEX_ARGS="--FemaleOnly=TRUE --sexCol=sex --FemaleCode=${code}"; SEX_SUFFIX="_FemaleOnly"
  else
    SEX_ARGS="--MaleOnly=TRUE --sexCol=sex --MaleCode=${code}"; SEX_SUFFIX="_MaleOnly"
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

# --useSparseGRMforVarRatio=TRUE: step 2 is given the sparse GRM, and a model
# whose variance-ratio file carries no `sparse` rows is REFUSED there (the
# All of Us models carry them). Without it every step-2 run fails on load.
# --tol is 0.02 (SAIGE's default) for BOTH trait types: the All of Us
# production runs used it for both, and the earlier 1e-5 for quantitative
# traits bought nothing measurable at a real cost in fit time.
TOL="0.02"
# Quantitative traits are ALWAYS inverse-rank normalised on import (a BRaVa
# requirement): --invNormalize=TRUE makes SAIGE replace the phenotype by
# qnorm((rank - 0.5) / n) on the final analysis set (after intersecting
# genotypes, GRM and complete covariates), before the covariates are fitted.
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
      --relatednessCutoff ${RELCUTOFF} \
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
      --isCovariateOffset TRUE \
      ${SEX_ARGS}"""

# A REFUSED fit (a gate, or any error after the output files were opened)
# leaves 0-byte <prefix>.rda and .varianceRatio.txt behind, which a later step
# 2 would load blindly and fail on uselessly. Remove them, so a refused trait
# leaves nothing that looks like a model; the log has the reason.
set +e
run_container; rc=$?
set -e
# a sex-specific fit's files carry SAIGE's _FemaleOnly/_MaleOnly suffix; step 2 expects <prefix>.rda
if [[ -n ${SEX_SUFFIX} ]]; then
  for ext in rda varianceRatio.txt; do
    [[ -f ${OUT}${SEX_SUFFIX}.${ext} ]] && mv -f "${OUT}${SEX_SUFFIX}.${ext}" "${OUT}.${ext}"
  done
fi
if (( rc != 0 )); then
  for f in "${OUT}.rda" "${OUT}.varianceRatio.txt"; do
    [[ -f ${f} && ! -s ${f} ]] && rm -f "${f}" && echo "removed empty ${f} (the fit was refused or failed; see the log)" >&2
  done
  exit ${rc}
fi
