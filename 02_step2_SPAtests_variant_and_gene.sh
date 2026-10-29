#!/bin/bash
# SAIGE step 2: single-variant tests, SAIGE-GENE+ group tests, or FlexRV group
# tests, for one chromosome. Genotypes are PLINK 2 (.pgen/.pvar/.psam,
# recommended) or PLINK 1 (.bed/.bim/.fam). There is NO VCF input: the slim
# SAIGE image has no VCF reader, and converting on every run would charge every
# chromosome of every phenotype for a job that is done once -- convert first:
#   resources/plink2 --vcf exome.chr20.vcf.gz --make-pgen --out exome.chr20
# The choices below are the ones the All of Us production runs used
# (astheeggeggs/saige-slim, docs/state/aou-saige-parameters.md).

source ./run_container.sh

POSITIONAL_ARGS=()

SINGULARITY=false
OUT="out"
TESTTYPE=""
PLINK=""
PGEN=""
VCF=""
MODELFILE=""
VARIANCERATIO=""
SPARSEGRM=""
SPARSEGRMID=""
RELCUTOFF="0.05"
GROUPFILE=""
ANNOTATIONS=""
SUBSAMPLES=""
CONDITION=""
FLEXRV_SCORE=""
FLEXRV_MAXMAF="0.001"
FLEXRV_LOFANNO="pLoF"
FLEXRV_PARTIAL_MASK="false"
DRYRUN=false
LOG10P="TRUE"
NEFF="TRUE"
SINGLE_IN_GROUP=""

WD=$(pwd)
HOME=$WD

while [[ $# -gt 0 ]]; do
  case $1 in
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
      SINGULARITY="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--plink)
      PLINK="$2"
      shift # past argument
      shift # past value
      ;;
    --pgen)
      PGEN="$2"
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
    --relatednessCutoff)
      RELCUTOFF="$2"
      shift # past argument
      shift # past value
      ;;
    --condition)
      CONDITION="$2"
      shift # past argument
      shift # past value
      ;;
    --flexRVscore)
      FLEXRV_SCORE="$2"
      shift # past argument
      shift # past value
      ;;
    --flexRVmaxMAF)
      FLEXRV_MAXMAF="$2"
      shift # past argument
      shift # past value
      ;;
    --flexRVlofAnno)
      FLEXRV_LOFANNO="$2"
      shift # past argument
      shift # past value
      ;;
    --flexRVpartialMask)
      FLEXRV_PARTIAL_MASK="true"
      shift # past argument
      ;;
    --log10P)
      LOG10P="$2"
      shift # past argument
      shift # past value
      ;;
    --nEff)
      NEFF="$2"
      shift # past argument
      shift # past value
      ;;
    --singleInGroup)
      SINGLE_IN_GROUP="$2"
      shift # past argument
      shift # past value
      ;;
    --dryRun)
      DRYRUN=true
      shift # past argument
      ;;
    -h|--help)
      echo "usage: 02_step2_SPAtests_variant_and_gene.sh
  required:
    --testType: type of test {variant,group}.
    --pgen: plink 2 filename prefix of pgen/pvar/psam files (RECOMMENDED). This must be relative to, and contained within, the current working directory.
    -p,--plink: plink 1 filename prefix of bim/bed/fam files, the alternative to --pgen.
    --modelFile: filename of the model file output from step 1. This must be relative to, and contained within, the current working directory.
    --varianceRatio: filename of the varianceRatio file output from step 1. This must be relative to, and contained within, the current working directory.
    --sparseGRM: filename of the sparseGRM .mtx file. This must be relative to, and contained within, the current working directory.
    --sparseGRMID: filename of the sparseGRM ID file. This must be relative to, and contained within, the current working directory.
    --chr: chromosome to test (spelled as in the .pvar/.bim, e.g. chr20 or 20; a spelling the file does not use is refused).
  optional:
    -o,--outputPrefix: output prefix of the SAIGE step 2 output. The results are <prefix>.txt; group tests also write
      <prefix>.txt.singleAssoc.txt, .markerList.txt, .skatoMethod.txt (the p-value method of every SKAT-O cell), and the
      sidecars .pooledTests.txt, .skatFailures.txt, .spaFallbacks.txt and, from image 152ffd8c, .exactByWeight.txt (weighted
      cells one variant carries, given the exact test instead of the saddlepoint), .exactByWeightAboveCap.txt (such cells
      with too many carriers to enumerate, left on the saddlepoint), .spaPinned.txt (binary: tests whose saddlepoint left
      out samples the model holds certain) and .stretchGate.txt (binary weighted cells the variance ratio stretches; the
      most stretched report the exact convolution) and, from image b5456dee, .spaRegionNA.txt (burden cells whose
      variance-rescale saddlepoint did not converge or was non-finite, given the normal-approximation p; LEDGER #181)
      and, from image 66b80710, .sigmaGuard.txt (tests on the sparse-GRM band whose variance from the GRM solve was not
      positive or fell below half the no-GRM variance; they report the no-GRM variance instead; LEDGER #183).
      Each only when there is something to report; see the README.
    -s,--isSingularity (default: false): is singularity (or apptainer) available? If not, it is assumed that docker is available.
    -g,--groupFile: required if group test is selected. Filename of the annotation file used for group tests. This must be relative to, and contained within, the current working directory.
    --annotations: required if group test is selected. Comma separated list of annotation masks to test (':' joins labels INTO one mask, ',' separates masks). Please use
      'pLoF,damaging_missense,pLoF:damaging_missense,pLoF:damaging_missense:spliceai2_gt040,synonymous,other_missense,pLoF:damaging_missense:spliceai2_gt040:other_missense:synonymous'
      on the October 2026 BRaVa group files, or on older BRaVa files (other label names)
      'pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous'
      Every label (and every --flexRVlofAnno label) must be on an 'anno' line of the group file, or the run is refused.
    --relatednessCutoff (default: 0.05): MUST equal the cutoff steps 0, 1 and 3 used; nothing in SAIGE checks it.
    --condition: comma separated list of SNPs to condition on. This must be in order of the SNP occurrence in the dosage file.
    --subSampleFile: single-column file of sample IDs to restrict the test to.
    --log10P (default: TRUE): group tests also write log10 p-value columns (log10_Pvalue, log10_Pvalue_Burden,
      log10_Pvalue_SKAT; on FlexRV log10_Pvalue_Burden), computed in log space, so p-values below the double-precision
      floor (1e-308) are still estimated; the pooled Cauchy row is combined in log space too. They sit right after
      SE_Burden, so every later column moves: pass FALSE if downstream code reads columns by position.
    --nEff (default: TRUE): group tests also write nEff_Burden (after the log10 columns), the effective number of
      carriers behind each cell's burden test; a FlexRV hit with a small value (below ~10) rests on a few people.
      NA on the pooled Cauchy row. Image 49178058 or later.
    --singleInGroup (default: TRUE for SAIGE-GENE+, FALSE for FlexRV): write the per-variant results file
      <prefix>.txt.singleAssoc.txt. Off for FlexRV and for the am_moderate run, as in All of Us: their variants are in
      the main SAIGE-GENE+ run's file already, and on FlexRV it is most of the run time.
    --dryRun: print the SAIGE command instead of running it.
  FlexRV (Schwartzentruber et al. 2025; a group test, burden statistic, 16 score x up to 12 MAF transforms per gene):
    --flexRVscore NAME: run FlexRV on the group file's 'score:NAME' line, which every region must have, or the run is
      refused (see 04_flexrv_groupfile.sh). Implies --testType group,
      ONE annotation mask (default 'pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering', the older
      BRaVa labels; on the October 2026 files pass --annotations pLoF:damaging_missense:other_missense:spliceai2_gt040),
      ONE max MAF (--flexRVmaxMAF, default 0.001) and r.corr = 1.
    --flexRVlofAnno (default: pLoF): the label(s) the 'lof' score transform keys on; must be the labels used when the score line was built.
    --flexRVpartialMask: allow a mask that leaves out variants carrying a score > 0 on the score line. Without it such a run is
      refused, naming the labels left out: SAIGE drops them before reading their scores (e.g. missense variants labelled with
      the SpliceAI category when the mask lists only pLoF and the missense labels). Pass it only for a deliberate subset.
  a VCF is refused: convert once with plink2 (--vcf FILE --make-pgen --out PREFIX) and pass --pgen.
      "
      exit 0
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
if [[ ${VCF} != "" ]]; then
  echo "ERROR: --vcf is not accepted. The SAIGE image reads PLINK 2 (.pgen/.pvar/.psam) or PLINK 1 (.bed/.bim/.fam) only,"
  echo "and a conversion on every run would be paid once per chromosome per phenotype. Convert ONCE, then pass --pgen:"
  echo "    resources/plink2 --vcf ${VCF} --make-pgen --out <prefix>      (bash download_resources.sh --plink2 fetches plink2)"
  exit 1
fi

for f in "--log10P:${LOG10P}" "--nEff:${NEFF}" "--singleInGroup:${SINGLE_IN_GROUP:-TRUE}"; do
  if [[ ${f#*:} != "TRUE" && ${f#*:} != "FALSE" ]]; then
    echo "${f%%:*} takes TRUE or FALSE, not '${f#*:}'"
    exit 1
  fi
done
# Per-variant results of a group test: on for SAIGE-GENE+, off for FlexRV, as in the All of Us runs (2026-10-08): on a
# FlexRV run the file was 80% of the time and repeats the SAIGE-GENE+ run's per-variant results.
[[ -z ${SINGLE_IN_GROUP} ]] && { [[ ${FLEXRV_SCORE} != "" ]] && SINGLE_IN_GROUP="FALSE" || SINGLE_IN_GROUP="TRUE"; }
if [[ ${FLEXRV_SCORE} != "" ]]; then
  if [[ ${TESTTYPE} != "" && ${TESTTYPE} != "group" ]]; then
    echo "--flexRVscore is a group test; --testType ${TESTTYPE} contradicts it"
    exit 1
  fi
  TESTTYPE="group"
  [[ ${ANNOTATIONS} == "" ]] && ANNOTATIONS="pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering"
  if [[ ${ANNOTATIONS} == *,* ]]; then
    echo "FlexRV tests ONE annotation mask (labels joined with ':'); --annotations '${ANNOTATIONS}' names several"
    exit 1
  fi
  if [[ ${FLEXRV_MAXMAF} == *,* ]]; then
    echo "FlexRV tests ONE max MAF; --flexRVmaxMAF '${FLEXRV_MAXMAF}' names several"
    exit 1
  fi
fi

if [[ ${TESTTYPE} == "" ]]; then
  echo "Test type not set"
  exit 1
fi

if [[ ${PLINK} == "" ]] && [[ ${PGEN} == "" ]]; then
  echo "genotypes not set: pass --pgen <prefix> (recommended) or --plink <prefix>"
  exit 1
fi

if [[ ${PLINK} != "" ]] && [[ ${PGEN} != "" ]]; then
  echo "pass ONE of --pgen and --plink"
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
  echo "variance ratio file not set"
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

# Every label the masks name (and, for FlexRV, every --flexRVlofAnno label) must
# be on an `anno` line of the group file. SAIGE checks nothing of the kind for
# the masks, and only warns for the lof labels (saige-slim LEDGER #170), so a
# typo, or BRaVa labels against a group file written with other ones (AoU's say
# damaging_missense where BRaVa's say damaging_missense_or_protein_altering),
# would test a smaller mask, or an empty one, under the name asked for.
if [[ ${TESTTYPE} == "group" ]]; then
  if [[ ! -s ${GROUPFILE} ]]; then
    echo "ERROR: group file '${GROUPFILE}' not found or empty"
    exit 1
  fi
  case ${GROUPFILE} in *.gz|*.bgz) read_group="gzip -cd" ;; *) read_group="cat" ;; esac
  asked="${ANNOTATIONS}"
  [[ ${FLEXRV_SCORE} != "" ]] && asked="${asked},${FLEXRV_LOFANNO}"
  # line 1: the labels asked for that the file lacks; line 2: the labels it has;
  # line 3: regions (var lines) and, for FlexRV, score:<NAME> lines
  labels=$(${read_group} "${GROUPFILE}" | awk -v asked="${asked}" -v score="score:${FLEXRV_SCORE}" '
    $2 == "var" { nvar++ }
    $2 == score { nscore++ }
    $2 == "anno" { for (i = 3; i <= NF; i++) if (!($i in seen)) { seen[$i]; order[++n] = $i } }
    END {
      m = split(asked, a, /[,:;]/)
      for (i = 1; i <= m; i++) if (a[i] != "" && !(a[i] in seen) && !(a[i] in told)) { told[a[i]]; miss = miss " " a[i] }
      for (i = 1; i <= n; i++) have = have " " order[i]
      print miss; print have; print nvar + 0, nscore + 0
    }')
  missing=$(echo "${labels}" | sed -n 1p)
  if [[ -n ${missing# } ]]; then
    echo "ERROR: annotation label(s) asked for but on no variant in the group file:${missing}"
    echo "  group file:         ${GROUPFILE}"
    echo "  labels it has:     $(echo "${labels}" | sed -n 2p)"
    echo "  --annotations:      ${ANNOTATIONS}"
    [[ ${FLEXRV_SCORE} != "" ]] && echo "  --flexRVlofAnno:    ${FLEXRV_LOFANNO}"
    if [[ ${FLEXRV_SCORE} != "" ]]; then
      echo "Pass --annotations (ONE mask, labels joined with ':') and --flexRVlofAnno with labels the file has to the FlexRV"
      echo "call (in templates/step_2_template.sh: section 3, which uses the defaults), or pass the group file they were written for."
    else
      echo "Change --annotations in your step-2 call to labels the file has (in templates/step_2_template.sh: the 'annots'"
      echo "line), or pass the group file those labels were written for."
    fi
    exit 1
  fi
  # FlexRV reads its weights from a score:<NAME> line, which every region needs
  # (one value per variant on its var line; 04_flexrv_groupfile.sh writes them).
  # The AoU step-2 wrapper refuses a file that lacks them before any compute.
  if [[ ${FLEXRV_SCORE} != "" ]]; then
    read -r nvar nscore <<< "$(echo "${labels}" | sed -n 3p)"
    if (( nscore != nvar )); then
      echo "ERROR: --flexRVscore ${FLEXRV_SCORE} needs one 'score:${FLEXRV_SCORE}' line per region; ${GROUPFILE}"
      echo "  has ${nscore} for ${nvar} region(s). Build the file with 04_flexrv_groupfile.sh --name ${FLEXRV_SCORE}, or pass the"
      echo "  name of the score line the file does carry."
      exit 1
    fi
    # The mask must hold every variant the score line weights: SAIGE drops markers matching no requested
    # label before reading any score, so a mask that omits a scored label tests fewer variants under the
    # same name, silently. Scores > 0 mark the variants meant to be tested (04_flexrv_groupfile.sh writes
    # --other, default 0, for the rest). Found in the BRaVa 2026-10 rebuild: missense variants labelled with
    # the SpliceAI category carried AlphaMissense scores the mask pLoF:damaging_missense:other_missense
    # never reached. The score line follows its region's anno line.
    if [[ ${FLEXRV_PARTIAL_MASK} != "true" ]]; then
      outside=$(${read_group} "${GROUPFILE}" | awk -v mask="${ANNOTATIONS}" -v score="score:${FLEXRV_SCORE}" '
        BEGIN { n = split(mask, a, /[:;]/); for (i = 1; i <= n; i++) inmask[a[i]] }
        $2 == "anno" { split("", lab); for (i = 3; i <= NF; i++) lab[i] = $i }
        $2 == score { for (i = 3; i <= NF; i++) if ($i + 0 > 0 && !(lab[i] in inmask)) { if (!(lab[i] in cnt)) order[++k] = lab[i]; cnt[lab[i]]++ } }
        END { for (j = 1; j <= k; j++) printf " %s (%d)", order[j], cnt[order[j]] }')
      if [[ -n ${outside} ]]; then
        echo "ERROR: the FlexRV mask leaves out variants that carry a score on the score:${FLEXRV_SCORE} line:${outside}"
        echo "  (label, number of variants with score > 0). SAIGE drops them before reading their scores."
        echo "  --annotations:      ${ANNOTATIONS}"
        echo "  group file:         ${GROUPFILE}"
        echo "Add the label(s) to the mask (labels joined with ':'), or pass --flexRVpartialMask if leaving them out is"
        echo "intended (a deliberate subset of the scored variants)."
        exit 1
      fi
    fi
  fi
fi

if [[ ${CHR:-} == "" ]]; then
  echo "--chr not set"
  exit 1
fi

# --chr must be spelled as column 1 of the .pvar/.bim spells it ("20" vs "chr20"):
# SAIGE matches it exactly, and a mismatch is not an error there. The AoU
# step-2 wrapper refuses it before any compute; so does this.
if [[ ${PGEN} != "" ]]; then varfile="${PGEN}.pvar"; else varfile="${PLINK}.bim"; fi
if [[ -r ${varfile} ]] && ! awk -v c="${CHR}" '!/^#/ && $1 == c { found = 1; exit } END { exit !found }' "${varfile}"; then
  echo "ERROR: --chr ${CHR} is not in column 1 of ${varfile}, which has: $(awk '!/^#/ { if (!($1 in s)) { s[$1]; printf "%s ", $1; if (++n == 5) exit } }' "${varfile}")"
  echo "  Pass --chr spelled as the file spells it (e.g. '20' or 'chr20')."
  exit 1
fi

if [[ $SUBSAMPLES != "" ]]; then
  SUBSAMPLES="${HOME}/${SUBSAMPLES}"
fi

if [[ $OUT = "out" ]]; then
  echo "Warning: outputPrefix not set, setting outputPrefix to 'out'. Check that this will not overwrite existing files."
fi

echo "OUT               = ${OUT}"
echo "TESTTYPE          = ${TESTTYPE}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "PGEN              = ${PGEN:+${PGEN}.{pgen/pvar/psam}}"
echo "PLINK             = ${PLINK:+${PLINK}.{bim/bed/fam}}"
echo "MODELFILE         = ${MODELFILE}"
echo "VARIANCERATIO     = ${VARIANCERATIO}"
echo "GROUPFILE         = ${GROUPFILE}"
echo "ANNOTATIONS       = ${ANNOTATIONS}"
echo "SPARSEGRM         = ${SPARSEGRM}"
echo "SPARSEGRMID       = ${SPARSEGRMID}"
echo "RELCUTOFF         = ${RELCUTOFF}"
echo "CONDITION         = ${CONDITION}"
echo "FLEXRV            = ${FLEXRV_SCORE:-off}${FLEXRV_SCORE:+ (maxMAF ${FLEXRV_MAXMAF}, lof labels ${FLEXRV_LOFANNO})}"
echo "LOG10P            = ${LOG10P}"
echo "NEFF              = ${NEFF}"
echo "SINGLE_IN_GROUP   = ${SINGLE_IN_GROUP}"

# For debugging
set -eo pipefail

## Set up directories
WD=$( pwd )

if [[ "$TESTTYPE" = "variant" ]]; then
  echo "variant testing"
  min_mac="4"
  GROUPFILE=""
  GROUP_ARGS=""
else
  echo "gene testing"
  min_mac="0.5"
  GROUPFILE="${HOME}/${GROUPFILE}"
  # log10 p-value columns (saige-slim image b5456dee and later; on by default since 2026-10-16, as in the AoU runs):
  # p-values past the double floor are estimated, not printed as 0 (AoU: GPT for ALT, ~1e-146 on its own and
  # beyond the floor in the pooled row)
  GROUP_ARGS="--groupFile=${GROUPFILE} --annotation_in_groupTest=${ANNOTATIONS} --is_output_markerList_in_groupTest=TRUE --is_single_in_groupTest=${SINGLE_IN_GROUP} --is_output_log10P=${LOG10P} --is_output_nEff=${NEFF}"
  if [[ ${FLEXRV_SCORE} != "" ]]; then
    GROUP_ARGS="${GROUP_ARGS} --maxMAF_in_groupTest=${FLEXRV_MAXMAF} --flexRV_score=${FLEXRV_SCORE} --flexRV_maxMAF=${FLEXRV_MAXMAF} --flexRV_lofAnno=${FLEXRV_LOFANNO} --r.corr=1"
  else
    GROUP_ARGS="${GROUP_ARGS} --maxMAF_in_groupTest=0.0001,0.001,0.01"
  fi
fi

if [[ ${PGEN} != "" ]]; then
  GENO_ARGS="--pgenFile=${HOME}/${PGEN}.pgen --pvarFile=${HOME}/${PGEN}.pvar --psamFile=${HOME}/${PGEN}.psam"
else
  GENO_ARGS="--bedFile=${HOME}/${PLINK}.bed --bimFile=${HOME}/${PLINK}.bim --famFile=${HOME}/${PLINK}.fam"
fi

# Firth off and fastTest off are the All of Us production choices: fastTest is
# a two-stage screen on the variant path and inert on the region path, and
# Firth was most of the CPU of a binary scan for effect sizes the tests do not
# use. Missingness and imputation are the build's defaults (0.15, best_guess).
cmd="step2_SPAtests.R \
        ${GENO_ARGS} \
        ${GROUP_ARGS} \
        --chrom=${CHR} \
        --minMAF=0 \
        --minMAC=${min_mac} \
        --GMMATmodelFile=${HOME}/${MODELFILE} \
        --varianceRatioFile=${HOME}/${VARIANCERATIO} \
        --sparseGRMFile=${HOME}/${SPARSEGRM} \
        --sparseGRMSampleIDFile=${HOME}/${SPARSEGRMID} \
        --relatednessCutoff=${RELCUTOFF} \
        --LOCO=FALSE \
        --is_Firth_beta=FALSE \
        --is_fastTest=FALSE \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${HOME}/${OUT}.txt"
[[ ${SUBSAMPLES} != "" ]] && cmd="${cmd} --subSampleFile=${SUBSAMPLES}"
[[ ${CONDITION}  != "" ]] && cmd="${cmd} --condition=${CONDITION}"

run_container
