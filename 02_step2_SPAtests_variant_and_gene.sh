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
DRYRUN=false

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
    --chr: chromosome to test (spelled as in the .pvar/.bim, e.g. chr20 or 20).
  optional:
    -o,--outputPrefix: output prefix of the SAIGE step 2 output. The results are <prefix>.txt; group tests also write
      <prefix>.txt.singleAssoc.txt, .markerList.txt, .skatoMethod.txt (the p-value method of every SKAT-O cell), and the
      sidecars .pooledTests.txt, .skatFailures.txt, .spaFallbacks.txt (each only when there is something to report).
    -s,--isSingularity (default: false): is singularity available? If not, it is assumed that docker is available.
    -g,--groupFile: required if group test is selected. Filename of the annotation file used for group tests. This must be relative to, and contained within, the current working directory.
    --annotations: required if group test is selected. Comma separated list of annotation masks to test (':' joins labels INTO one mask, ',' separates masks). Please use
      'pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous'
    --relatednessCutoff (default: 0.05): MUST equal the cutoff step 1 fitted under; nothing in SAIGE checks it.
    --condition: comma separated list of SNPs to condition on. This must be in order of the SNP occurrence in the dosage file.
    --subSampleFile: single-column file of sample IDs to restrict the test to.
    --dryRun: print the SAIGE command instead of running it.
  FlexRV (Schwartzentruber et al. 2025; a group test, burden statistic, 16 score x up to 12 MAF transforms per gene):
    --flexRVscore NAME: run FlexRV on the group file's 'score:NAME' line (see 04_flexrv_groupfile.sh). Implies --testType group,
      ONE annotation mask (default 'pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering'),
      ONE max MAF (--flexRVmaxMAF, default 0.001) and r.corr = 1.
    --flexRVlofAnno (default: pLoF): the label(s) the 'lof' score transform keys on; must be the labels used when the score line was built.
  a VCF is refused: convert once with plink2 (--vcf FILE --make-pgen --out PREFIX) and pass --pgen.
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
if [[ ${VCF} != "" ]]; then
  echo "ERROR: --vcf is not accepted. The SAIGE image reads PLINK 2 (.pgen/.pvar/.psam) or PLINK 1 (.bed/.bim/.fam) only,"
  echo "and a conversion on every run would be paid once per chromosome per phenotype. Convert ONCE, then pass --pgen:"
  echo "    resources/plink2 --vcf ${VCF} --make-pgen --out <prefix>      (bash download_resources.sh --plink2 fetches plink2)"
  exit 1
fi

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

if [[ ${CHR:-} == "" ]]; then
  echo "--chr not set"
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
  GROUP_ARGS="--groupFile=${GROUPFILE} --annotation_in_groupTest=${ANNOTATIONS} --is_output_markerList_in_groupTest=TRUE --is_single_in_groupTest=TRUE"
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
