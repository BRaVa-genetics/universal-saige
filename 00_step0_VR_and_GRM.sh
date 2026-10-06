#!/bin/bash

source ./run_container.sh

POSITIONAL_ARGS=()

SINGULARITY=false
generate_grm=false
generate_plink_for_vr=false
RELCUTOFF="0.05"

WD=$(pwd)
HOME=$WD

sample_lines(){   # $1 lines of stdin at random: a stand-in for shuf -n, which macOS lacks; seeded, so a rerun picks the same markers
    awk 'BEGIN {srand(20230116)} {printf "%.9f\t%s\n", rand(), $0}' | sort | head -n "$1" | cut -f2-
}

subset_variants(){
    echo "Subsetting genetic data for GRM / VR"

    # PLINK 1 is used as is; PLINK 2 is converted to .bed here, because the
    # merge, pruning and counts below are plink 1.9 (--max-alleles 2: a .bed
    # cannot hold a multiallelic variant, and none is needed for a GRM)
    if [[ $GENETIC_DATA_FORMAT == "pgen" ]]; then
        for file in $(ls ${GENETIC_DATA_DIR}/*.pgen); do
            prefix="${file%.pgen}"
            ./resources/plink2 --pfile "${prefix}" --max-alleles 2 --make-bed --out "${TMPD}/$(basename "${prefix}")"
            echo "${TMPD}/$(basename "${prefix}")" >> ${TMPD}/merge_list.txt
        done
    elif [[ $GENETIC_DATA_FORMAT == "plink" ]]; then
        for file in $(ls ${GENETIC_DATA_DIR}/*.bed); do
            echo "${file%.bed}" >> ${TMPD}/merge_list.txt
        done
    fi

    if [[ -n "$SAMPLEIDS" ]]; then
      # keep by IID alone: a PLINK 2 .psam usually has no FID (it is 0 in the
      # converted .fam), and a PLINK 1 FID need not equal the IID
      while read -r prefix; do cat "${prefix}.fam"; done < ${TMPD}/merge_list.txt \
        | awk 'NR == FNR {want[$1]; next} ($2 in want) {print $1, $2}' "$SAMPLEIDS" - \
        | sort -u > ${TMPD}/keep.txt
      ./resources/plink --merge-list ${TMPD}/merge_list.txt --make-bed --out ${TMPD}/merged --keep ${TMPD}/keep.txt
    else
      ./resources/plink --merge-list ${TMPD}/merge_list.txt --make-bed --out ${TMPD}/merged 
    fi

}

generate_GRM(){
    echo "LD pruning file for GRM generation"

    numRandomMarkerforSparseKin=5000

    ./resources/plink \
        --bfile "${TMPD}/merged" \
        --indep-pairwise 50 5 0.05 \
        --out "${TMPD}/merged"

    # Extract set of pruned variants and export to bfile
    ./resources/plink \
        --bfile "${TMPD}/merged" \
        --extract "${TMPD}/merged.prune.in" \
        --make-bed \
        --out "${HOME}/${OUT}.plink_for_grm"

    cmd="createSparseGRM.R \
        --plinkFile="${HOME}/${OUT}.plink_for_grm" \
        --nThreads=$(ncpu) \
        --outputPrefix="${HOME}/${OUT}" \
        --numRandomMarkerforSparseKin=$numRandomMarkerforSparseKin \
        --relatednessCutoff=${RELCUTOFF}"

    variant_count=$(wc -l < "${HOME}/${OUT}.plink_for_grm.bim")
    if [[ $variant_count -ge $numRandomMarkerforSparseKin ]]; then
      run_container
    else
      echo "Error: ${variant_count} variants found in ${OUT}.plink_for_grm, which is less than the required ${numRandomMarkerforSparseKin} variants."
      exit 1
    fi
    
    echo "GRM generated!"
}

generate_plink_for_vr(){
    # get count of variants in merged plink file:

    ./resources/plink \
        --bfile "${TMPD}/merged" \
        --freq counts \
        --out "${TMPD}/merged"

    variants_lessthan_20_MAC=2000
    variants_greaterthan_20_MAC=2000

    # .frq.counts columns: CHR SNP A1 A2 C1 C2 G0, where C1 and C2 are the two allele counts,
    # so the minor allele count is the smaller of C1 and C2
    cat <(
        tail -n +2 "${TMPD}/merged.frq.counts" \
        | awk '{mac = ($5 < $6) ? $5 : $6} mac >= 10 && mac < 20 {print $2}' \
        | sample_lines $variants_lessthan_20_MAC ) \
    <( \
        tail -n +2 "${TMPD}/merged.frq.counts" \
        | awk '{mac = ($5 < $6) ? $5 : $6} mac >= 20 {print $2}' \
        | sample_lines $variants_greaterthan_20_MAC \
        ) > "${TMPD}/merged.markerid.list"

    actual_variants_lessthan_20_MAC=$(awk 'NR > 1 {mac = ($5 < $6) ? $5 : $6} NR > 1 && mac >= 10 && mac < 20' "${TMPD}/merged.frq.counts" | wc -l)
    actual_variants_greaterthan_20_MAC=$(awk 'NR > 1 {mac = ($5 < $6) ? $5 : $6} NR > 1 && mac >= 20' "${TMPD}/merged.frq.counts" | wc -l)

    if [[ $variants_lessthan_20_MAC -gt $actual_variants_lessthan_20_MAC ]]; then
        echo "Error: ${actual_variants_lessthan_20_MAC} variants (MAC<20) found - less than the required ${variants_lessthan_20_MAC} variants."
        exit 1
    elif [[ $variants_greaterthan_20_MAC -gt $actual_variants_greaterthan_20_MAC ]]; then
        echo "Error: ${actual_variants_greaterthan_20_MAC} variants (MAC>20) found - less than the required ${variants_greaterthan_20_MAC} variants."
        exit 1
    fi

    # Extract markers from the large PLINK file
    ./resources/plink \
        --bfile "${TMPD}/merged" \
        --extract "${TMPD}/merged.markerid.list" \
        --make-bed \
        --out "${OUT}.plink_for_var_ratio"
}

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|--outputPrefix)
      OUT="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--isSingularity)
      SINGULARITY="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--geneticDataDirectory)
      GENETIC_DATA_DIR="$2"
      shift # past argument
      shift # past value
      ;;
    --geneticDataFormat)
      GENETIC_DATA_FORMAT="$2"
      shift # past argument
      shift # past value
      ;;
    --geneticDataType)
      GENETIC_DATA_TYPE="$2"
      shift # past argument
      shift # past value
      ;;
    --generate_GRM)
      generate_grm=true
      shift # past argument
      ;;
    --generate_plink_for_vr)
      generate_plink_for_vr=true
      shift # past argument
      ;;
    --relatednessCutoff)
      RELCUTOFF="$2"
      shift # past argument
      shift # past value
      ;;
    --sampleIDs)
      SAMPLEIDS="$2" 
      shift
      shift
      ;; 
    -h|--help)
      echo "usage: 00_step0_VR_and_GRM.sh
            required:
                --geneticDataDirectory: directory containing the genetic data (genotype/WES/WGS data in PLINK 1 or PLINK 2 format)
                --geneticDataFormat: format of the genetic data {plink,pgen}: PLINK 1 .bed/.bim/.fam or PLINK 2 .pgen/.pvar/.psam.
                  A VCF is refused; convert it once with plink2 (--vcf FILE --make-pgen --out PREFIX).
                --geneticDataType: type of the genetic data {WES,WGS,genotype}.
                -o,--outputPrefix: output prefix of the SAIGE step 0 output. This must be relative to, and contained within, the current working directory.
            optional:
                -s,--isSingularity (default: false): is singularity (or apptainer) available? If not, it is assumed that docker is available.
                --generate_GRM (default: false): generate GRM for the genetic data.
                --generate_plink_for_vr (default: false): generate plink file for vr.
                --relatednessCutoff (default: 0.05): GRM entries below it are dropped. Pass the SAME value to steps 1 and 2
                  (All of Us used 0.05, and 0.125 for its admixed amr cohort, whose GRM was too dense to fit at 0.05). The output
                  is <outputPrefix>_relatednessCutoff_<value>_5000_randomMarkersUsed.sparseGRM.mtx.
                --sampleIDs: path to a file containing sampleIDs (as a single column) to be used to define the GRM.
                Note that if nothing is passed, then all of the samples in the plink/pgen files will be used.
                Samples are matched on IID.
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

# check if either generate_GRM or generate_plink_for_vr:
if [[ ${generate_grm} = false ]] && [[ ${generate_plink_for_vr} = false ]]; then
  echo "Error: either generate_GRM or generate_plink_for_vr must be set to true"
  exit 1
fi

# PLINK 1 or PLINK 2 only, as in steps 1 and 2. Convert a VCF once with plink2.
if [[ ${GENETIC_DATA_FORMAT} != "plink" ]] && [[ ${GENETIC_DATA_FORMAT} != "pgen" ]]; then
  echo "geneticDataFormat must be in {plink,pgen}: PLINK 1 (.bed/.bim/.fam) or PLINK 2 (.pgen/.pvar/.psam)."
  echo "Convert a VCF once:  resources/plink2 --vcf FILE.vcf.gz --make-pgen --out PREFIX"
  exit 1
fi
if [[ ${GENETIC_DATA_FORMAT} == "pgen" && ! -x resources/plink2 ]]; then
  echo "resources/plink2 missing: bash download_resources.sh --plink2"
  exit 1
fi
if [[ ! -x resources/plink ]]; then
  echo "resources/plink (1.9) missing: bash download_resources.sh --plink"
  exit 1
fi

if [[ ${GENETIC_DATA_TYPE} != "WES" ]] && [[ ${GENETIC_DATA_TYPE} != "WGS" ]] && [[ ${GENETIC_DATA_TYPE} != "genotype" ]]; then
  echo "geneticDataType must be in {WES,WGS,genotype}"
  exit 1
fi

# check if genetic data directory exists and if files with the correct extension are present:
if [[ ! -d ${GENETIC_DATA_DIR} ]]; then
  echo "geneticDataDirectory does not exist"
  exit 1
fi

if [[ ${OUT:-} == "" ]]; then
  echo "--outputPrefix is required"
  exit 1
fi

check_relcutoff "${RELCUTOFF}"

echo "OUT               = ${OUT}"
echo "RELCUTOFF         = ${RELCUTOFF}"
echo "SINGULARITY       = ${SINGULARITY}"
echo "GENETIC DATA      = ${GENETIC_DATA_DIR}/*.{${GENETIC_DATA_FORMAT}}"
echo "SAMPLEIDS         = ${SAMPLEIDS}"

check_container_env $SINGULARITY

# Scratch space for this run only. A shared /tmp let one run read another's
# files (stale merge lists and conversions), and two cohorts run at once
# overwrote each other's merged data. Set TMPDIR to put it on a bigger disk.
TMPD=$(mktemp -d "${TMPDIR:-/tmp}/universal_saige_step0.XXXXXX")
trap 'rm -rf "${TMPD}"' EXIT

# For debugging
set -exo pipefail

## Set up directories
WD=$( pwd )

subset_variants

if [[ ${generate_grm} = true ]]; then
  echo "generating GRM"
  generate_GRM
fi

if [[ ${generate_plink_for_vr} = true ]]; then
  echo "generating plink for vr"
  generate_plink_for_vr
fi

# The GRM's density is the LAST thing step 0 prints, out of the trace, because
# this is where the cutoff is chosen: a dense GRM is fixed by rerunning step 0
# with a higher --relatednessCutoff, before any step-1 fit is paid for.
set +x
if [[ ${generate_grm} = true ]]; then
  echo
  grm_density_check "${HOME}/${OUT}_relatednessCutoff_${RELCUTOFF}_5000_randomMarkersUsed.sparseGRM.mtx" "${RELCUTOFF}" || true
fi
