#!/bin/bash
# Fetch the tools universal-saige needs. Nothing here is required more than once.
#
#   bash download_resources.sh --saige-image --plink2            # Docker
#   bash download_resources.sh --saige-image --plink2 --singularity
#   bash download_resources.sh --plink                            # plink 1.9, step 0 only
#   bash download_resources.sh --alphamissense [--isoforms]       # FlexRV weights (AlphaMissense release)
#
# THE IMAGE is the slim SAIGE build (https://github.com/astheeggeggs/saige-slim),
# pulled from Docker Hub, which is free to pull from anywhere. It reads PLINK 1
# (.bed/.bim/.fam) and PLINK 2 (.pgen/.pvar/.psam) genotypes and NOTHING ELSE:
# there is no VCF reader. Convert once with plink2 (see the README) -- PLINK 2
# format is the recommended input for every step.
set -euo pipefail

SAIGE_IMAGE="${SAIGE_IMAGE:-astheeggeggs/saige-slim}"        # Docker Hub namespace/repository
SAIGE_VERSION="${SAIGE_VERSION:-1.5.2-dev-3e92d89d}"          # the tag; pin it, and record it with your results
PLINK2_DATE="20260919"                                        # plink2 alpha 6 build, all platforms
AM_RECORD="https://zenodo.org/records/8208688/files"          # AlphaMissense release (Cheng et al. 2023)

GET_IMAGE=false; GET_PLINK=false; GET_PLINK2=false; GET_AM=false; GET_ISO=false; SINGULARITY=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --saige-image)  GET_IMAGE=true; shift ;;
    --singularity)  SINGULARITY=true; shift ;;
    --plink)        GET_PLINK=true; shift ;;
    --plink2)       GET_PLINK2=true; shift ;;
    --alphamissense) GET_AM=true; shift ;;
    --isoforms)     GET_ISO=true; shift ;;
    -h|--help)      sed -n '2,15p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
mkdir -p resources/
machine=$(uname); arch=$(uname -m)

if [[ ${GET_IMAGE} = true ]]; then
  ref="${SAIGE_IMAGE}:${SAIGE_VERSION}"
  if [[ ${SINGULARITY} = true ]]; then
    # apptainer is singularity's renamed successor; either pulls the same .sif
    sing=$(command -v singularity || command -v apptainer) || { echo "neither singularity nor apptainer found on PATH" >&2; exit 1; }
    if [[ ! -s resources/saige.sif ]]; then
      # The pull downloads and unpacks the image's layers. By default that is in
      # /tmp plus a cache in $HOME, both often small on clusters; here it is a
      # scratch directory next to the .sif (or APPTAINER_TMPDIR/SINGULARITY_TMPDIR
      # if set), removed afterwards, with no cache. Running the .sif needs neither.
      pull_tmp=$(mktemp -d "${PWD}/resources/.pull_tmp.XXXXXX")
      trap 'rm -rf "${pull_tmp}"' EXIT
      tmp="${APPTAINER_TMPDIR:-${SINGULARITY_TMPDIR:-${pull_tmp}}}"
      APPTAINER_TMPDIR="${tmp}" SINGULARITY_TMPDIR="${tmp}" \
        "${sing}" pull --disable-cache "resources/saige.sif" "docker://${ref}"
    fi
  else
    docker pull "${ref}"
    docker save -o "resources/saige.tar" "${ref}"
  fi
  echo "${ref}" > resources/saige.image
  echo "SAIGE image: ${ref}"
fi

if [[ ${GET_PLINK2} = true ]]; then
  # plink2 is what converts a VCF, subsets samples/variants and writes .pgen;
  # the linux x86_64 build runs on any x86 CPU (the avx2 one is faster, if yours has it)
  case "${machine}-${arch}" in
    Linux-x86_64)  f="plink2_linux_x86_64_${PLINK2_DATE}.zip" ;;
    Darwin-arm64)  f="plink2_mac_arm64_${PLINK2_DATE}.zip" ;;
    Darwin-*)      f="plink2_mac_${PLINK2_DATE}.zip" ;;
    *) echo "no plink2 build known for ${machine}-${arch}; see https://www.cog-genomics.org/plink/2.0/" >&2; exit 1 ;;
  esac
  wget -nc "https://s3.amazonaws.com/plink2-assets/alpha6/${f}" -P resources/
  unzip -o "resources/${f}" plink2 -d resources/ >/dev/null
  chmod a+x resources/plink2 && echo "plink2: $(resources/plink2 --version 2>/dev/null | head -1)"
fi

if [[ ${GET_PLINK} = true ]]; then
  # plink 1.9: used by step 0 to merge, LD-prune and count (its flags are 1.9 syntax)
  case "${machine}" in
    Darwin) f="plink_mac_20230116.zip" ;;
    Linux)  f="plink_linux_x86_64_20230116.zip" ;;
    *) echo "operating system not supported" >&2; exit 1 ;;
  esac
  wget -nc "https://s3.amazonaws.com/plink1-assets/${f}" --no-check-certificate -P resources/
  unzip -o "resources/${f}" -d resources/ >/dev/null && echo "plink 1.9 in resources/"
fi

if [[ ${GET_AM} = true ]]; then
  # canonical transcripts (~600 MB); --isoforms adds the non-canonical complement (~1.2 GB),
  # which the release marks as less evaluated -- only 04_flexrv_groupfile.sh --isoforms uses it
  wget -nc "${AM_RECORD}/AlphaMissense_hg38.tsv.gz" -P resources/
  [[ ${GET_ISO} = true ]] && wget -nc "${AM_RECORD}/AlphaMissense_isoforms_hg38.tsv.gz" -P resources/
  echo "AlphaMissense in resources/"
fi
