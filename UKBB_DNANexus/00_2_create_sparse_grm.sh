#!/bin/bash 
#
# Based on https://saigegit.github.io/SAIGE-doc/docs/UK_Biobank_WES_analysis.html
#

set -e # Stop job if any command fails

readonly bfile=$1
readonly out=$2

createSparseGRM.R       \
    --plinkFile="${bfile}" \
    --nThreads=72  \
    --outputPrefix="${out}"       \
    --numRandomMarkerforSparseKin=5000      \
    --relatednessCutoff=0.05