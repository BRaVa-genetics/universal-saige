#!/bin/bash
# Build the FlexRV group file for one chromosome: the BRaVa group file plus a
# 'GENE score:NAME v1 v2 ...' line per gene, one value in [0, 1] per variant.
# FlexRV reads ONE score line per file, so one file per (chromosome, weight set):
#
#   AlphaMissense (the first weight set; Cheng et al. 2023):
#     bash 04_flexrv_groupfile.sh --group in/brava.chr20.txt --chr 20 --name AM --out in/brava.chr20.flexrv_AM.txt
#   a score carried as a column of the BRaVa long-form annotation table
#   (the consortium's second weight set arrives this way, or as an
#   AlphaMissense-shaped per-variant table passed to --am):
#     bash 04_flexrv_groupfile.sh --group in/brava.chr20.txt --annoTable in/brava.chr20.long.tsv.gz \
#          --scoreColumn MYSCORE --name MYSCORE --out in/brava.chr20.flexrv_MYSCORE.txt
#
# Then: 02_step2_SPAtests_variant_and_gene.sh --flexRVscore NAME --groupFile <out> ...
# The scoring is flexrv_score_from_alphamissense.py (standard-library Python 3;
# `python3 flexrv_score_from_alphamissense.py --help` for every option and
# `--selftest` for its controls). LoF variants score 1.0 by annotation;
# missense variants take the AlphaMissense pathogenicity joined on
# (CHROM, POS, REF, ALT); a missense variant with no score takes the gene's
# mean (--missing); every other variant a placeholder that is never tested.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
GROUP=""; OUT=""; CHR=""; NAME="AM"; AM="resources/AlphaMissense_hg38.tsv.gz"; ISO=""; ANNO=""; MANE=""; GTF=""; SCORECOL=""
LOF="pLoF"; MISSENSE="damaging_missense_or_protein_altering,other_missense_or_protein_altering"; DROP="non_coding"; EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --group)       GROUP="$2"; shift 2 ;;
    --out)         OUT="$2"; shift 2 ;;
    --chr)         CHR="${2#chr}"; shift 2 ;;
    --name)        NAME="$2"; shift 2 ;;
    --am)          AM="$2"; shift 2 ;;
    --isoforms)    ISO="$2"; shift 2 ;;
    --annoTable)   ANNO="$2"; shift 2 ;;
    --mane)        MANE="$2"; shift 2 ;;
    --gtf)         GTF="$2"; shift 2 ;;
    --scoreColumn) SCORECOL="$2"; shift 2 ;;
    --lofAnno)     LOF="$2"; shift 2 ;;
    --missenseAnno) MISSENSE="$2"; shift 2 ;;
    --dropAnno)    DROP="$2"; shift 2 ;;
    --missing)     EXTRA+=(--missing "$2"); shift 2 ;;
    -h|--help)     sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
[[ -n ${GROUP} && -n ${OUT} ]] || { echo "--group and --out are required" >&2; exit 1; }
args=(--group "${GROUP}" --out "${OUT}" --name "${NAME}" --lof-anno "${LOF}" --missense-anno "${MISSENSE}")
[[ -n ${DROP} ]] && args+=(--drop-anno "${DROP}")
if [[ -n ${SCORECOL} ]]; then
  [[ -n ${ANNO} ]] || { echo "--scoreColumn needs --annoTable" >&2; exit 1; }
  args+=(--anno-table "${ANNO}" --score-column "${SCORECOL}")
else
  [[ -s ${AM} ]] || { echo "AlphaMissense file not found: ${AM} (bash download_resources.sh --alphamissense)" >&2; exit 1; }
  args+=(--am "${AM}")
  [[ -n ${CHR} ]] && args+=(--chrom "${CHR}")
  [[ -n ${ISO} ]] && args+=(--isoforms "${ISO}")
  [[ -n ${ANNO} ]] && args+=(--anno-table "${ANNO}" --prefer-anno-transcript)
  # MANE summary: a fallback row on a transcript MANE puts in another gene is never used
  [[ -n ${MANE} ]] && args+=(--mane "${MANE}")
  # GENCODE GTF of the annotation's release (v39 for AoU v8): never a row on another gene's transcript
  [[ -n ${GTF} ]] && args+=(--gtf "${GTF}")
fi
# ${EXTRA[@]+...}: an empty array is "unbound" under set -u in bash < 4.4 (macOS ships 3.2)
python3 "${HERE}/flexrv_score_from_alphamissense.py" "${args[@]}" ${EXTRA[@]+"${EXTRA[@]}"}
