#!/bin/bash
# End-to-end test of every driver on a simulated cohort with a known answer
# (tests/simulate_cohort.py), then tests/check_results.py to score it.
#
#   bash tests/run_e2e.sh             # simulate (once), run every step, check
#   bash tests/run_e2e.sh --check     # re-score the last run only
#   bash tests/run_e2e.sh --fresh     # re-simulate as well
#   bash tests/run_e2e.sh --resume    # re-run only the runs that did not succeed (and every expected refusal)
#   SINGULARITY=true bash tests/run_e2e.sh   # the same, with every step on Singularity (resources/saige.sif)
#
# Needs Docker (or Singularity) and resources/ from
#   bash download_resources.sh --saige-image --plink --plink2 [--singularity]
# Writes tests/work/ (in/, out/, logs/, status.tsv). Takes 20-40 minutes (a group test 1-6 minutes). If it is
# much slower, check Docker Desktop's backend CPU: restarting Docker fixed a 10x slowdown once.
# Each run's log is tests/work/logs/<name>.log; status.tsv has its exit code.
set -uo pipefail
cd "$(dirname "$0")/.."

W=tests/work; IN=${W}/in; OUT=${W}/out; LOG=${W}/logs; STATUS=${W}/status.tsv
FRESH=false; CHECK_ONLY=false; RESUME=false; MAXJOBS=${MAXJOBS:-1}; SING=${SINGULARITY:-false}
for a in "$@"; do
  case "$a" in
    --fresh) FRESH=true ;;
    --check) CHECK_ONLY=true ;;
    --resume) RESUME=true ;;
    *) echo "unknown option: $a" >&2; exit 1 ;;
  esac
done
if [[ ${CHECK_ONLY} = true ]]; then
  exec python3 tests/check_results.py "${W}"
fi
if [[ ${SING} = true ]]; then
  IMAGE_FILES="resources/saige.sif"; GET="--singularity"
  command -v singularity > /dev/null || command -v apptainer > /dev/null \
    || { echo "SINGULARITY=true but neither singularity nor apptainer is on PATH" >&2; exit 1; }
else
  IMAGE_FILES="resources/saige.tar resources/saige.image"; GET=""
fi
for f in resources/plink resources/plink2 ${IMAGE_FILES}; do
  [[ -s ${f} ]] || { echo "${f} missing: bash download_resources.sh --saige-image --plink --plink2 ${GET}" >&2; exit 1; }
done
echo "container runtime: $([[ ${SING} = true ]] && echo Singularity || echo Docker)"

run () {   # run NAME ok|fail COMMAND...: stdin closed, output to the log, exit code recorded
  local name=$1 expect=$2; shift 2
  local t0=${SECONDS} rc
  if [[ ${RESUME} = true ]] && grep -q "^${name}"$'\t' "${STATUS}"; then
    printf "%-28s done in an earlier run\n" "${name}"
    return
  fi
  "$@" < /dev/null > "${LOG}/${name}.log" 2>&1; rc=$?
  printf "%s\t%s\t%d\t%d\n" "${name}" "${expect}" "${rc}" $(( SECONDS - t0 )) >> "${STATUS}"
  printf "%-28s exit %-3d (expected %-4s) %4ds\n" "${name}" "${rc}" "${expect}" $(( SECONDS - t0 ))
}

if [[ ${RESUME} = true ]]; then
  [[ ${FRESH} = false && -s ${STATUS} && -s ${IN}/truth.json ]] || { echo "--resume needs an earlier run (and no --fresh)" >&2; exit 1; }
  # keep the runs that succeeded; the rest go again, including the expected
  # refusals, whose non-zero exit alone does not show they failed for the right reason
  awk -F'\t' '$2 == "ok" && $3 == 0' "${STATUS}" > "${STATUS}.keep" && mv "${STATUS}.keep" "${STATUS}"
else
  if [[ ${FRESH} = true || ! -s ${IN}/truth.json ]]; then
    rm -rf "${IN}"
    python3 tests/simulate_cohort.py --out "${IN}" || exit 1
  fi
  rm -rf "${OUT}" "${LOG}"; : > "${STATUS}"
fi
mkdir -p "${OUT}" "${LOG}" "${IN}/array_pgen"

# ---- inputs as PLINK 2 too: the exome, and the array for step 0's PLINK 2 arm,
# with an IID-only .psam (the usual PLINK 2 layout: FID becomes 0 in a .bed)
resources/plink2 --bfile "${IN}/exome/chr7" --make-pgen --out "${IN}/exome/chr7" > /dev/null || exit 1
for c in 1 2 3 4; do
  resources/plink2 --bfile "${IN}/array/chr${c}" --make-pgen --out "${IN}/array_pgen/chr${c}" > /dev/null || exit 1
  awk 'BEGIN {FS = OFS = "\t"} NR == 1 {$2 = "#IID"} {$1 = ""; sub(/^\t/, ""); print}' \
    "${IN}/array_pgen/chr${c}.psam" > "${IN}/array_pgen/chr${c}.psam.tmp" && mv "${IN}/array_pgen/chr${c}.psam.tmp" "${IN}/array_pgen/chr${c}.psam"
done

GRM=${OUT}/step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx
GRMID=${GRM}.sampleIDs.txt
VR=${OUT}/step0.plink_for_var_ratio
COV="age,age2,sex,PC1,PC2,PC3,PC4"
MASKS="pLoF,damaging_missense_or_protein_altering,other_missense_or_protein_altering,synonymous,pLoF:damaging_missense_or_protein_altering,pLoF:damaging_missense_or_protein_altering:other_missense_or_protein_altering:synonymous"
S1=(bash 01_step1_fitNULLGLMM.sh --isSingularity "${SING}" --phenoFile "${IN}/pheno.tsv" --sparseGRM "${GRM}" --sparseGRMID "${GRMID}"
    --sampleIDs "${IN}/sample_ids.txt" --sampleIDCol IID)
S2=(bash 02_step2_SPAtests_variant_and_gene.sh --isSingularity "${SING}" --chr 7 --sparseGRM "${GRM}" --sparseGRMID "${GRMID}")

echo "== refusals and dry runs (no container)"
run s2_refuse_vcf          fail "${S2[@]}" --testType variant --vcf "${IN}/exome/chr7.vcf.gz" --modelFile m --varianceRatio v
run s2_refuse_both_geno    fail "${S2[@]}" --testType variant --pgen "${IN}/exome/chr7" --plink "${IN}/exome/chr7" --modelFile m --varianceRatio v
run s2_refuse_flex_masks   fail "${S2[@]}" --flexRVscore AM --annotations pLoF,synonymous --pgen "${IN}/exome/chr7" --groupFile g --modelFile m --varianceRatio v
run s2_refuse_flex_mafs    fail "${S2[@]}" --flexRVscore AM --flexRVmaxMAF 0.001,0.01 --pgen "${IN}/exome/chr7" --groupFile g --modelFile m --varianceRatio v
run s2_refuse_flex_variant fail "${S2[@]}" --flexRVscore AM --testType variant --pgen "${IN}/exome/chr7" --groupFile g --modelFile m --varianceRatio v
run s2_dryrun              ok   "${S2[@]}" --testType group --pgen "${IN}/exome/chr7" --groupFile g --annotations pLoF --modelFile m --varianceRatio v --dryRun
run s1_refuse_both_geno    fail "${S1[@]}" --traitType quantitative --phenoCol Q_pos --genotypePlink "${VR}" --genotypePgen "${VR}" --dryRun
run s1_sex_mf_ok           ok   "${S1[@]}" --traitType quantitative --phenoFile "${IN}/pheno_mf.tsv" --phenoCol Q_f --sex F --genotypePlink "${VR}" --dryRun
run s1_sex_mf_wrong        fail "${S1[@]}" --traitType quantitative --phenoFile "${IN}/pheno_mf.tsv" --phenoCol Q_f --sex M --genotypePlink "${VR}" --dryRun
run s1_sex_mixed           fail "${S1[@]}" --traitType quantitative --phenoCol Q_null --sex F --genotypePlink "${VR}" --dryRun
run s1_sex_numeric_ok      ok   "${S1[@]}" --traitType quantitative --phenoCol Q_female --sex F --genotypePlink "${VR}" --dryRun
run s0_refuse_vcf          fail bash 00_step0_VR_and_GRM.sh --isSingularity "${SING}" --geneticDataDirectory "${IN}/array" --geneticDataFormat vcf \
  --geneticDataType genotype --outputPrefix "${OUT}/never" --generate_GRM
run s0_refuse_no_out       fail bash 00_step0_VR_and_GRM.sh --isSingularity "${SING}" --geneticDataDirectory "${IN}/array" --geneticDataFormat plink \
  --geneticDataType genotype --generate_GRM

echo "== step 0"
run s0_plink ok bash 00_step0_VR_and_GRM.sh --isSingularity "${SING}" --geneticDataDirectory "${IN}/array" --geneticDataFormat plink \
  --geneticDataType genotype --outputPrefix "${OUT}/step0" --sampleIDs "${IN}/sample_ids.txt" --generate_GRM --generate_plink_for_vr
run s0_pgen ok bash 00_step0_VR_and_GRM.sh --isSingularity "${SING}" --geneticDataDirectory "${IN}/array_pgen" --geneticDataFormat pgen \
  --geneticDataType genotype --outputPrefix "${OUT}/step0_pgen" --sampleIDs "${IN}/sample_ids.txt" --generate_GRM
resources/plink --bfile "${VR}" --freq counts --out "${OUT}/vr_freq" > /dev/null
resources/plink2 --bfile "${VR}" --make-pgen --out "${VR}" > /dev/null

echo "== step 1"
run s1_Q_pos         ok   "${S1[@]}" --traitType quantitative --phenoCol Q_pos  --covarColList "${COV}" --categCovarColList sex --genotypePlink "${VR}" --outputPrefix "${OUT}/Q_pos"
run s1_Q_pos_pgen    ok   "${S1[@]}" --traitType quantitative --phenoCol Q_pos  --covarColList "${COV}" --categCovarColList sex --genotypePgen "${VR}" --outputPrefix "${OUT}/Q_pos_pgen"
run s1_Q_null        ok   "${S1[@]}" --traitType quantitative --phenoCol Q_null --covarColList "${COV}" --categCovarColList sex --genotypePlink "${VR}" --outputPrefix "${OUT}/Q_null"
run s1_Q_perm        ok   "${S1[@]}" --traitType quantitative --phenoCol Q_perm --covarColList "${COV}" --categCovarColList sex --genotypePlink "${VR}" --outputPrefix "${OUT}/Q_perm"
run s1_B_pos         ok   "${S1[@]}" --traitType binary --phenoCol B_pos  --covarColList "${COV}" --categCovarColList sex --genotypePlink "${VR}" --outputPrefix "${OUT}/B_pos"
run s1_B_rare_gated  fail "${S1[@]}" --traitType binary --phenoCol B_rare --covarColList "${COV}" --categCovarColList sex --genotypePlink "${VR}" --outputPrefix "${OUT}/B_rare_gated"
run s1_B_rare_ungated ok  env SAIGE_FIT_GATES=0 "${S1[@]}" --traitType binary --phenoCol B_rare --covarColList "${COV}" --categCovarColList sex --genotypePlink "${VR}" --outputPrefix "${OUT}/B_rare"
run s1_Q_female      ok   "${S1[@]}" --traitType quantitative --phenoCol Q_female --sex F --covarColList "age,age2,PC1,PC2,PC3,PC4" --genotypePlink "${VR}" --outputPrefix "${OUT}/Q_female"

echo "== step 4 (FlexRV group file)"
run s4_selftest ok python3 flexrv_score_from_alphamissense.py --selftest
run s4_flexrv   ok bash 04_flexrv_groupfile.sh --group "${IN}/group.chr7.txt" --chr 7 --name AM --am "${IN}/am_dummy.tsv.gz" \
  --out "${OUT}/group.chr7.flexrv_AM.txt"
# negative control for the score wiring: GENE_FLEX's scores shuffled among its variants
python3 - "${OUT}/group.chr7.flexrv_AM.txt" "${OUT}/group.chr7.flexrv_AM_shuffled.txt" << 'EOF'
import random, sys
rng = random.Random(1)
with open(sys.argv[1]) as fh, open(sys.argv[2], "w") as out:
    for line in fh:
        p = line.split()
        if p[0] == "GENE_FLEX" and p[1].startswith("score"):
            vals = p[2:]
            rng.shuffle(vals)
            line = " ".join(p[:2] + vals) + "\n"
        out.write(line)
EOF

echo "== step 2 (${MAXJOBS} at a time)"
# independent single-threaded runs (the image sets OMP, OpenBLAS and RcppParallel threads to 1), so
# MAXJOBS can go up to the core count; it is 1 because Docker Desktop on macOS bogged down under parallel runs
throttle () { while (( $(jobs -rp | wc -l) >= MAXJOBS )); do sleep 2; done; }
m () { echo --modelFile "${OUT}/$1.rda" --varianceRatio "${OUT}/$1.varianceRatio.txt"; }
throttle; run s2_Q_pos_variant     ok   "${S2[@]}" $(m Q_pos)  --testType variant --pgen "${IN}/exome/chr7" --outputPrefix "${OUT}/Q_pos.variant" &
throttle; run s2_Q_pos_group       ok   "${S2[@]}" $(m Q_pos)  --testType group --pgen "${IN}/exome/chr7" --groupFile "${IN}/group.chr7.txt" --annotations "${MASKS}" --outputPrefix "${OUT}/Q_pos.group" &
throttle; run s2_Q_pos_group_plink ok   "${S2[@]}" $(m Q_pos)  --testType group --plink "${IN}/exome/chr7" --groupFile "${IN}/group.chr7.txt" --annotations "${MASKS}" --outputPrefix "${OUT}/Q_pos.group_plink" &
throttle; run s2_Q_pos_flexrv      ok   "${S2[@]}" $(m Q_pos)  --flexRVscore AM --pgen "${IN}/exome/chr7" --groupFile "${OUT}/group.chr7.flexrv_AM.txt" --outputPrefix "${OUT}/Q_pos.flexrv_AM" &
throttle; run s2_Q_pos_flexrv_shuf ok   "${S2[@]}" $(m Q_pos)  --flexRVscore AM --pgen "${IN}/exome/chr7" --groupFile "${OUT}/group.chr7.flexrv_AM_shuffled.txt" --outputPrefix "${OUT}/Q_pos.flexrv_AM_shuffled" &
throttle; run s2_Q_null_variant    ok   "${S2[@]}" $(m Q_null) --testType variant --pgen "${IN}/exome/chr7" --outputPrefix "${OUT}/Q_null.variant" &
throttle; run s2_Q_null_group      ok   "${S2[@]}" $(m Q_null) --testType group --pgen "${IN}/exome/chr7" --groupFile "${IN}/group.chr7.txt" --annotations "${MASKS}" --outputPrefix "${OUT}/Q_null.group" &
throttle; run s2_Q_perm_group      ok   "${S2[@]}" $(m Q_perm) --testType group --pgen "${IN}/exome/chr7" --groupFile "${IN}/group.chr7.txt" --annotations "${MASKS}" --outputPrefix "${OUT}/Q_perm.group" &
throttle; run s2_B_pos_variant     ok   "${S2[@]}" $(m B_pos)  --testType variant --pgen "${IN}/exome/chr7" --outputPrefix "${OUT}/B_pos.variant" &
throttle; run s2_B_pos_group       ok   "${S2[@]}" $(m B_pos)  --testType group --pgen "${IN}/exome/chr7" --groupFile "${IN}/group.chr7.txt" --annotations "${MASKS}" --outputPrefix "${OUT}/B_pos.group" &
throttle; run s2_B_rare_gated      fail "${S2[@]}" $(m B_rare) --testType variant --pgen "${IN}/exome/chr7" --outputPrefix "${OUT}/B_rare.variant" &
wait

echo "== step 3 (Nglmm)"
run s3_nglmm ok bash 03_estimate_nGlmm.sh --isSingularity "${SING}" --contPhenos Q_pos --binaryPhenos "B_pos B_rare" --phenoFile "${IN}/pheno.tsv" \
  --covarList "${COV}" --sparseGRM "${GRM}" --sparseGRMID "${GRMID}" --outputFile "${OUT}/neff.csv"
# inputs outside the working directory, by absolute path through SAIGE_EXTRA_MOUNTS:
# the same Nglmm with the mount, "does not exist" without it
EXT=$(mktemp -d "${TMPDIR:-/tmp}/universal_saige_e2e.XXXXXX"); trap 'rm -rf "${EXT}"' EXIT
cp "${IN}/pheno.tsv" "${GRM}" "${GRMID}" "${EXT}/"
S3X=(bash 03_estimate_nGlmm.sh --isSingularity "${SING}" --phenoFile "${EXT}/pheno.tsv" --covarList "${COV}"
     --sparseGRM "${EXT}/${GRM##*/}" --sparseGRMID "${EXT}/${GRMID##*/}")
run s3_nglmm_extmount   ok   env SAIGE_EXTRA_MOUNTS="${EXT}" "${S3X[@]}" --contPhenos Q_pos --binaryPhenos "B_pos B_rare" \
  --outputFile "${OUT}/neff_extmount.csv"
# Docker sees only what is mounted, but apptainer/singularity also binds /tmp and
# the site's bind paths (apptainer.conf), where ${EXT} may well be: then the file
# is reachable without the mount and the negative control cannot fail. Ask the
# container itself, through run_container, and skip the control when it can see ${EXT}.
if ( source ./run_container.sh; WD=$(pwd); SINGULARITY=${SING}; SAIGE_EXTRA_MOUNTS=; DRYRUN=false
     cmd="test -d ${EXT}"; run_container ) < /dev/null > /dev/null 2>&1; then
  echo "${EXT}" > "${OUT}/s3_nglmm_nomount.skipped"
  printf "%-28s skipped: the container sees %s without a mount\n" s3_nglmm_nomount "${EXT}"
else
  run s3_nglmm_nomount  fail env SAIGE_EXTRA_MOUNTS= "${S3X[@]}" --contPhenos Q_pos --outputFile "${OUT}/neff_nomount.csv"
fi
run s3_refuse_bad_mount fail env SAIGE_EXTRA_MOUNTS=/nonexistent_universal_saige "${S3X[@]}" --contPhenos Q_pos \
  --outputFile "${OUT}/neff_badmount.csv" --dryRun

echo "== checks"
python3 tests/check_results.py "${W}"
