# Run $cmd inside the SAIGE image, with the working directory mounted at the
# same path (every file passed to a step must be under it).
#
# The slim build's FIT GATES are on: a binary trait fitted on fewer than 100
# cases, a categorical covariate level with fewer than 10 cases or controls, a
# separated covariate model, or a fit that did not converge, is REFUSED with a
# message that names the gate (step 1 checks the fit; step 2 re-checks the
# model it loads). That refusal is correct: the run would not have been a
# valid analysis. SAIGE_FIT_GATES=0 in your environment turns them into
# warnings -- not recommended, and the log then says so on every run.
#
# Directories outside it can be bound too, read-only and at the same path: list
# them, colon-separated, in SAIGE_EXTRA_MOUNTS. Only step 3 takes absolute paths
# into them; steps 0-2 read every file relative to the working directory.
#
# DRYRUN=true prints the command instead of running it.
ncpu () {   # portable core count (nproc is Linux-only)
  nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1
}
singularity_bin () {   # singularity, or apptainer (its renamed successor; many clusters now ship only apptainer)
  command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null
}
check_relcutoff () {   # $1 = a --relatednessCutoff value; stops unless it is a number in (0, 1)
  # The same cutoff must reach steps 0, 1 and 2: step 0 keeps GRM entries at or
  # above it, and steps 1 and 2 thin the GRM they load to it again. Nothing in
  # SAIGE checks that they agree (All of Us: 0.05, but 0.125 for amr, whose
  # admixed GRM had ~46M entries at 0.05 and a fit that never finished).
  awk -v x="$1" 'BEGIN { exit !(x ~ /^(0?\.[0-9]+|0)$/ && x + 0 > 0) }' \
    || { echo "--relatednessCutoff must be a number between 0 and 1 (got '$1')" >&2; exit 1; }
}
grm_density_check () {   # $1 = sparse GRM .mtx, $2 = the cutoff it is used at; warns LOUDLY when it is dense
  # Mean relatives per sample (2 x off-diagonal pairs at or above the cutoff / N),
  # counted as SAIGE will use the matrix (it drops entries below the cutoff). A
  # GRM from 5,000 random markers is noisy, and in an admixed cohort ancestry
  # also reads as relatedness, so a low cutoff can flood it. All of Us amr
  # (N 52,192): ~885 relatives per sample at 0.05 (46.2M entries), a step-1 fit
  # that never finished and tau collapsed to 0 in ~95% of models; at 0.125, ~5.7
  # and a fit in 58 s (saige-slim handoff 2026-09-25). Real families give a few
  # (the e2e pedigree ~2.4). The warning threshold, 100, is a judgement between.
  local mtx=$1 cutoff=$2 stats n pairs rel
  [[ -r ${mtx} ]] || return 0
  stats=$(awk -v c="${cutoff}" '
    NR == 1 { sym = ($0 ~ /symmetric/) }
    /^%/ { next }
    !dims { n = $1; dims = 1; next }
    $1 != $2 && $3 >= c { pairs++ }
    END { if (!sym) pairs /= 2; printf "%d %d %.1f\n", n, pairs, n ? 2 * pairs / n : 0 }' "${mtx}")
  read -r n pairs rel <<< "${stats}"
  echo "sparse GRM at --relatednessCutoff ${cutoff}: ${n} samples, ${pairs} related pairs, ${rel} relatives per sample"
  if awk -v r="${rel}" 'BEGIN { exit !(r > 100) }'; then
    {
      echo "################################################################################"
      echo "WARNING: THE SPARSE GRM IS DENSE: ${rel} relatives per sample on average"
      echo "  (${pairs} pairs among ${n} samples at --relatednessCutoff ${cutoff})."
      echo "  At this density the step-1 fit can run for hours or never finish, and the"
      echo "  random effect can collapse to 0. Most of these 'relatives' are noise or"
      echo "  shared ancestry, not family. RAISE --relatednessCutoff -- All of Us used"
      echo "  0.125 for its admixed amr cohort: ~885 relatives per sample at 0.05, ~5.7"
      echo "  at 0.125, and a fit that never finished took 58 s -- and pass the SAME"
      echo "  value to steps 0, 1 and 2."
      echo "################################################################################"
    } >&2
  fi
  return 0
}
check_container_env () {   # $1 = true for singularity; stops early when the runtime or the image is missing
  if [[ ${1:-false} = true ]]; then
    singularity_bin >/dev/null || { echo "neither singularity nor apptainer found on PATH" >&2; exit 1; }
    [[ -s resources/saige.sif ]] || { echo "resources/saige.sif missing: bash download_resources.sh --saige-image --singularity" >&2; exit 1; }
  else
    command -v docker >/dev/null || { echo "docker not found on PATH" >&2; exit 1; }
    [[ -s resources/saige.tar ]] || { echo "resources/saige.tar missing: bash download_resources.sh --saige-image" >&2; exit 1; }
  fi
}
extra_mounts () {   # one directory per line, from SAIGE_EXTRA_MOUNTS
  local d IFS=:
  for d in ${SAIGE_EXTRA_MOUNTS:-}; do [[ -n ${d} ]] && echo "${d}"; done
}
run_container () {
  local d mounts=() binds=() volumes=()
  while IFS= read -r d; do
    # docker would create a missing directory (as root) rather than complain
    [[ -d ${d} ]] || { echo "SAIGE_EXTRA_MOUNTS: ${d} is not a directory" >&2; return 1; }
    mounts+=("${d}"); binds+=(--bind "${d}:${d}:ro"); volumes+=(-v "${d}:${d}:ro")
  done < <(extra_mounts)
  if [[ ${DRYRUN:-false} = true ]]; then
    echo "DRY RUN -- would run inside the SAIGE image:"; echo "  ${cmd}"
    (( ${#mounts[@]} > 0 )) && echo "  with read-only mounts: ${mounts[*]}"
    return 0
  fi
  local gates=()
  [[ -n ${SAIGE_FIT_GATES:-} ]] && gates=(SAIGE_FIT_GATES="${SAIGE_FIT_GATES}")
  if [[ ${SINGULARITY} = true ]]; then
    env ${gates[@]+"${gates[@]}"} "$(singularity_bin)" exec \
      --home "${WD}" --pwd "${WD}" \
      --bind "${WD}:${WD}" ${binds[@]+"${binds[@]}"} \
      "resources/saige.sif" ${cmd}
  else
    local ref image_id=""
    ref=$(cat resources/saige.image 2>/dev/null || true)
    # load the tar only when the image is not there already (~5 s per call otherwise)
    if [[ -z ${ref} ]] || ! docker image inspect "${ref}" > /dev/null 2>&1; then
      docker load -i "resources/saige.tar" > /dev/null
    fi
    [[ -n ${ref} ]] && image_id=$(docker images --filter=reference="${ref}" --format "{{.ID}}" | head -n 1)
    [[ -n ${image_id} ]] || image_id=$(docker images --format "{{.Repository}}:{{.Tag}} {{.ID}}" | awk '$1 ~ /saige-slim:/ {print $2; exit}')
    [[ -n ${image_id} ]] || { echo "no saige-slim image loaded; run download_resources.sh --saige-image" >&2; return 1; }
    local envs=(-e "HOME=${WD}")
    [[ -n ${SAIGE_FIT_GATES:-} ]] && envs+=(-e "SAIGE_FIT_GATES=${SAIGE_FIT_GATES}")
    docker run --rm "${envs[@]}" -v "${WD}/:${WD}/" ${volumes[@]+"${volumes[@]}"} "${image_id}" ${cmd}
  fi
}
