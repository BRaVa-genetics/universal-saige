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
# DRYRUN=true prints the command instead of running it.
ncpu () {   # portable core count (nproc is Linux-only)
  nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1
}
check_container_env () {   # $1 = true for singularity; stops early when the runtime or the image is missing
  if [[ ${1:-false} = true ]]; then
    command -v singularity >/dev/null || { echo "singularity not found on PATH" >&2; exit 1; }
    [[ -s resources/saige.sif ]] || { echo "resources/saige.sif missing: bash download_resources.sh --saige-image --singularity" >&2; exit 1; }
  else
    command -v docker >/dev/null || { echo "docker not found on PATH" >&2; exit 1; }
    [[ -s resources/saige.tar ]] || { echo "resources/saige.tar missing: bash download_resources.sh --saige-image" >&2; exit 1; }
  fi
}
run_container () {
  if [[ ${DRYRUN:-false} = true ]]; then
    echo "DRY RUN -- would run inside the SAIGE image:"; echo "  ${cmd}"; return 0
  fi
  local gates=()
  [[ -n ${SAIGE_FIT_GATES:-} ]] && gates=(SAIGE_FIT_GATES="${SAIGE_FIT_GATES}")
  if [[ ${SINGULARITY} = true ]]; then
    env "${gates[@]}" singularity exec \
      --home "${WD}" --pwd "${WD}" \
      --bind "${WD}:${WD}" \
      "resources/saige.sif" ${cmd}
  else
    docker load -i "resources/saige.tar" > /dev/null
    local ref image_id=""
    ref=$(cat resources/saige.image 2>/dev/null || true)
    [[ -n ${ref} ]] && image_id=$(docker images --filter=reference="${ref}" --format "{{.ID}}" | head -n 1)
    [[ -n ${image_id} ]] || image_id=$(docker images --format "{{.Repository}}:{{.Tag}} {{.ID}}" | awk '$1 ~ /saige-slim:/ {print $2; exit}')
    [[ -n ${image_id} ]] || { echo "no saige-slim image loaded; run download_resources.sh --saige-image" >&2; return 1; }
    local envs=(-e "HOME=${WD}")
    [[ -n ${SAIGE_FIT_GATES:-} ]] && envs+=(-e "SAIGE_FIT_GATES=${SAIGE_FIT_GATES}")
    docker run --rm "${envs[@]}" -v "${WD}/:${WD}/" "${image_id}" ${cmd}
  fi
}
