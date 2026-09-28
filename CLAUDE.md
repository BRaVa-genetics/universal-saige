# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Bash drivers that run SAIGE (null GLMM fit, single-variant, SAIGE-GENE+ and FlexRV group tests) inside a container, for BRaVa analysts across biobanks. The SAIGE build is `astheeggeggs/saige-slim` (Docker Hub), pinned by `SAIGE_IMAGE`/`SAIGE_VERSION` in `download_resources.sh`. `docs/saige-slim-migration.md` records the current parameter choices, why each was made, and how they were checked.

## Commands

There is no build or linter. Setup and checks:

```bash
bash download_resources.sh --saige-image --plink --plink2 [--singularity]   # image -> resources/saige.{tar,sif}, resources/saige.image
bash download_resources.sh --alphamissense                                   # FlexRV weights source (~600 MB)

bash tests/run_e2e.sh            # the end-to-end test: every driver on a simulated cohort, then PASS/FAIL per check (Docker, 20-40 min)
bash tests/run_e2e.sh --check    # re-score the last run without re-running it
bash tests/run_e2e.sh --resume   # re-run only what did not succeed (e.g. after Docker died mid-run)
python3 flexrv_score_from_alphamissense.py --selftest   # controls for the FlexRV score builder alone; no downloads
bash 01_step1_fitNULLGLMM.sh ... --dryRun               # print the SAIGE command instead of running it (same flag on step 2)
```

`tests/simulate_cohort.py` simulates a cohort where the right answer is known: a family pedigree the GRM should recover, planted gene, FlexRV and common-variant effects, null and permuted traits, a binary trait below the 100-case gate, and untidy sample-ID overlaps whose intersection is the expected N. `tests/check_results.py` compares every step's output with `tests/work/in/truth.json`, and its docstring states each expected effect size. When you change a driver, add a check there. Output goes to `tests/work/` (gitignored), and each run's log is `tests/work/logs/<name>.log`. The Singularity path is not covered.

Performance notes, both measured:
- Step 2 is single-threaded, because the image sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `RCPP_PARALLEL_NUM_THREADS` to 1. Don't pass the host's values through: HPC schedulers often set `OMP_NUM_THREADS` and would override that.
- On Docker Desktop for macOS, bind-mount I/O leaves `com.docker.backend` burning CPU after the containers exit. A long test run can slow down about 10× this way. Restarting Docker fixes it.

## Architecture

**Every driver must be run from the repo root.** Each one does `source ./run_container.sh`, sets `WD=$(pwd)` and `HOME=$WD`, parses flags into globals, builds a SAIGE command string in `$cmd`, and calls `run_container`. That function mounts `$WD` at the same path inside Docker or Singularity and runs `$cmd` unquoted. So:
- every file argument must be a relative path inside the working directory. Drivers prefix it with `${HOME}/`, and symlinks pointing outside the directory break. The one exception is step 3, which also takes absolute paths under directories that `run_container` binds read-only: `/mnt/project` (UKB RAP) whenever it exists, and anything in `SAIGE_EXTRA_MOUNTS`;
- `$cmd` is word-split, so paths and column names cannot contain spaces;
- the image is found through `resources/saige.image` (Docker: `docker load` of `resources/saige.tar` on each call) or `resources/saige.sif`.

Pipeline and hand-offs (outputs of one step are flags of the next):

| Step | Where it runs | Produces |
|---|---|---|
| `00_step0_VR_and_GRM.sh` | host plink 1.9 (merge, LD-prune, seeded MAC-stratified marker sample; scratch in a per-run `mktemp -d` under `$TMPDIR`), then `createSparseGRM.R` in the container | `<out>.plink_for_var_ratio.{bed,bim,fam}`, `<out>_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx` (+ `.sampleIDs.txt`) |
| `01_step1_fitNULLGLMM.sh` | `step1_fitNULLGLMM.R` | `<out>.rda`, `<out>.varianceRatio.txt` |
| `02_step2_SPAtests_variant_and_gene.sh` | `step2_SPAtests.R`, one chromosome | `<out>.txt` (+ `.singleAssoc.txt`, `.markerList.txt`, sidecars for group tests) |
| `03_estimate_nGlmm.sh` | `extractNglmm.R`, once per phenotype | `neff.csv` (`--outputFile`), where Nglmm = 1ᵀK⁻¹1 on the GRM (× 4Pn(1−Pn) for binary traits) |
| `04_flexrv_groupfile.sh` | host `python3` wrapper around `flexrv_score_from_alphamissense.py` (stdlib only) | group file with a `score:NAME` line per gene |

Step numbering is not execution order: step 4 builds the FlexRV group file that step 2 `--flexRVscore NAME` then reads.

## Invariants that are easy to break

- **The parameters are the All of Us production choices.** Their authoritative record is `docs/state/aou-saige-parameters.md` in `saige-slim`. Don't change them without a reason recorded there: step 1 `--tol 0.02` for both trait types, IRNT for quantitative traits; step 2 Firth off, fastTest off, `--minMAC 4` for variant tests and `0.5` for group tests, max MAF `0.0001,0.001,0.01`, `--LOCO=FALSE`.
- **Step 1 `--useSparseGRMforVarRatio=TRUE` is required.** Without it the variance-ratio file has no `sparse` rows, and step 2 refuses the model when it is given the GRM.
- `--relatednessCutoff` (0.05) must be the same in steps 1 and 2. SAIGE does not check this.
- **PLINK 1 or PLINK 2 input only, in every step.** The image has no VCF reader, and steps 0 and 2 refuse a VCF on purpose, with the one-line plink2 conversion in the message. Step 0 takes `--geneticDataFormat {plink,pgen}` and converts `.pgen` to `.bed` in its scratch directory, because its merge, pruning and counts are plink 1.9. Steps 1 and 2 take `--pgen` (recommended) or `--plink`/`--genotypePlink`, exactly one of the two.
- **Fit gates are on by default.** A binary trait with fewer than 100 cases, a categorical covariate level with fewer than 10 cases or controls, a separated model, or a fit that did not converge is refused. `SAIGE_FIT_GATES=0` turns refusals into warnings and `run_container` passes it through. After a refused fit, step 1 deletes the 0-byte `.rda`/`.varianceRatio.txt` it leaves behind. Keep that cleanup.
- `--chr` must match the chromosome spelling in the `.pvar`/`.bim` exactly (`20` vs `chr20`).
- Sample IDs are matched on IID throughout. FID is ignored: plink2 writes FID `0` when a `.psam` has only IIDs.
- **The drivers must run on macOS as well as Linux.** macOS ships bash 3.2, where an empty array is "unbound" under `set -u`; use `${arr[@]+"${arr[@]}"}`. Avoid GNU-only tools such as `shuf` and `nproc`, and use `ncpu()` from `run_container.sh` for the core count.
- Never edit `tests/run_e2e.sh`, or a driver it is running, while a test run is in progress. Bash reads scripts from disk as it executes them.
- **FlexRV:** a group file carries one `score:NAME` line per gene, so each weight set needs its own file and its own step-2 run. The score line goes directly after the `anno` line and before any `weight` line. It needs one value in [0, 1] per variant on the `var` line, including variants the masks later drop. FlexRV mode forces one mask, one max MAF and `--r.corr=1`. The pooled p-value per gene is the `Group == Cauchy` row.
- `flexrv_score_from_alphamissense.py` was copied from `saige-slim` `harness/tools/`. Paths in its docstring (`harness/tools/…`, `docs/design/flexrv-design.md`) refer to that repo.

## Keeping docs in sync

Driver flags are documented in four places: each script's `-h` text, `README.md`, `templates/step_{0,1,2}_template.sh`, and `walkthrough/README.md`. A flag change has to be made in all four. The templates are fill-in-the-blanks for analysts, not runnable scripts: `[free text]` marks a value to replace and `{a,b}` a choice.

Code style: argument parsing is a `while/case` block with a `-h` usage string. Comments explain *why* a value was chosen and cite where it came from (AoU record, LEDGER number, paper). Match that when adding parameters.
