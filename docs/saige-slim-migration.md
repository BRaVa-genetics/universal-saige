# The move to `saige-slim` and the All of Us parameters (branch `saige-slim-aou-parameters`, 2026-09-28)

What changed, why, and how it was checked -- for whoever picks this repository
up next. The parameter record it implements is
`docs/state/aou-saige-parameters.md` in the `saige-slim` repository
(github.com/astheeggeggs/saige-slim), which names the source of every value.

## Decisions (owner, 2026-09-28)

1. **Image**: the slim SAIGE build, `astheeggeggs/saige-slim:<tag>` from Docker
   Hub (free to pull anywhere; gcr stays for All of Us, whose workers reach only
   Google registries). Consequences carried into the drivers and README:
   - **no VCF reader**: step 2 refuses `--vcf` with the one-line plink2
     conversion; converting per run would charge every chromosome of every
     phenotype for a one-off job. Step 0 still reads a VCF through plink 1.9
     for the GRM.
   - **PLINK 2 (.pgen) is the recommended input**, and a `--pgen` option exists
     beside `--plink` in steps 1 and 2. Measured on the SAIGE fixture: the two
     arms give byte-identical variance ratios and byte-identical FlexRV results.
   - **the fit gates are on**: a binary trait under 100 cases, a covariate level
     under 10 cases or controls, a separated model or a non-converged fit is
     refused, with the gate named. `SAIGE_FIT_GATES=0` in the environment
     downgrades it (the runner passes it through); the README says not to.
2. **Step-2 defaults = the All of Us production choices**: Firth off, fastTest
   off, `--minMAC 4` for single-variant tests (0.5 for group tests), the build's
   missingness defaults (0.15, best_guess), `--relatednessCutoff` passed
   explicitly and required to equal step 1's.
3. **Step 1**: `--tol 0.02` for both trait types (was 1e-5 for quantitative);
   `--useSparseGRMforVarRatio=TRUE` ADDED -- without it the variance-ratio file
   has no `sparse` rows and step 2, given the GRM, refuses the model on the slim
   build (found running the drivers end to end; the All of Us models carry
   those rows).
4. **FlexRV** as a first-class step-2 mode (`--flexRVscore NAME`): one mask,
   one max MAF (0.001), burden statistic, the group file's `score:NAME` line.
   SAIGE reads ONE score line per file, so one file and one run per weight set.
   The weights come from `04_flexrv_groupfile.sh` around
   `flexrv_score_from_alphamissense.py` (copied from `saige-slim`
   `harness/tools/`; standard-library Python; `--selftest`): AlphaMissense first
   (LoF 1.0 by annotation, missense the AlphaMissense pathogenicity, missing
   missense the gene mean), a second consortium weight set through the same door
   as a column of the BRaVa long-form table (`--scoreColumn`) or an
   AlphaMissense-shaped per-variant table (`--am`).

## Other fixes made on the way

- `nproc` is Linux-only; `ncpu()` in `run_container.sh` is portable.
- `check_container_env` was called by step 0 but defined nowhere; it now checks
  the runtime and the image.
- A refused step-1 fit leaves 0-byte `.rda`/`.varianceRatio.txt` (SAIGE creates
  them before the eligibility check; `saige-slim` LEDGER #171); the step-1 driver
  removes them on a non-zero exit so a refused trait leaves no false model.
- `03_estimate_nGlmm.sh` ran the old image by hand, patched R source inside it
  and reinstalled the package; the slim image has `extractNglmm.R` on its PATH.
- `--dryRun` on steps 1 and 2 prints the SAIGE command.

## How it was checked (2026-09-28, local Docker, image `saige-slim:dev` = `1.5.2-dev-3e92d89d`)

On the SAIGE fixture (`extdata/input` of `saige-slim`: 1,000 samples, 98 cases;
the 0.05 sparse GRM; `group_new_snpid_flexrv.txt` with its `score:REVEL` line):
- gate control: step 1 with the gates on REFUSED the fixture (`exit 1`,
  `cases in the fitted set: 98 of 1000`), and left no file behind;
- step 1, `--genotypePlink` and `--genotypePgen`, `SAIGE_FIT_GATES=0`: both fit,
  variance-ratio files byte-identical, `sparse`/`null`/`null_noXadj` rows present;
- step 2: variant (`--pgen`), SKAT-O group (`--annotations lof,missense,lof:missense`)
  with the sidecars, FlexRV (`--flexRVscore REVEL --annotations lof:missense
  --flexRVlofAnno lof --flexRVmaxMAF 0.01`) on `--pgen` and on `--plink`,
  byte-identical (112 transform sets at N = 1,000, LEDGER #142);
- `--dryRun` for both drivers; every refusal path of step 2 (VCF, both genotype
  kinds, several masks or max MAFs under FlexRV, `--testType variant` with
  FlexRV);
- `python3 flexrv_score_from_alphamissense.py --selftest`.

To repeat: `docker save <image> -o resources/saige.tar; echo <image> > resources/saige.image`,
copy the fixture files into a working directory, run the commands above from it.
