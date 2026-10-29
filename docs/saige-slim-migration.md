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
     phenotype for a one-off job. (Step 0 then still read a VCF through plink
     1.9; since the second pass below it refuses one too.)
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
   Since 2026-10-21 `--amScale percentile` gives the FlexRV paper's score (a missense
   variant's percentile among all possible missense variants of its gene; README,
   step 4), which the BRaVa release ships as `flexrv_AMPCT` and `flexrv_SCALLION`.

## Other fixes made on the way

- `nproc` is Linux-only; `ncpu()` in `run_container.sh` is portable.
- `check_container_env` was called by step 0 but defined nowhere; it now checks
  the runtime and the image.
- A refused step-1 fit leaves 0-byte `.rda`/`.varianceRatio.txt` (SAIGE creates
  them before the eligibility check; `saige-slim` LEDGER #171); the step-1 driver
  removes them on a non-zero exit so a refused trait leaves no false model.
  saige-slim fixed #171 on 2026-10-04 and the pinned `928f95ad` includes it;
  the driver's cleanup stays, a no-op now, against an older image.
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

## Second pass (2026-09-28): an end-to-end test, and what it found

`tests/run_e2e.sh` runs every driver on a simulated cohort where the right
answer is known: a pedigree for the GRM to recover, planted gene, FlexRV and
common-variant effects, null and permuted traits, a trait below the case
gate, and sample-ID overlaps whose intersection is the expected N. It is
scored by `tests/check_results.py` (see `CLAUDE.md`). What it and the review
around it changed:

- **Step 0 takes PLINK 1 or PLINK 2 only** (`--geneticDataFormat {plink,pgen}`),
  like steps 1 and 2. A VCF is refused with the plink2 conversion. `.pgen`
  is converted to `.bed` in the scratch directory, because the merge, pruning
  and counts are plink 1.9.
- **Step 0 `--sampleIDs` matches on IID.** `--keep "ID ID"` silently dropped
  every sample whose FID differed from its IID, which covers any PLINK 2 input
  (plink2 writes FID 0).
- **Step 0 MAC bins:** the minor allele count was computed as C2 - C1 from
  `.frq.counts`, whose C2 is the other allele count, not the total. Near-50%
  variants were binned as MAC 10-20 (1% of that bin in the test).
- Step 0: `shuf` (absent on macOS) replaced by a seeded sampler; scratch in a
  per-run `mktemp -d` instead of shared `/tmp` names; `--outputPrefix` required.
- Step 1 `--sex` read unset variables (so it read stdin) and could not stop the
  run from inside `| while`. It then required one value of `sex` among samples
  with a phenotype (replaced in the third pass below).
- Step 3 mounted the working directory at the host `$HOME` and took `$?` from
  `tee`, so it only ran from `$HOME` and hid failures. It now uses
  `run_container` (Singularity, fit-gate pass-through), takes `--outputFile`,
  and exits non-zero for a phenotype without an Nglmm. Inputs outside the working
  directory go through read-only binds in `run_container`, listed in
  `SAIGE_EXTRA_MOUNTS`; step 3 takes absolute paths under them. Nglmm
  reconciles exactly with 1'K^-1 1 on the GRM.
- Step 3 never passed `--relatednessCutoff`, so `extractNglmm.R` thinned the GRM
  at its own default, 0.125, while steps 1 and 2 used 0.05: Nglmm described a
  sparser GRM than the model was fitted on. It now passes the same cutoff as
  steps 0-2 (default 0.05). The e2e check had mirrored the 0.125 and so could
  not see this; it now recomputes 1'K^-1 1 at 0.05. Nglmm changes slightly
  against earlier BRaVa submissions, which all used 0.125.
- `04_flexrv_groupfile.sh` failed under macOS bash 3.2 (an empty array under `set -u`).
- Templates: `--t` -> `--traitType`; the GRM file is `<out>_relatednessCutoff_...`.

## Third pass (2026-10-06/07): AoU parity checked against its job logs, and guards

The AoU step-1 wrapper (`saige_step1.sh`) and three of its v4 job logs, and the
step-2 wrapper, were read flag by flag against the drivers; the record in
saige-slim was completed from them. The end-to-end test grew from 111 to 200
checks; it also passes on a SLURM cluster under Apptainer.

- **Runtime**: Apptainer is used when `singularity` is not on PATH; the `.sif`
  is pulled without a cache and re-pulled when the pin changes. The UKB RAP
  `/mnt/project` bind is gone (`SAIGE_EXTRA_MOUNTS` remains). Image
  `1.5.2-dev-928f95ad` (LEDGER #171, #172 and #176 fixed, and a
  `.stretchGate.txt` with fixed columns; results byte-identical to
  `152ffd8c`, which came first); plink2 alpha 7.11 and
  plink 1.9.0 stable.
- **Step 1**: IRNT for quantitative traits, as AoU (now in the record, from the
  logs). `--relatednessCutoff` is a flag in steps 0, 1 and 3, as in step 2: AoU
  fitted amr at 0.125. Phenotype and covariate names must be plain R names
  (SAIGE pastes them into a formula, LEDGER #68; the old check accepted
  everything). An empty `--sampleIDs` is refused (it swallowed the next flag).
- **Dense GRM**: steps 0 and 1 print relatives per sample from the `.mtx`
  header; step 1 refuses above 100 unless `--forceDenseGRM` (AoU amr at 0.05:
  ~644, fits that never finished; at 0.125 ~3.9; other cohorts ~0.6).
- **`--sex`**: passes `--FemaleOnly/--MaleOnly --sexCol sex` as AoU did, with
  BRaVa's coding (0 = female, 1 = male; SAIGE's default and the AoU file are the
  opposite) and `--femaleCode/--maleCode`. Codes that select nobody or put every
  case in the dropped sex are refused; suspicious splits warned; genetic sex in
  the `.fam`/`.psam` checks the codes (over half discordant refused as flipped).
  Verified: `--sex F` gives byte-identical variance ratios to males set NA.
- **Step 2**: refused before SAIGE runs: mask labels absent from the group file,
  a FlexRV score line missing from any region, a `--chr` spelled unlike the
  `.pvar`/`.bim`. The four new sidecars of `152ffd8c` are documented.
- **Image `1.5.2-dev-b5456dee`** (2026-10-16, replacing `928f95ad`): LEDGER #181,
  a non-finite saddlepoint in a burden cell's variance rescale killed the chromosome
  (AoU HTN); the cell now takes the normal-approximation p, disclosed in the log and in
  `.spaRegionNA.txt`. Also LEDGER #179 (faster burden / single-in-group runs, output
  byte-identical) and SAIGE's `--is_output_log10P`, which step 2 passes on group tests by default
  (`--log10P`, default TRUE, from 2026-10-16): log10 p-value columns after `SE_Burden`.
- **Image `1.5.2-dev-49178058`** (2026-10-16, replacing `b5456dee`): SAIGE's `--is_output_nEff`, passed
  by default (`--nEff`): `nEff_Burden`, Kish's effective carrier count of each cell's burden vector, after
  the log10 columns; added after FlexRV height hits in All of Us rested on fewer than five effective
  carriers. Step 2 also writes per-variant results only on SAIGE-GENE+ runs by default (`--singleInGroup`,
  FALSE for FlexRV), as the All of Us runs did from 2026-10-08.
- **Image `1.5.2-dev-41d78ae3`** (2026-10-29, replacing `49178058`): LEDGER #183. On the sparse-GRM band a
  thinned GRM's Sigma can be indefinite, and a converged solve then returns a negative or tiny variance (an
  All of Us amr random-score cell printed p 1e-73 one step before the p 1 / BETA NA collapse). Every step-2 run
  now compares each such variance against the no-GRM variance of the same vector and, when the solve is not
  positive or below half of it, reports the no-GRM variance instead, naming the test in `.sigmaGuard.txt`
  (image `66b80710`; measured on the rebuilt All of Us amr fleet: 0 flagged tests on 672 real-arm jobs). The
  image also carries `aou_sigma_pd_check.R`, and the new `05_sigma_check.sh` runs it ONCE per GRM and cutoff,
  after step 0 and before any step 1 (`--grmOnly`: if the thinned GRM is positive semidefinite, every model's Sigma
  is positive definite; exit 2 and a banner otherwise: rebuild the GRM, fit nothing on it). `--modelFile` checks one
  model's own Sigma, a diagnosis.
  The All of Us amr GRM failed it and was rebuilt from KING-robust kinship (`king_to_mtx.py`); its step 1 and
  step 2 were rerun. The step-2 template gains run 5: the `flexrv_RAND` twin on every phenotype, FlexRV's
  file-level null, read per trait.
- **Step 3**: Nglmm at the same cutoff as the fit (default 0.05; it was
  `extractNglmm.R`'s own 0.125, which every earlier BRaVa submission used).
- **Docs**: templates use 20 PCs, as AoU did; walkthrough and templates checked
  against the drivers (an example that merged the exome into the GRM, outdated
  mask labels, a lambda_GC computed from p-values); every driver's `-h` now exits
  after printing.
