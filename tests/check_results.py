#!/usr/bin/env python3
"""Score a tests/run_e2e.sh run against the simulation's truth.json.

    python3 tests/check_results.py [tests/work]

Prints one line per check (PASS/FAIL, what was observed, what was expected)
and exits non-zero if any check fails. Expectations, written down before the
first run and derived from the simulation (N ~ 1,800 analysed):
  * positive controls: GENE_BURDEN pLoF burden z ~ 0.8 * sqrt(~60 carriers) ~ 6
    (p ~ 1e-9); GENE_FLEX causal-only z ~ 1.5 * sqrt(~30 copies) ~ 8 before
    the relatedness correction, unweighted ~ 1/sqrt(2) of that (half the
    copies are null), so score weighting should win by >= 100x in p and must
    not once the scores are shuffled; the common variant
    z ~ 0.25 * sqrt(2 * 0.3 * 0.7 * N) ~ 7. Pass bar p < 1e-6 (quantitative),
    with the sign of the effect checked;
  * negative controls: no null gene or variant below 1e-6 (about 1% chance
    over ~15k tests), lambda_GC of MAC >= 20 variants within 0.85-1.15, and
    GENE_BURDEN p > 1e-3 when the phenotype is permuted or has no effect;
  * step 0 keeps the --sampleIDs samples by IID: the array .fam has FID =
    family and the PLINK 2 arm's .psam has no FID, and both must give 1,900;
  * the sparse GRM recovers the simulated pedigree: first-degree pairs near
    0.5 and no unrelated pair above 0.12 (noise SD ~ 1/sqrt(5000) = 0.014);
  * the same analysis computed two ways (PLINK 1 vs PLINK 2 input) is
    byte-identical, as it was on the saige-slim fixture; step 3's Nglmm
    equals 1'K^-1 1 recomputed here from the GRM file.
"""

import json
import math
import os
import statistics
import sys

W = sys.argv[1] if len(sys.argv) > 1 else "tests/work"
IN, OUT, LOGS = (os.path.join(W, d) for d in ("in", "out", "logs"))
T = json.load(open(os.path.join(IN, "truth.json")))
RESULTS = []
status = {}
for line in open(os.path.join(W, "status.tsv")):
    name, expect, rc, secs = line.rstrip("\n").split("\t")
    status[name] = (expect, int(rc), int(secs))


SECTION = []


def section(*runs):
    """The runs whose output the checks that follow read."""
    SECTION[:] = runs


SKIPPED = []


def skip(name, why):
    """A check that cannot be informative on this machine: reported, not counted."""
    SKIPPED.append((name, why))


def check(name, ok, observed, expected, runs=None):
    """If a run the check reads did not end as expected, its output may be
    partial (a killed run leaves half a file), so the check fails whatever it says."""
    bad = [r for r in (SECTION if runs is None else runs) if not completed(r)]
    if bad:
        ok, observed = False, "run %s did not complete as expected" % ", ".join(bad)
    RESULTS.append((bool(ok), name, str(observed), str(expected)))


def completed(name):
    if name not in status:
        return False
    expect, rc, _ = status[name]
    return (rc == 0) == (expect == "ok")


def out(p):
    return os.path.join(OUT, p)


def log(name):
    try:
        return open(os.path.join(LOGS, name + ".log"), errors="replace").read()
    except OSError:
        return ""


def table(path):
    """Rows of a whitespace-delimited file with a header, as dicts."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path) as fh:
        header = fh.readline().split()
        return [dict(zip(header, line.split())) for line in fh if line.strip()]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def same_bytes(a, b):
    try:
        return open(a, "rb").read() == open(b, "rb").read()
    except OSError:
        return False


Z = statistics.NormalDist()
CHI2_MEDIAN = 0.454936423119572


def lambda_gc(ps):
    chis = [Z.inv_cdf(p / 2) ** 2 for p in ps if 0 < p <= 1]
    return statistics.median(chis) / CHI2_MEDIAN if chis else float("nan")


def fmt_p(p):
    return "%.2g" % p


# ------------------------------------------------------------ exit codes
section()
for name, (expect, rc, secs) in status.items():
    ok = (int(rc) == 0) == (expect == "ok")
    check("exit: " + name, ok, "exit %s (%ss)" % (rc, secs), "exit 0" if expect == "ok" else "non-zero")
for name in ("s2_B_pos_variant", "s2_B_pos_group", "s3_nglmm", "s2_Q_null_group", "s2_Q_perm_group",
             "s3_nglmm_extmount"):
    if name not in status:
        check("exit: " + name, False, "never ran", "exit 0")

# ------------------------------------------------------------ refusal messages
for name, text in [
    ("s2_refuse_vcf", "--vcf is not accepted"),
    ("s2_refuse_both_geno", "pass ONE of --pgen and --plink"),
    ("s2_refuse_flex_masks", "FlexRV tests ONE annotation mask"),
    ("s2_refuse_flex_mafs", "FlexRV tests ONE max MAF"),
    ("s2_refuse_flex_variant", "contradicts"),
    ("s2_dryrun", "DRY RUN"),
    ("s2_refuse_anno_labels", "on no variant in the group file: damaging_missense other_missense"),
    ("s2_refuse_lof_label", "on no variant in the group file: LoF"),
    ("s2_refuse_no_score", "has 0 for 300 region(s)"),
    ("s2_refuse_partial_score", "has 299 for 300 region(s)"),
    ("s2_refuse_chr_name", "--chr chr7 is not in column 1 of"),
    ("s1_refuse_both_geno", "Pass ONE of --genotypePlink and --genotypePgen"),
    ("s1_sex_mf_ok", "--FemaleOnly=TRUE --sexCol=sex --FemaleCode=F"),
    ("s1_sex_mf_wrong", "no sample with sex == M has a non-missing Q_f"),
    ("s1_sex_mixed", "More samples with a Q_null are of the sex being dropped (1013) than kept (897)"),
    ("s1_sex_code_absent", "is coded 0, which is not a value of the 'sex' column"),
    ("s1_sex_binary_flipped", "cases of B_f have sex == 0, the sex being dropped: the codes look flipped"),
    ("s1_sex_binary_both", "cases of B_pos have sex == 1, the sex being dropped"),
    ("s1_sex_genetic_flipped", "100.0% of samples disagree with their genetic sex"),
    ("s1_sex_genetic_unknown", "has no usable genetic sex (0 male, 0 female"),
    ("s1_sex_genetic_partial", "94 of 1873 samples (5.0%) have a 'sex' that disagrees with their genetic sex"),
    ("s1_sex_genetic_pgen", "plink_for_var_ratio.psam) agrees with the phenotype file's for 1873 of 1873 samples"),
    ("s1_Q_female", "plink_for_var_ratio.fam) agrees with the phenotype file's for 1873 of 1873 samples"),
    ("s1_sex_numeric_ok", "--FemaleOnly=TRUE --sexCol=sex --FemaleCode=0"),
    ("s1_sex_flipped_code", "no sample with sex == 1 has a non-missing Q_female"),
    ("s1_refuse_sex_covar", "'sex' is also a covariate"),
    ("s0_refuse_vcf", "geneticDataFormat must be in {plink,pgen}"),
    ("s0_refuse_no_out", "--outputPrefix is required"),
    ("s0_refuse_bad_relcut", "--relatednessCutoff must be a number between 0 and 1"),
    ("s1_refuse_bad_relcut", "--relatednessCutoff must be a number between 0 and 1"),
    ("s1_relcut_dryrun", "--relatednessCutoff 0.125"),
    ("s3_nglmm_nomount", "does not exist"),
    ("s3_refuse_bad_mount", "is not a directory"),
    ("grm_dense_005", "THE SPARSE GRM IS DENSE: 399.0 relatives per sample"),
    ("grm_relcut_mismatch", "built at --relatednessCutoff 0.05 and this step uses 0.125"),
    ("s1_refuse_dense", "REFUSED: the sparse GRM is too dense to fit"),
    ("s1_dense_override", "--forceDenseGRM: fitting on the dense GRM anyway"),
    ("s0_plink", "1900 samples, 2267 related pairs, 2.4 relatives per sample (built at --relatednessCutoff 0.05)"),
    ("s0_relcut_0125", "1900 samples, 2258 related pairs, 2.4 relatives per sample (built at --relatednessCutoff 0.125)"),
    ("s1_Q_pos", "relatives per sample"),
]:
    if os.path.exists(out(name + ".skipped")):
        skip("message: " + name, "the container sees %s without a mount (apptainer binds /tmp and the "
             "site's bind paths), so the run cannot fail; the SAIGE_EXTRA_MOUNTS check below then shows only "
             "that binding a visible directory is harmless" % open(out(name + ".skipped")).read().strip())
        continue
    check("message: " + name, text in log(name), "found" if text in log(name) else "absent", repr(text))

# the dense-GRM banner: only where the GRM is dense at the cutoff in use
for name in ("grm_relcut_mismatch", "s0_plink", "s0_relcut_0125", "s1_Q_pos"):
    check("no dense-GRM warning: " + name, "THE SPARSE GRM IS DENSE" not in log(name),
          "absent" if "THE SPARSE GRM IS DENSE" not in log(name) else "PRESENT", "absent", runs=[name])
for name in ("s0_plink", "s0_relcut_0125", "s1_Q_pos"):
    check("no cutoff-mismatch warning: " + name, "this step uses" not in log(name),
          "absent" if "this step uses" not in log(name) else "PRESENT", "absent", runs=[name])

# step 0 ends on the GRM's density, where it cannot be missed
for name in ("s0_plink", "s0_pgen", "s0_relcut_0125"):
    last = ([l for l in log(name).splitlines() if l.strip()] or [""])[-1]
    check("step0 log ends on the GRM density: " + name, "relatives per sample" in last, last[-60:], "the density line",
          runs=[name])

# ------------------------------------------------------------ step 0
kept_n = T["n_kept"]
section("s0_plink")
fam = out("step0.plink_for_var_ratio.fam")
n_fam = sum(1 for _ in open(fam)) if os.path.exists(fam) else 0
check("step0: VR plink samples = include list x genotyped", n_fam == kept_n, n_fam, kept_n)
freq = table(out("vr_freq.frq.counts"))
if freq:
    macs = [min(int(r["C1"]), int(r["C2"])) for r in freq]
    lo = sum(1 for m in macs if 10 <= m < 20)
    hi = sum(1 for m in macs if m >= 20)
    check("step0: VR markers MAC 10-20 / MAC >= 20", (lo, hi) == (2000, 2000), (lo, hi), (2000, 2000))
else:
    check("step0: VR markers MAC 10-20 / MAC >= 20", False, "no frq file", (2000, 2000))


def grm_check(prefix, label, run, cutoff="0.05"):
    """The GRM step 0 wrote at --relatednessCutoff `cutoff`; returns the number of unrelated pairs it kept."""
    section(run)
    ids_path = out(prefix + "_relatednessCutoff_%s_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" % cutoff)
    mtx_path = out(prefix + "_relatednessCutoff_%s_5000_randomMarkersUsed.sparseGRM.mtx" % cutoff)
    if not (os.path.exists(ids_path) and os.path.exists(mtx_path)):
        check(label + ": GRM files exist", False, "missing", "present")
        return None
    ids = [l.strip() for l in open(ids_path) if l.strip()]
    check(label + ": GRM samples = include list x genotyped", len(ids) == kept_n and len(set(ids)) == kept_n,
          len(ids), kept_n)
    entries = {}
    with open(mtx_path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            dims = line.split()
            break
        for line in fh:
            i, j, v = line.split()
            i, j = int(i) - 1, int(j) - 1
            entries[(min(i, j), max(i, j))] = float(v)
    idx = {s: k for k, s in enumerate(ids)}
    diag = [v for (i, j), v in entries.items() if i == j]
    check(label + ": GRM diagonal mean ~ 1", 0.9 < statistics.mean(diag) < 1.1,
          "%.3f (n=%d)" % (statistics.mean(diag), len(diag)), "0.9-1.1")
    ped = set()
    for a, b in T["first_degree_pairs"]:
        i, j = idx[a], idx[b]
        ped.add((min(i, j), max(i, j)))
    found = [entries.get(p) for p in ped]
    good = sum(1 for v in found if v is not None and 0.35 <= v <= 0.65)
    check(label + ": first-degree pairs recovered at 0.35-0.65", good >= 0.95 * len(ped),
          "%d of %d (median %.3f)" % (good, len(ped), statistics.median(v for v in found if v is not None)),
          ">= 95%")
    spurious = [v for (i, j), v in entries.items() if i != j and (i, j) not in ped]
    check(label + ": no unrelated pair above 0.12", not spurious or max(spurious) < 0.12,
          "%d unrelated pairs kept, max %.3f" % (len(spurious), max(spurious) if spurious else 0), "< 0.12")
    low = [v for (i, j), v in entries.items() if i != j and v < float(cutoff)]
    check(label + ": no pair below the cutoff %s" % cutoff, not low,
          "%d below, min %.3f" % (len(low), min(low)) if low else "none", "none")
    return len(spurious)


n_spur_005 = grm_check("step0", "step0 plink", "s0_plink")
grm_check("step0_pgen", "step0 pgen (IID-only .psam)", "s0_pgen")
# --relatednessCutoff 0.125 (All of Us amr): the pedigree survives (checked above, at
# 0.35-0.65), and the unrelated pairs the 0.05 GRM keeps (max ~0.06) are dropped
n_spur_0125 = grm_check("step0_rc0125", "step0 at --relatednessCutoff 0.125", "s0_relcut_0125", "0.125")
section("s0_plink", "s0_relcut_0125")
check("step0: 0.125 drops the unrelated pairs 0.05 keeps", bool(n_spur_005) and n_spur_0125 == 0,
      "%s at 0.05, %s at 0.125" % (n_spur_005, n_spur_0125), "> 0 at 0.05, 0 at 0.125")

# ------------------------------------------------------------ step 1
for trait in ("Q_pos", "Q_pos_pgen", "Q_null", "Q_perm", "B_pos", "B_rare", "Q_female", "Q_null_female", "Q_null_malesNA"):
    section("s1_B_rare_ungated" if trait == "B_rare" else "s1_" + trait)
    vr = out(trait + ".varianceRatio.txt")
    rda = out(trait + ".rda")
    has = os.path.exists(rda) and os.path.getsize(rda) > 0 and os.path.exists(vr) and os.path.getsize(vr) > 0
    sparse = has and "sparse" in open(vr).read()
    check("step1 %s: model + variance ratio with 'sparse' rows" % trait, has and sparse,
          "ok" if has and sparse else ("no sparse rows" if has else "missing/empty"), "present")
    # quantitative traits are IRNT'd on import (a BRaVa requirement), binary ones never;
    # SAIGE prints this line (its spelling) only when it applies the transform
    irnt = "Perform the inverse nomalization" in log("s1_B_rare_ungated" if trait == "B_rare" else "s1_" + trait)
    quant = trait.startswith("Q_")
    check("step1 %s: IRNT %s" % (trait, "applied" if quant else "not applied"), irnt == quant,
          "applied" if irnt else "not applied", "applied" if quant else "not applied")
section("s1_Q_female", "s1_Q_null_female")
section("s1_Q_null_female", "s1_Q_null_malesNA")
check("step1 --sex F: SAIGE dropping the males == males NA in the file (VR bytes)",
      same_bytes(out("Q_null_female.varianceRatio.txt"), out("Q_null_malesNA.varianceRatio.txt")), "", "identical")
n_used = lambda name: next((l.split()[0] for l in log(name).splitlines() if "samples will be used for analysis" in l), None)
check("step1 --sex F: same samples as males NA", n_used("s1_Q_null_female") == n_used("s1_Q_null_malesNA") is not None,
      "%s vs %s" % (n_used("s1_Q_null_female"), n_used("s1_Q_null_malesNA")), "equal")
for name in ("s1_sex_numeric_ok", "s1_sex_binary_ok", "s1_Q_female", "s1_sex_genetic_pgen"):
    check("no sex-coding warning: " + name, "CHECK THE SEX CODING" not in log(name),
          "absent" if "CHECK THE SEX CODING" not in log(name) else "PRESENT", "absent", runs=[name])
section("s1_Q_female", "s1_Q_null_female")
for trait in ("Q_female", "Q_null_female"):
    left = [f for f in os.listdir(OUT) if f.startswith(trait + "_FemaleOnly")]
    check("step1 --sex F: SAIGE's _FemaleOnly files renamed: " + trait, not left, left or "none", "none")
section("s1_Q_pos", "s1_Q_pos_pgen")
check("step1: VR from --genotypePlink == --genotypePgen (bytes)",
      same_bytes(out("Q_pos.varianceRatio.txt"), out("Q_pos_pgen.varianceRatio.txt")), "", "identical")
section("s1_B_rare_gated")
left = [f for f in os.listdir(OUT) if f.startswith("B_rare_gated")]
check("step1 gate: refused fit leaves no files", not left, left or "none", "none")
gate_log = log("s1_B_rare_gated")
want = "has %d cases in the fitted set, below --minCaseCount=" % T["cases"]["B_rare"]
check("step1 gate: refused for the case count", want in gate_log,
      next((l.strip() for l in gate_log.splitlines() if "below --minCaseCount" in l), "no gate message")[:90],
      repr(want))
section("s2_B_rare_gated")
s2_gate = log("s2_B_rare_gated")
want = "was fitted on %d cases, below --minCaseCount=" % T["cases"]["B_rare"]
check("step2 gate: B_rare model (fitted with gates off) refused on load", want in s2_gate,
      next((l.strip() for l in s2_gate.splitlines() if "below --minCaseCount" in l), "no gate message")[:90],
      repr(want))

# ------------------------------------------------------------ step 4
section("s4_flexrv")
fg = out("group.chr7.flexrv_AM.txt")
if os.path.exists(fg):
    lines = [l.split() for l in open(fg) if l.strip()]
    genes = {}
    order_ok = True
    for k, parts in enumerate(lines):
        if parts[1].startswith("score"):
            order_ok &= lines[k - 1][1] == "anno" and lines[k - 1][0] == parts[0]
        genes.setdefault(parts[0], {})[parts[1]] = parts[2:]
    n_score = sum(1 for g in genes.values() if "score:AM" in g)
    check("step4: every gene has a score:AM line after anno", n_score == T["n_genes"] and order_ok,
          "%d genes, order %s" % (n_score, "ok" if order_ok else "WRONG"), T["n_genes"])
    lens_ok = all(len(g["score:AM"]) == len(g["var"]) for g in genes.values() if "score:AM" in g)
    vals = [fnum(v) for g in genes.values() for v in g.get("score:AM", [])]
    check("step4: one score per variant, all in [0, 1]", lens_ok and all(0 <= v <= 1 for v in vals),
          "lengths %s, range %.3g-%.3g" % ("ok" if lens_ok else "WRONG", min(vals), max(vals)), "")

    def scores(gene, ids):
        g = genes[gene]
        m = dict(zip(g["var"], g["score:AM"]))
        return [fnum(m[i]) for i in ids if i in m]

    b = scores("GENE_BURDEN", T["variants"]["GENE_BURDEN_pLoF"])
    fc = scores("GENE_FLEX", T["variants"]["GENE_FLEX_causal"])
    fn = scores("GENE_FLEX", T["variants"]["GENE_FLEX_null"])
    check("step4: pLoF = 1; planted missense > 0.85; null missense < 0.3",
          b and all(v == 1 for v in b) and fc and min(fc) > 0.85 and fn and max(fn) < 0.3,
          "pLoF %s; causal min %.3f; null max %.3f" % (sorted(set(b)), min(fc), max(fn)), "")
else:
    check("step4: FlexRV group file", False, "missing", "present")

# ------------------------------------------------------------ step 2
N = T["expected_n"]
C = T["common"]


def variant_checks(trait, positive):
    section("s2_%s_variant" % trait)
    rows = table(out(trait + ".variant.txt"))
    if not rows:
        check("step2 %s variant: results" % trait, False, "missing", "present")
        return
    if "N" in rows[0]:
        ns = {r["N"] for r in rows}
        check("step2 %s variant: N = expected join" % trait, ns == {str(N[trait])}, sorted(ns)[:3], N[trait])
    else:
        nc = {r["N_case"] for r in rows}
        nt = {str(int(r["N_case"]) + int(r["N_ctrl"])) for r in rows}
        check("step2 %s variant: N_case, N = expected join" % trait,
              nc == {str(T["cases"][trait])} and nt == {str(N[trait])},
              (sorted(nc)[:3], sorted(nt)[:3]), (T["cases"][trait], N[trait]))
    common_rows = [r for r in rows if r["MarkerID"] == C["id"]]
    ps20 = [fnum(r["p.value"]) for r in rows
            if min(fnum(r["AC_Allele2"]), 2 * N[trait] - fnum(r["AC_Allele2"])) >= 20
            and r["MarkerID"] != C["id"]]
    lam = lambda_gc(ps20)
    check("step2 %s variant: lambda_GC (MAC >= 20)" % trait, 0.85 <= lam <= 1.15,
          "%.3f over %d variants" % (lam, len(ps20)), "0.85-1.15")
    if positive:
        if common_rows:
            r = common_rows[0]
            p, beta = fnum(r["p.value"]), fnum(r["BETA"])
            sign = 1 if r["Allele2"] == C["alt"] else -1
            bar = 1e-6 if trait.startswith("Q") else 1e-4
            check("step2 %s variant: planted common variant" % trait, p < bar and sign * beta > 0,
                  "p %s, BETA %.3f on %s (ALT %s)" % (fmt_p(p), beta, r["Allele2"], C["alt"]),
                  "p < %g, effect on ALT > 0" % bar)
        else:
            check("step2 %s variant: planted common variant" % trait, False, "not in output", "present")
    else:
        pmin = min(fnum(r["p.value"]) for r in rows)
        check("step2 %s variant: no variant below 1e-6" % trait, pmin > 1e-6,
              "min p %s over %d" % (fmt_p(pmin), len(rows)), "> 1e-6")


def group_rows(prefix):
    rows = table(out(prefix + ".txt"))
    return rows or []


def gene_p(rows, gene, group, maf, col="Pvalue_Burden"):
    for r in rows:
        if r["Region"] == gene and r["Group"] == group and abs(fnum(r["max_MAF"]) - maf) < 1e-9:
            return fnum(r[col]), fnum(r.get("BETA_Burden"))
    return float("nan"), float("nan")


def null_gene_min(rows, exclude):
    ps = [fnum(r["Pvalue"]) for r in rows if r["Region"] not in exclude and r["Group"] != "Cauchy"]
    ps = [p for p in ps if not math.isnan(p)]
    return (min(ps), len(ps)) if ps else (float("nan"), 0)


PLANTED = {"GENE_BURDEN", "GENE_FLEX"}


def group_complete(trait, rows):
    genes = {r["Region"] for r in rows}
    check("step2 %s group: every gene tested" % trait, len(genes) == T["n_genes"], len(genes), T["n_genes"])
variant_checks("Q_pos", True)
variant_checks("B_pos", True)
variant_checks("Q_null", False)

for trait, bar in (("Q_pos", 1e-6), ("B_pos", 1e-6)):
    section("s2_%s_group" % trait)
    rows = group_rows(trait + ".group")
    group_complete(trait, rows)
    p, beta = gene_p(rows, "GENE_BURDEN", "pLoF", 0.01)
    check("step2 %s group: GENE_BURDEN pLoF burden" % trait, p < bar and beta > 0,
          "p %s, BETA %.3f" % (fmt_p(p), beta), "p < %g, BETA > 0" % bar)
    pmin, n = null_gene_min(rows, PLANTED)
    check("step2 %s group: no null gene below 1e-6" % trait, pmin > 1e-6, "min p %s over %d" % (fmt_p(pmin), n), "> 1e-6")

for trait in ("Q_null", "Q_perm"):
    section("s2_%s_group" % trait)
    rows = group_rows(trait + ".group")
    group_complete(trait, rows)
    p, _ = gene_p(rows, "GENE_BURDEN", "pLoF", 0.01)
    check("step2 %s group: GENE_BURDEN not significant" % trait, p > 1e-3, "p %s" % fmt_p(p), "> 1e-3")
    pmin, n = null_gene_min(rows, set())
    check("step2 %s group: no gene below 1e-6" % trait, pmin > 1e-6, "min p %s over %d" % (fmt_p(pmin), n), "> 1e-6")

section("s2_Q_pos_group", "s2_Q_pos_group_plink")
for suffix in (".txt", ".txt.singleAssoc.txt"):
    a, b = out("Q_pos.group" + suffix), out("Q_pos.group_plink" + suffix)
    check("step2: group%s from --pgen == --plink (bytes)" % suffix, same_bytes(a, b), "", "identical")

section("s2_Q_pos_flexrv")
flex = group_rows("Q_pos.flexrv_AM")
p_flex, _ = gene_p(flex, "GENE_FLEX", "Cauchy", 0.001) if flex else (float("nan"), 0)
if math.isnan(p_flex):          # the pooled row's max_MAF may not be the run's max MAF
    p_flex = next((fnum(r["Pvalue_Burden"]) for r in flex if r["Region"] == "GENE_FLEX" and r["Group"] == "Cauchy"),
                  float("nan"))
check("step2 Q_pos FlexRV: GENE_FLEX pooled p", p_flex < 1e-6, "p %s" % fmt_p(p_flex), "< 1e-6")
null_flex = [fnum(r["Pvalue_Burden"]) for r in flex if r["Group"] == "Cauchy" and r["Region"] not in PLANTED]
null_flex = [p for p in null_flex if not math.isnan(p)]
check("step2 Q_pos FlexRV: no null gene below 1e-6", null_flex and min(null_flex) > 1e-6,
      "min p %s over %d genes" % (fmt_p(min(null_flex)) if null_flex else "nan", len(null_flex)), "> 1e-6")


def score_advantage(rows):
    """(unweighted p, best score-weighted p) for GENE_FLEX, both without a MAF weight.
    Transform sets are named <score>_<scoreTransform>_<mafTransform>."""
    unw, weighted = float("nan"), []
    for r in rows:
        if r["Region"] != "GENE_FLEX" or r["Group"] == "Cauchy":
            continue
        parts = r.get("Weight", "").split("_")
        if len(parts) < 3 or parts[-1] != "unw":
            continue
        if parts[-2] == "unw":
            unw = fnum(r["Pvalue_Burden"])
        else:
            weighted.append(fnum(r["Pvalue_Burden"]))
    return unw, min(weighted) if weighted else float("nan")


u, wbest = score_advantage(flex)
check("step2 Q_pos FlexRV: score weighting beats unweighted >= 100x (scores reach SAIGE aligned)",
      wbest * 100 < u, "unweighted %s, best weighted %s" % (fmt_p(u), fmt_p(wbest)), "weighted < unweighted / 100")
section("s2_Q_pos_flexrv_shuf")
u, wbest = score_advantage(group_rows("Q_pos.flexrv_AM_shuffled"))
check("step2 Q_pos FlexRV, scores shuffled: advantage gone",
      not math.isnan(u) and not wbest * 100 < u, "unweighted %s, best weighted %s" % (fmt_p(u), fmt_p(wbest)),
      "weighted >= unweighted / 100")
section("s2_Q_pos_group", "s2_Q_pos_flexrv")
grp = group_rows("Q_pos.group")
sg = min((gene_p(grp, "GENE_FLEX", g, 0.001)[0] for g in
          ("damaging_missense_or_protein_altering", "other_missense_or_protein_altering")), default=float("nan"))
check("info: GENE_FLEX best SAIGE-GENE+ missense-mask burden p (MAF 0.001)", True, fmt_p(sg),
      "for comparison with FlexRV %s" % fmt_p(p_flex))

# ------------------------------------------------------------ step 3
# extractNglmm.R's Nglmm is 1' K^-1 1 on the sparse GRM K (entries <= 0.125
# zeroed) over the trait's analysed samples, times 4 Pn (1 - Pn) for a binary
# trait: it does not depend on the fitted variance components. Recompute it
# here, block by block, from the .mtx.


def solve(A, b):
    m = len(A)
    M = [A[r][:] + [b[r]] for r in range(m)]
    for c in range(m):
        piv = max(range(c, m), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        for r in range(m):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [a - f * bb for a, bb in zip(M[r], M[c])]
    return [M[r][m] / M[r][r] for r in range(m)]


def nglmm_direct(trait, binary):
    g = out("step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx")
    ids = [l.strip() for l in open(g + ".sampleIDs.txt") if l.strip()]
    ph = {r["IID"]: r for r in table(os.path.join(IN, "pheno.tsv"))}
    use = {k for k, s in enumerate(ids) if s in ph and ph[s][trait] != "NA"}
    val, adj = {}, {k: set() for k in use}
    with open(g) as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for line in fh:
            i, j, v = line.split()
            i, j, v = int(i) - 1, int(j) - 1, float(v)
            # getsubGRM's drop0(tol = cutoff) keeps entries ABOVE the cutoff, 0.05 here as in steps 1-2
            if i in use and j in use and v > 0.05:
                val[(i, j)] = val[(j, i)] = v
                if i != j:
                    adj[i].add(j)
                    adj[j].add(i)
    seen, total = set(), 0.0
    for s in use:
        if s in seen:
            continue
        comp, stack = [], [s]
        seen.add(s)
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj[x] - seen:
                seen.add(y)
                stack.append(y)
        total += sum(solve([[val.get((a, b), 0.0) for b in comp] for a in comp], [1.0] * len(comp)))
    if binary:
        pn = sum(1 for k in use if ph[ids[k]][trait] == "1") / len(use)
        total *= 4 * pn * (1 - pn)
    return total


section("s3_nglmm")
path = out("neff.csv")
vals = {}
if os.path.exists(path):
    for line in open(path).read().split("\n")[1:]:
        if "," in line:
            k, v = line.split(",", 1)
            vals[k] = fnum(v)
for trait in ("Q_pos", "B_pos", "B_rare"):
    v = vals.get(trait, float("nan"))
    try:
        d = nglmm_direct(trait, trait.startswith("B_"))
    except (OSError, TypeError):
        d = float("nan")
    check("step3 %s: Nglmm = 1'K^-1 1 recomputed from the GRM" % trait, abs(v - d) <= 1e-4 * abs(d),
          "%.4f vs %.4f (N %d)" % (v, d, N[trait]), "equal to 1e-4")

section("s3_nglmm", "s3_nglmm_extmount")
check("step3: absolute paths through SAIGE_EXTRA_MOUNTS give the same csv (bytes)",
      same_bytes(out("neff.csv"), out("neff_extmount.csv")), "", "identical")

# ------------------------------------------------------------ report
width = max(len(r[1]) for r in RESULTS)
for ok, name, obs, exp in RESULTS:
    print("%s  %-*s  %s%s" % ("PASS" if ok else "FAIL", width, name, obs, ("   [expected %s]" % exp) if exp else ""))
for name, why in SKIPPED:
    print("SKIP  %-*s  %s" % (width, name, why))
n_fail = sum(1 for r in RESULTS if not r[0])
print("\n%d checks, %d failed%s" % (len(RESULTS), n_fail, ", %d skipped" % len(SKIPPED) if SKIPPED else ""))
sys.exit(1 if n_fail else 0)
