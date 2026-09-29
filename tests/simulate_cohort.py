#!/usr/bin/env python3
"""Simulate a small cohort with a KNOWN answer, for tests/run_e2e.sh.

Standard-library Python 3 only. Writes, under --out (default tests/work/in):

  array/chr{1..4}.{bed,bim,fam}  "genotyping array" for step 0 (common SNPs for
                                 the GRM, MAC 10-20 SNPs for the variance ratio);
                                 FID is the family, so FID != IID
  exome/chr7.{bed,bim,fam}       "exome" for step 2: genes of rare variants
  group.chr7.txt                 BRaVa-style SAIGE group file (var + anno lines)
  am_dummy.tsv.gz                AlphaMissense-shaped per-variant scores
  pheno.tsv, pheno_mf.tsv        phenotypes and covariates
  sample_ids.txt                 the step 0 / step 1 sample include list
  truth.json                     what the pipeline should find

Design (all numbers fixed by --seed):
  * 500 families of 4 (two unrelated founders, two full sibs) = 2,000 samples;
    loci unlinked. First-degree pairs have kinship 0.25 (GRM 0.5): the known
    answer for the sparse GRM.
  * Sample-ID bookkeeping is deliberately untidy: the phenotype file drops 30
    genotyped samples, adds 50 that were never genotyped and is shuffled; the
    include list drops 100 and adds 10 phantoms; each trait has ~5% NA. The
    expected N per trait is the intersection, so a wrong join shows as a
    wrong N.
  * Planted effects (the positive controls), in SD units of the liability:
      GENE_BURDEN  8 pLoF variants, MAF 5e-4..3.5e-3, beta 0.8 per allele
      GENE_FLEX    60 missense with AlphaMissense > 0.85, MAF 5e-5..2.5e-4,
                   beta 1.5; 60 more with AlphaMissense < 0.3 and no effect;
                   the damaging/other label is a coin flip, so only the
                   score tells them apart. Half the carriers are null, so a
                   score-weighted burden should beat the unweighted one
                   (z ~ sqrt(2) larger), and must not once the scores are
                   shuffled
      COMMON       one variant in GENE_0100, MAF 0.3, beta 0.25 on the ALT allele
  * Traits:
      Q_pos, B_pos  covariates + polygenic (h2 0.2) + planted effects;
                    B_pos is the top 20% of its liability
      Q_null        the same without planted effects      (negative control)
      Q_perm        Q_pos shuffled across samples          (negative control)
      B_rare        top 2.5%: about 45 cases, below the image's 100-case gate
      Q_female      Q_null, non-missing for sex 0 only (the --sex path)
"""

import argparse
import gzip
import json
import math
import os
import random

N_FAM = 500
FAM_SIZE = 4                      # father, mother, two children
N = N_FAM * FAM_SIZE
N_FOUNDER_HAPS = N_FAM * 4        # 2 founders x 2 haplotypes per family
BASES = "ACGT"


def binomial(rng, n, p):
    f = getattr(rng, "binomialvariate", None)     # Python >= 3.12
    if f is not None:
        return f(n, p)
    return sum(1 for _ in range(n) if rng.random() < p)


def simulate_variant(rng, maf):
    """ALT dosage (0/1/2) for every sample; samples 4f..4f+3 are family f."""
    k = binomial(rng, N_FOUNDER_HAPS, maf)
    haps = rng.sample(range(N_FOUNDER_HAPS), k)
    carried = bytearray(N_FOUNDER_HAPS)           # haplotypes 4f,4f+1 father; 4f+2,4f+3 mother
    g = bytearray(N)
    fams = set()
    for h in haps:
        carried[h] = 1
        f = h // 4
        g[4 * f + (h % 4) // 2] += 1
        fams.add(f)
    for f in fams:
        for child in (2, 3):
            g[4 * f + child] = (carried[4 * f + rng.getrandbits(1)]
                                + carried[4 * f + 2 + rng.getrandbits(1)])
    return g


# .bed: SNP-major, 2 bits per sample, first sample in the low bits. With the
# .bim's A1 = ALT and A2 = REF: 00 = ALT/ALT, 10 = het, 11 = REF/REF.
CODE = {2: 0b00, 1: 0b10, 0: 0b11}
PACK = [0] * 81
for a in range(3):
    for b in range(3):
        for c in range(3):
            for d in range(3):
                PACK[a + 3 * b + 9 * c + 27 * d] = (
                    CODE[a] | CODE[b] << 2 | CODE[c] << 4 | CODE[d] << 6)


def pack(g):
    return bytes(PACK[g[i] + 3 * g[i + 1] + 9 * g[i + 2] + 27 * g[i + 3]]
                 for i in range(0, N, 4))


class Plink1Writer:
    def __init__(self, prefix, iids, family_fid=False):
        self.bed = open(prefix + ".bed", "wb")
        self.bed.write(bytes([0x6C, 0x1B, 0x01]))
        self.bim = open(prefix + ".bim", "w")
        with open(prefix + ".fam", "w") as fam:
            for i, iid in enumerate(iids):
                sex = 1 if SEX[i] == 1 else 2        # plink: 1 male, 2 female
                fid = "FAM%03d" % (i // FAM_SIZE) if family_fid else iid
                fam.write("%s %s 0 0 %d -9\n" % (fid, iid, sex))

    def add(self, chrom, vid, pos, ref, alt, g):
        self.bim.write("%s\t%s\t0\t%d\t%s\t%s\n" % (chrom, vid, pos, alt, ref))
        self.bed.write(pack(g))

    def close(self):
        self.bed.close()
        self.bim.close()


def alleles(rng):
    ref = rng.choice(BASES)
    alt = rng.choice([b for b in BASES if b != ref])
    return ref, alt


def positions(rng, n, start, span):
    return sorted(rng.sample(range(start, start + span), n))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="tests/work/in")
    ap.add_argument("--seed", type=int, default=20260928)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = args.out
    for d in ("array", "exome"):
        os.makedirs(os.path.join(out, d), exist_ok=True)

    iids = ["S%04d" % i for i in range(N)]
    global SEX
    SEX = [rng.getrandbits(1) for _ in range(N)]     # 0 = female, 1 = male
    age = [rng.uniform(40, 70) for _ in range(N)]
    pcs = [[rng.gauss(0, 1) for _ in range(N)] for _ in range(4)]

    # ------------------------------------------------------------ array data
    # 7,000 common SNPs (GRM, polygenic score) + 5,000 at MAF 0.002-0.005,
    # which lands most of them at MAC 10-20 (the step 0 variance-ratio bin).
    specs = ([("common", rng.uniform(0.05, 0.5)) for _ in range(7000)]
             + [("low", rng.uniform(0.002, 0.005)) for _ in range(5000)])
    rng.shuffle(specs)
    by_chrom = {c: [] for c in (1, 2, 3, 4)}
    for s in specs:
        by_chrom[rng.choice((1, 2, 3, 4))].append(s)
    poly = [0.0] * N
    n_poly, h2 = 1000, 0.2
    poly_left = n_poly
    for c, cs in by_chrom.items():
        # FID = family, not IID: step 0 must match --sampleIDs on IID alone
        w = Plink1Writer(os.path.join(out, "array", "chr%d" % c), iids, family_fid=True)
        for (kind, maf), pos in zip(cs, positions(rng, len(cs), 1_000_000, 100_000_000)):
            g = simulate_variant(rng, maf)
            ref, alt = alleles(rng)
            w.add(c, "chr%d:%d:%s:%s" % (c, pos, ref, alt), pos, ref, alt, g)
            if kind == "common" and poly_left > 0:
                poly_left -= 1
                beta = rng.gauss(0, math.sqrt(h2 / n_poly))
                sd = math.sqrt(2 * maf * (1 - maf))
                for i in range(N):
                    poly[i] += beta * (g[i] - 2 * maf) / sd
        w.close()

    # ------------------------------------------------------------ exome data
    ANNOS = [("pLoF", 0.10), ("damaging_missense_or_protein_altering", 0.25),
             ("other_missense_or_protein_altering", 0.35), ("synonymous", 0.25),
             ("non_coding", 0.05)]
    MISSENSE = ("damaging_missense_or_protein_altering",
                "other_missense_or_protein_altering")

    def draw_anno():
        u, acc = rng.random(), 0.0
        for a, p in ANNOS:
            acc += p
            if u < acc:
                return a
        return ANNOS[-1][0]

    genes = ["GENE_%04d" % i for i in range(300)]
    genes[40], genes[160] = "GENE_BURDEN", "GENE_FLEX"
    w = Plink1Writer(os.path.join(out, "exome", "chr7"), iids)
    group_lines, am_rows = [], []
    burden_dose = [0] * N
    flex_dose = [0] * N
    common_dose, common = None, None
    truth_variants = {"GENE_BURDEN_pLoF": [], "GENE_FLEX_causal": [], "GENE_FLEX_null": []}
    start = 10_000_000
    for gene in genes:
        # (annotation, maf, role, am score or None)
        vs = []
        if gene == "GENE_BURDEN":
            vs += [("pLoF", rng.uniform(5e-4, 3.5e-3), "burden", None) for _ in range(8)]
            vs += [(a, m, None, None) for a, m in
                   ((draw_anno(), math.exp(rng.uniform(math.log(2.5e-4), math.log(1e-2))))
                    for _ in range(35)) if a != "pLoF"]
        elif gene == "GENE_FLEX":
            # MAF well under FlexRV's 0.001: a family carries a founder copy
            # about twice, so a planted 5e-4 lands above the cutoff often
            vs += [(rng.choice(MISSENSE), rng.uniform(5e-5, 2.5e-4), "flex", rng.uniform(0.85, 1.0))
                   for _ in range(60)]
            vs += [(rng.choice(MISSENSE), rng.uniform(5e-5, 2.5e-4), "flexnull", rng.uniform(0.0, 0.3))
                   for _ in range(60)]
            vs += [("synonymous", math.exp(rng.uniform(math.log(2.5e-4), math.log(1e-2))), None, None)
                   for _ in range(10)]
        else:
            vs += [(draw_anno(), math.exp(rng.uniform(math.log(2.5e-4), math.log(1e-2))), None, None)
                   for _ in range(rng.randint(30, 50))]
        vs += [(rng.choice(("synonymous", "other_missense_or_protein_altering")),
                rng.uniform(0.05, 0.45), None, None) for _ in range(3)]
        if gene == "GENE_0100":
            vs.append(("synonymous", 0.3, "common", None))
        rng.shuffle(vs)
        ids, annos = [], []
        for (anno, maf, role, am), pos in zip(vs, positions(rng, len(vs), start, 50_000)):
            g = simulate_variant(rng, maf)
            ref, alt = alleles(rng)
            vid = "chr7:%d:%s:%s" % (pos, ref, alt)
            w.add(7, vid, pos, ref, alt, g)
            ids.append(vid)
            annos.append(anno)
            if role == "burden":
                truth_variants["GENE_BURDEN_pLoF"].append(vid)
                for i in range(N):
                    burden_dose[i] += g[i]
            elif role == "flex":
                truth_variants["GENE_FLEX_causal"].append(vid)
                for i in range(N):
                    flex_dose[i] += g[i]
            elif role == "flexnull":
                truth_variants["GENE_FLEX_null"].append(vid)
            elif role == "common":
                common_dose, common = g, {"id": vid, "ref": ref, "alt": alt, "gene": gene}
            if anno in MISSENSE:
                if am is None:
                    am = (rng.uniform(0.564, 1.0) if anno == MISSENSE[0]
                          else rng.uniform(0.0, 0.34))
                # 5% of ordinary missense have no AlphaMissense row (gene-mean imputation)
                if role in ("flex", "flexnull") or rng.random() > 0.05:
                    am_rows.append(("chr7", pos, ref, alt, am))
        group_lines.append("%s var %s" % (gene, " ".join(ids)))
        group_lines.append("%s anno %s" % (gene, " ".join(annos)))
        start += 200_000
    w.close()
    with open(os.path.join(out, "group.chr7.txt"), "w") as fh:
        fh.write("\n".join(group_lines) + "\n")

    # decoy rows the join must ignore: unmatched positions and another chromosome
    for _ in range(300):
        ref, alt = alleles(rng)
        am_rows.append(("chr7", rng.randint(1, 9_000_000), ref, alt, rng.random()))
    ref, alt = alleles(rng)
    am_rows.append(("chr8", 12345, ref, alt, 0.5))
    am_rows.sort(key=lambda r: (r[0], r[1]))
    with gzip.open(os.path.join(out, "am_dummy.tsv.gz"), "wt") as fh:
        fh.write("# dummy AlphaMissense table for tests/run_e2e.sh\n")
        fh.write("#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\t"
                 "protein_variant\tam_pathogenicity\tam_class\n")
        for c, pos, ref, alt, s in am_rows:
            cls = ("likely_pathogenic" if s > 0.564 else
                   "ambiguous" if s >= 0.34 else "likely_benign")
            fh.write("%s\t%d\t%s\t%s\thg38\tP00000\tENST00000000001.1\tA1B\t%.4f\t%s\n"
                     % (c, pos, ref, alt, s, cls))

    # ------------------------------------------------------------ phenotypes
    age_mean = sum(age) / N
    age_sd = math.sqrt(sum((a - age_mean) ** 2 for a in age) / N)
    base = [0.3 * (age[i] - age_mean) / age_sd + 0.2 * SEX[i] + poly[i] for i in range(N)]
    c_maf = 0.3
    planted = [0.8 * burden_dose[i] + 1.5 * flex_dose[i] + 0.25 * (common_dose[i] - 2 * c_maf)
               for i in range(N)]
    noise = lambda: rng.gauss(0, math.sqrt(0.7))
    q_null = [base[i] + noise() for i in range(N)]
    q_pos = [base[i] + planted[i] + noise() for i in range(N)]
    lb = [base[i] + 1.5 * planted[i] + noise() for i in range(N)]
    lr = [base[i] + noise() for i in range(N)]

    def top(vals, frac):
        cut = sorted(vals)[int(len(vals) * (1 - frac))]
        return [1 if v >= cut else 0 for v in vals]

    b_pos, b_rare = top(lb, 0.20), top(lr, 0.025)

    # the file: 30 genotyped samples dropped, 50 never-genotyped added, shuffled
    geno = set(iids)
    dropped = set(rng.sample(iids, 30))
    rows = []
    for i in range(N):
        if iids[i] in dropped:
            continue
        rows.append(dict(IID=iids[i], age=age[i], sex=SEX[i], pcs=[p[i] for p in pcs],
                         Q_null=q_null[i], Q_pos=q_pos[i], B_pos=b_pos[i], B_rare=b_rare[i]))
    for j in range(50):
        s = rng.getrandbits(1)
        rows.append(dict(IID="X%04d" % j, age=rng.uniform(40, 70), sex=s,
                         pcs=[rng.gauss(0, 1) for _ in range(4)], Q_null=rng.gauss(0, 1),
                         Q_pos=rng.gauss(0, 1), B_pos=rng.getrandbits(1), B_rare=0))
    perm = [r["Q_pos"] for r in rows]
    rng.shuffle(perm)
    for r, v in zip(rows, perm):
        r["Q_perm"] = v
        r["Q_female"] = r["Q_null"] if r["sex"] == 0 else None
    traits = ["Q_null", "Q_pos", "Q_perm", "B_pos", "B_rare", "Q_female"]
    for t in traits:                                   # ~5% missing per trait
        for r in rows:
            if rng.random() < 0.05:
                r[t] = None
    rng.shuffle(rows)

    def fmt(v):
        if v is None:
            return "NA"
        return str(v) if isinstance(v, int) else "%.6f" % v

    cols = ["IID", "age", "age2", "sex", "PC1", "PC2", "PC3", "PC4"] + traits
    with open(os.path.join(out, "pheno.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = ([r["IID"], "%.3f" % r["age"], "%.3f" % (r["age"] ** 2), str(r["sex"])]
                    + ["%.6f" % p for p in r["pcs"]] + [fmt(r[t]) for t in traits])
            fh.write("\t".join(vals) + "\n")

    # a few rows coded M/F, for the --sex checks (they never reach SAIGE)
    with open(os.path.join(out, "pheno_mf.tsv"), "w") as fh:
        fh.write("IID\tsex\tQ_f\tQ_both\n")
        for i in range(20):
            s = "F" if i % 2 == 0 else "M"
            fh.write("%s\t%s\t%s\t%.3f\n" % (iids[i], s, "%.3f" % rng.gauss(0, 1) if s == "F" else "NA",
                                             rng.gauss(0, 1)))

    # include list: 100 genotyped samples left out, 10 phantoms added
    left_out = set(rng.sample(iids, 100))
    include = [s for s in iids if s not in left_out] + ["Z%04d" % j for j in range(10)]
    rng.shuffle(include)
    with open(os.path.join(out, "sample_ids.txt"), "w") as fh:
        fh.write("\n".join(include) + "\n")

    # ------------------------------------------------------------ truth
    kept = geno - left_out                            # step 0 --keep: GRM and VR samples
    expected_n, cases = {}, {}
    for t in traits:
        ids = [r["IID"] for r in rows if r[t] is not None and r["IID"] in kept]
        expected_n[t] = len(ids)
        if t.startswith("B_"):
            cases[t] = sum(1 for r in rows if r[t] == 1 and r["IID"] in kept)
    pairs = []
    for f in range(N_FAM):
        fam = [iids[4 * f + k] for k in range(4)]
        for a, b in ((0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
            if fam[a] in kept and fam[b] in kept:
                pairs.append([fam[a], fam[b]])
    truth = dict(
        n_genotyped=N, n_kept=len(kept), expected_n=expected_n, cases=cases,
        first_degree_pairs=pairs, common=common, genes=dict(burden="GENE_BURDEN", flex="GENE_FLEX"),
        variants=truth_variants, n_genes=len(genes), seed=args.seed,
    )
    with open(os.path.join(out, "truth.json"), "w") as fh:
        json.dump(truth, fh, indent=1)
    print("simulated %d samples (%d kept), %d genes; N per trait %s; cases %s"
          % (N, len(kept), len(genes), expected_n, cases))


if __name__ == "__main__":
    main()
