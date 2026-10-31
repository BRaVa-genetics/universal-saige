#!/usr/bin/env python3
"""A sparse GRM for SAIGE from PC-Relate kinship, computed natively (numpy, blocked BLAS).

Why: KING-robust is biased on admixed pairs of differing ancestry -- on AoU amr at 0.05,
third-degree links joined 37k people into one indefinite component (LEDGER #183) -- and the
alternative to a sparse GRM is a dense fit that costs a thousandfold more. PC-Relate
(Conomos et al. 2016) estimates kinship CONDITIONAL on ancestry: each marker is regressed on
the PCs over an unrelated set, giving individual-specific allele frequencies mu_im, and

    kinship_ij = sum_m (g_im - 2 mu_im)(g_jm - 2 mu_jm) / (4 sum_m sqrt(mu_im (1-mu_im) mu_jm (1-mu_jm)))

over the markers both samples carry with 0 < mu < 1. GENESIS's pcrelate() computes this in an
R-level block loop and returns every pair (3.2e9 at N = 80k: a day). The estimator is two
matrix products, R R' and A A' with R = (G - 2 mu) and A = sqrt(mu (1 - mu)), masked entries
zero -- so it is streamed here in marker blocks against sample tiles with numpy's threaded
BLAS in float32, only the pairs at or above --cutoff (on SAIGE's GRM scale, 2 x kinship) are
kept, and the output is the Matrix Market file SAIGE reads, same conventions as
king_to_mtx.py convert (off-diagonal = 2 x kinship, symmetric coordinate, lower triangle,
1-based) with the diagonal 1 (default) or 1 + the PC-Relate inbreeding estimate.

  run       --bed PREFIX (PLINK 1 .bed/.bim/.fam, SNP-major) or --pgen PREFIX (PLINK 2
            .pgen/.pvar/.psam, hard calls; needs the Pgenlib python package, `pip install Pgenlib`,
            or convert with `plink2 --pfile X --make-bed`) --pcs FILE (plink2 .eigenvec or any
            header-first table with an IID / sample.id / ID column and PC1.. columns; PC-AiR's
            output is the recommended source) --out PREFIX [--npcs K] [--unrelated IDS: the
            samples the regressions are fitted on; default every sample with PCs] [--cutoff 0.05]
            [--maf 0.01: markers below this allele frequency on the regression set are skipped]
            [--diag one|inbreeding] [--block-markers 2048] [--block-samples 4096] [--threads T]
            Writes <out>.mtx, <out>.mtx.sampleIDs.txt (fam order, samples with PCs) and
            <out>.report.tsv (counts: samples, markers used, pairs kept, pairs per cutoff band).
            No pair count below --mask is printed (real data).
  selftest  a fabricated two-population cohort with admixed individuals and planted relatives
            (parent-offspring and full sibs at 0.25, half sibs at 0.125, first cousins at
            0.0625, cross-population parent-offspring, one inbred child of full sibs at f 0.25),
            PCs from the unrelated founders, every claim checked: relatives recovered (positive),
            unrelated pairs within and across populations at 0 with nothing above the cutoff
            (negative), the SAME estimator with the PCs removed misreads the cross-population
            pairs (the correction is load-bearing), a within-marker permutation destroys every
            pair (null), tiny blocks reproduce one block (two ways agree), 2 % missing genotypes
            change nothing material, and the .mtx round-trips.

Memory at N samples: two float32 accumulators of block_samples x N, sized by default to
--memory-gb (16 GB at N = 80k gives a tile of 25k, 4 passes), plus the marker block
(N x block_markers x 4 B x 3). The genotypes are decoded and prepared once per sample tile, so
a larger budget means fewer passes. Flops 4 N^2 M: at N = 80k, M = 100k that is ~2.6e15, tens of
minutes on a 32-core BLAS. --threads sets OPENBLAS/OMP/MKL_NUM_THREADS before numpy loads; do
NOT source harness/env.sh first (it pins one thread for the goldens).

Usage:
  pcrelate_to_mtx.py run --bed pruned --pcs pcair.tsv --npcs 10 --unrelated unrel.txt --out cohort_pcrelate
  pcrelate_to_mtx.py run --pgen pruned --pcs pcair.tsv --npcs 10 --unrelated unrel.txt --out cohort_pcrelate
  pcrelate_to_mtx.py selftest     (the PGEN arm runs when Pgenlib and plink2 are available, else says so)
"""
import argparse, gzip, os, sys, tempfile, time


def set_threads(t):
    for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = str(t)


def openr(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


# ----------------------------------------------------------------------------- inputs

def read_fam(prefix):
    with open(prefix + ".fam") as fh:
        ids = [l.split()[1] for l in fh if l.strip()]
    if len(set(ids)) != len(ids):
        raise SystemExit(".fam has a repeated IID (SAIGE keys on IID alone)")
    return ids


def count_bim(prefix):
    with open(prefix + ".bim") as fh:
        return sum(1 for l in fh if l.strip())


def read_pcs(path, npcs=None):
    """{iid: np.array(PCs)} from a header-first table (plink2 .eigenvec, PC-AiR, ...)."""
    import numpy as np
    with openr(path) as fh:
        head = fh.readline().lstrip("#").split()
        up = [h.upper() for h in head]
        idc = next((i for i, h in enumerate(up) if h in ("IID", "SAMPLE.ID", "SAMPLE_ID", "ID", "SAMPLEID")), None)
        if idc is None:
            raise SystemExit("--pcs: no IID / sample.id / ID column in the header: " + " ".join(head))
        pcc = [i for i, h in enumerate(up) if h.startswith("PC") and h[2:].isdigit()]
        pcc.sort(key=lambda i: int(up[i][2:]))
        if not pcc:
            raise SystemExit("--pcs: no PC1.. columns in the header")
        if npcs is not None:
            if npcs > len(pcc):
                raise SystemExit(f"--npcs {npcs} but the file has {len(pcc)} PC columns")
            pcc = pcc[:npcs]
        out = {}
        for l in fh:
            f = l.split()
            if not f:
                continue
            out[f[idc]] = np.array([float(f[i]) for i in pcc], dtype=np.float64)
    if not out:
        raise SystemExit("--pcs: no rows")
    return out


def read_ids(p):
    with openr(p) as fh:
        return [l.strip().split()[0] for l in fh if l.strip()]


class Bed:
    """PLINK 1 .bed, SNP-major; decode a block of markers to an int8 (markers x samples) array
    with 0/1/2 copies of A2 and -1 missing (2-bit codes 00 hom A1, 01 missing, 10 het, 11 hom A2)."""

    def __init__(self, prefix, n_samples, n_markers):
        import numpy as np
        self.path = prefix + ".bed"
        self.n = n_samples
        self.m = n_markers
        self.bpm = (n_samples + 3) // 4
        with open(self.path, "rb") as fh:
            magic = fh.read(3)
        if magic != b"\x6c\x1b\x01":
            raise SystemExit(f"{self.path}: not a SNP-major PLINK 1 .bed")
        size = os.path.getsize(self.path)
        if size != 3 + self.bpm * self.m:
            raise SystemExit(f"{self.path}: {size} bytes, expected 3 + {self.bpm} x {self.m}")
        # byte -> 4 codes, sample order within the byte = low bits first
        lut = np.empty((256, 4), dtype=np.int8)
        code = np.array([0, -1, 1, 2], dtype=np.int8)   # 00, 01, 10, 11
        for b in range(256):
            for k in range(4):
                lut[b, k] = code[(b >> (2 * k)) & 3]
        self.lut = lut
        self.fh = open(self.path, "rb")

    def block(self, start, count):
        import numpy as np
        self.fh.seek(3 + start * self.bpm)
        raw = np.frombuffer(self.fh.read(count * self.bpm), dtype=np.uint8).reshape(count, self.bpm)
        g = self.lut[raw].reshape(count, self.bpm * 4)[:, : self.n]
        return g   # markers x samples, int8


def read_psam(prefix):
    """IIDs from a .psam (header #FID IID ... or #IID ...)."""
    ids = []
    with open(prefix + ".psam") as fh:
        head = None
        for l in fh:
            if l.startswith("##"):
                continue
            if head is None:
                if not l.startswith("#"):
                    raise SystemExit(prefix + ".psam: no header line")
                head = l.lstrip("#").split()
                up = [h.upper() for h in head]
                idc = up.index("IID") if "IID" in up else 0
                continue
            f = l.split()
            if f:
                ids.append(f[idc])
    if len(set(ids)) != len(ids):
        raise SystemExit(".psam has a repeated IID (SAIGE keys on IID alone)")
    return ids


def count_pvar(prefix):
    with open(prefix + ".pvar") as fh:
        return sum(1 for l in fh if l.strip() and not l.startswith("#"))


class Pgen:
    """PLINK 2 .pgen hard calls through pgenlib: 0/1/2 copies of ALT, -1 missing, markers x samples."""

    def __init__(self, prefix, n_samples, n_markers):
        try:
            import pgenlib
        except ImportError:
            raise SystemExit("--pgen needs the Pgenlib python package (pip install Pgenlib), "
                             "or convert with: plink2 --pfile PREFIX --make-bed")
        self.r = pgenlib.PgenReader((prefix + ".pgen").encode())
        self.n = n_samples
        self.m = n_markers
        if self.r.get_raw_sample_ct() != n_samples or self.r.get_variant_ct() != n_markers:
            raise SystemExit(f"{prefix}.pgen: {self.r.get_raw_sample_ct()} samples x {self.r.get_variant_ct()} "
                             f"variants, but .psam/.pvar say {n_samples} x {n_markers}")

    def block(self, start, count):
        import numpy as np
        buf = np.empty((count, self.n), dtype=np.int8)
        self.r.read_range(start, start + count, buf)
        buf[buf < 0] = -1          # pgenlib's missing is -9
        return buf


def open_genotypes(a):
    """(ids, M, reader) for --bed or --pgen"""
    if a.pgen:
        ids = read_psam(a.pgen)
        M = count_pvar(a.pgen)
        return ids, M, Pgen(a.pgen, len(ids), M)
    ids = read_fam(a.bed)
    M = count_bim(a.bed)
    return ids, M, Bed(a.bed, len(ids), M)


# ----------------------------------------------------------------------------- estimator

def pcrelate_blocks(get_block, n_markers, X, train, cutoff, maf=0.01, block_markers=2048,
                    block_samples=4096, log=print, keep_all=False):
    """PC-Relate over markers delivered by get_block(start, count) -> int8 (count x N).

    X: N x (k+1) design (intercept + PCs), train: bool mask of the regression samples.
    Returns (pairs, diag_f, n_used) where pairs is a dict {(i, j): 2 * kinship} for i > j at
    or above cutoff (or every pair when keep_all), diag_f the per-sample inbreeding estimate.
    Two ways of blocking must give the same numbers: the selftest checks it.
    """
    import numpy as np
    N = X.shape[0]
    Xt = X[train]
    XtX_inv = np.linalg.pinv(Xt.T @ Xt)
    proj = (XtX_inv @ Xt.T).astype(np.float64)   # (k+1) x n_train

    # per-sample inbreeding accumulators (one pass, independent of the tile loop)
    f_num = np.zeros(N, dtype=np.float64)
    f_den = np.zeros(N, dtype=np.float64)

    Xf = X.astype(np.float32)

    def prepared(start, count):
        """R, A (float32, N x count) for one marker block, plus the markers used.

        Everything after the regression coefficients is float32 and in place: this runs once per
        marker block per sample tile, and measured on 20k x 50k it cost as much as the two GEMMs
        when written with float64 temporaries (4.1 s against 3.9 s per block)."""
        g = get_block(start, count).T                                   # N x count int8
        miss = g < 0
        gt = g[train].astype(np.float32)
        mt = miss[train]
        cnt = (~mt).sum(axis=0)
        gt[mt] = 0
        mean_t = gt.sum(axis=0) / np.maximum(cnt, 1)
        p_t = mean_t / 2
        use = (cnt > 0) & (p_t >= maf) & (p_t <= 1 - maf)
        gt[mt] = np.broadcast_to(mean_t, gt.shape)[mt]                  # mean-imputed for the fit
        beta = (proj @ gt.astype(np.float64)).astype(np.float32)        # (k+1) x count
        del gt
        mu = Xf @ beta                                                  # N x count, 2 mu so far
        mu *= 0.5
        valid = ~miss
        valid &= mu > 0
        valid &= mu < 1
        valid &= use[None, :]
        R = g.astype(np.float32)
        R -= 2 * mu
        R[~valid] = 0
        A = mu * (1 - mu)
        np.sqrt(A, out=A)
        A[~valid] = 0
        return R, A, g, mu, valid, int(use.sum())

    pairs = {}
    n_used_total = None
    tiles = list(range(0, N, block_samples))
    t0 = time.time()
    for ti, ts in enumerate(tiles):
        te = min(ts + block_samples, N)
        num = np.zeros((te - ts, N - ts), dtype=np.float32)   # tile rows x columns from ts on
        den = np.zeros((te - ts, N - ts), dtype=np.float32)
        n_used = 0
        for ms in range(0, n_markers, block_markers):
            mc = min(block_markers, n_markers - ms)
            R, A, g, mu, valid, used = prepared(ms, mc)
            n_used += used
            num += R[ts:te] @ R[ts:].T
            den += A[ts:te] @ A[ts:].T
            if ti == 0:
                # inbreeding: f_i = sum (g^2 - (1 + 2 mu) g + 2 mu^2) / sum 2 mu (1 - mu)
                gm = np.where(valid, g, 0).astype(np.float64)
                mm = np.where(valid, mu, 0).astype(np.float64)
                f_num += (gm * gm - (1 + 2 * mm) * gm + 2 * mm * mm).sum(axis=1)
                f_den += (2 * mm * (1 - mm)).sum(axis=1)
        if n_used_total is None:
            n_used_total = n_used
        with np.errstate(divide="ignore", invalid="ignore"):
            kin = np.where(den > 0, num / (4 * den), 0)
        grm = 2 * kin
        # lower triangle only: row i (global ts + r) against column j (global ts + c), i > j
        rows, cols = np.nonzero(grm >= cutoff) if not keep_all else np.nonzero(np.ones_like(grm, dtype=bool))
        for r, c in zip(rows.tolist(), cols.tolist()):
            i, j = ts + r, ts + c
            if i == j:
                continue
            if i < j:
                i, j = j, i
            if (i, j) in pairs:
                continue
            pairs[(i, j)] = float(grm[r, c])
        log(f"tile {ti + 1}/{len(tiles)}: samples {ts}..{te - 1}, pairs kept so far {len(pairs)}, "
            f"{time.time() - t0:.0f}s")
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(f_den > 0, f_num / f_den, 0.0)
    return pairs, f, n_used_total


# ----------------------------------------------------------------------------- output

def write_mtx(path, n, ent):
    with open(path, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate real symmetric\n")
        fh.write(f"{n} {n} {len(ent)}\n")
        for (i, j) in sorted(ent):
            fh.write(f"{i + 1} {j + 1} {ent[(i, j)]!r}\n")


def read_mtx(path):
    ent = {}
    with openr(path) as fh:
        head = fh.readline()
        if "symmetric" not in head:
            raise SystemExit(f"{path}: not a symmetric Matrix Market file")
        n, n2, nnz = (int(x) for x in fh.readline().split())
        for l in fh:
            i, j, v = l.split()
            ent[(int(i) - 1, int(j) - 1)] = float(v)
    if len(ent) != nnz:
        raise SystemExit(f"{path}: header says {nnz} entries, file has {len(ent)}")
    return n, ent


def mask_count(k, mask):
    return f"<{mask}" if 0 < k < mask else str(k)


def run(a):
    set_threads(a.threads)
    import numpy as np
    t0 = time.time()
    fam, M, bed = open_genotypes(a)
    pcs = read_pcs(a.pcs, a.npcs)
    keep = [k for k, s in enumerate(fam) if s in pcs]
    if not keep:
        raise SystemExit("no .fam sample has PCs")
    ids = [fam[k] for k in keep]
    N = len(ids)
    X = np.column_stack([np.ones(N), np.vstack([pcs[s] for s in ids])])
    if a.unrelated:
        unrel = set(read_ids(a.unrelated))
        train = np.array([s in unrel for s in ids])
        missing = len(unrel - set(ids))
        if train.sum() < X.shape[1] * 10:
            raise SystemExit(f"--unrelated leaves {train.sum()} regression samples for {X.shape[1]} coefficients")
    else:
        train = np.ones(N, dtype=bool)
        missing = 0
    keep_arr = np.array(keep)

    def get_block(start, count):
        return bed.block(start, count)[:, keep_arr]

    if a.block_samples is None:
        a.block_samples = int(max(256, min(N, a.memory_gb * 1e9 / (8.0 * N))))
    print(f"sample tile {a.block_samples} ({-(-N // a.block_samples)} pass(es) over the genotypes, "
          f"accumulators {8.0 * a.block_samples * N / 1e9:.1f} GB)")
    print(f"samples: {N} of {len(fam)} in the {'.psam' if a.pgen else '.fam'} have PCs; regression set {int(train.sum())}"
          + (f" ({missing} --unrelated IDs not in the data)" if missing else "")
          + f"; PCs {X.shape[1] - 1}; markers {M}", flush=True)
    pairs, f, n_used = pcrelate_blocks(get_block, M, X, train, a.cutoff, a.maf, a.block_markers,
                                       a.block_samples, log=lambda s: print(s, flush=True))
    ent = dict(pairs)
    for k in range(N):
        ent[(k, k)] = 1.0 + float(f[k]) if a.diag == "inbreeding" else 1.0
    write_mtx(a.out + ".mtx", N, ent)
    with open(a.out + ".mtx.sampleIDs.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    vals = np.array(list(pairs.values())) if pairs else np.zeros(0)
    bands = [(a.cutoff, 0.125), (0.125, 0.25), (0.25, 0.4), (0.4, 0.75), (0.75, 10)]
    rows = [("samples", N), ("markers_total", M), ("markers_used", n_used), ("pcs", X.shape[1] - 1),
            ("regression_samples", int(train.sum())), ("cutoff_grm_scale", a.cutoff),
            ("pairs_kept", mask_count(len(pairs), a.mask)), ("diag", a.diag),
            ("inbreeding_median", f"{float(np.median(f)):.4g}"),
            ("inbreeding_max", f"{float(f.max()):.4g}" if N else "NA")]
    for lo, hi in bands:
        rows.append((f"pairs_{lo}_{hi}", mask_count(int(((vals >= lo) & (vals < hi)).sum()), a.mask)))
    rows.append(("seconds", f"{time.time() - t0:.0f}"))
    with open(a.out + ".report.tsv", "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")
            print(f"{k:24s} {v}")
    print(f"wrote {a.out}.mtx ({len(ent)} entries), {a.out}.mtx.sampleIDs.txt, {a.out}.report.tsv")


# ----------------------------------------------------------------------------- selftest

def simulate(rng, n_pop=120, n_adm=120, n_markers=16000, fst=0.12):
    """Two Balding-Nichols populations, admixed founders, planted relatives from haplotypes.

    Returns g (N x M int8), labels, and the planted pairs {(i, j): expected kinship}.
    """
    import numpy as np
    p_anc = rng.uniform(0.1, 0.9, n_markers)
    a = p_anc * (1 - fst) / fst
    b = (1 - p_anc) * (1 - fst) / fst
    p1 = np.clip(rng.beta(a, b), 0.01, 0.99)
    p2 = np.clip(rng.beta(a, b), 0.01, 0.99)
    # founders: pop1, pop2, admixed (ancestry proportion uniform)
    anc = np.concatenate([np.ones(n_pop), np.zeros(n_pop), rng.uniform(0, 1, n_adm)])
    nf = len(anc)
    pf = anc[:, None] * p1[None, :] + (1 - anc[:, None]) * p2[None, :]
    hap = (rng.uniform(size=(nf, 2, n_markers)) < pf[:, None, :]).astype(np.int8)   # founder haplotypes
    labels = ["pop1"] * n_pop + ["pop2"] * n_pop + ["adm"] * n_adm
    haps = [hap[i] for i in range(nf)]
    expected = {}

    def child(pa, pb):
        """one haplotype from each parent by independent transmission (unlinked loci)"""
        ha = haps[pa][rng.integers(0, 2, n_markers), np.arange(n_markers)]
        hb = haps[pb][rng.integers(0, 2, n_markers), np.arange(n_markers)]
        haps.append(np.stack([ha, hb]))
        return len(haps) - 1

    # 20 within-pop1 families: parents (pop1) -> two full sibs; 10 cross-population pairs
    # (pop1 x pop2 parents -> child, the pair KING misreads); 10 half-sib trios in adm;
    # 10 first-cousin pairs via pop2 sibs' children; one inbred child of full sibs
    fam_parents = []
    for k in range(20):
        pa, pb = 2 * k, 2 * k + 1                 # pop1 founders
        c1, c2 = child(pa, pb), child(pa, pb)
        labels += ["sib", "sib"]
        for c in (c1, c2):
            expected[(c, pa)] = 0.25
            expected[(c, pb)] = 0.25
        expected[(c2, c1)] = 0.25
        fam_parents.append((pa, pb, c1, c2))
    for k in range(10):
        pa, pb = 40 + k, n_pop + k                 # pop1 x pop2
        c = child(pa, pb)
        labels.append("cross_child")
        expected[(c, pa)] = 0.25
        expected[(c, pb)] = 0.25
    for k in range(10):
        pa, pb, pc = 2 * n_pop + 3 * k, 2 * n_pop + 3 * k + 1, 2 * n_pop + 3 * k + 2   # adm founders
        c1, c2 = child(pa, pb), child(pa, pc)
        labels += ["halfsib", "halfsib"]
        expected[(c2, c1)] = 0.125
        expected[(c1, pa)] = 0.25
        expected[(c2, pa)] = 0.25
    for k in range(10):
        pa, pb = n_pop + 20 + 2 * k, n_pop + 21 + 2 * k          # pop2 founders -> two sibs
        s1, s2 = child(pa, pb), child(pa, pb)
        labels += ["sib2", "sib2"]
        expected[(s2, s1)] = 0.25
        q1, q2 = n_pop + 60 + 2 * k, n_pop + 61 + 2 * k          # unrelated pop2 spouses
        k1, k2 = child(s1, q1), child(s2, q2)
        labels += ["cousin", "cousin"]
        expected[(k2, k1)] = 0.0625
        expected[(k1, s1)] = 0.25
        expected[(k1, s2)] = 0.125
    pa, pb, c1, c2 = fam_parents[0]
    inbred = child(c1, c2)
    labels.append("inbred")
    g = np.stack([h.sum(axis=0) for h in haps]).astype(np.int8)
    return g, labels, expected, inbred, nf


def pcs_from_unrelated(g, unrel, k=4):
    """PC-AiR-like: PCs from the unrelated founders, everyone projected on their loadings."""
    import numpy as np
    G = g.astype(np.float64)
    p = G[unrel].mean(axis=0) / 2
    sd = np.sqrt(2 * p * (1 - p))
    ok = sd > 0
    Z = (G[:, ok] - 2 * p[ok]) / sd[ok]
    U, S, Vt = np.linalg.svd(Z[unrel], full_matrices=False)
    scores = Z @ Vt[:k].T / S[:k]
    return scores


def selftest(a):
    set_threads(a.threads)
    import numpy as np
    rng = np.random.default_rng(20261031)
    g, labels, expected, inbred, nf = simulate(rng)
    N, M = g.shape
    unrel = np.zeros(N, dtype=bool)
    unrel[:nf] = True       # founders (relatives were appended after them)
    # the inbred child's parents are relatives of founders; founders themselves are unrelated
    pcs = pcs_from_unrelated(g, unrel, k=4)
    X = np.column_stack([np.ones(N), pcs])
    cutoff = 0.05
    fails = []

    def check(cond, msg):
        print(("ok   " if cond else "FAIL ") + msg)
        if not cond:
            fails.append(msg)

    def run_est(gm, Xd, train, n_markers=None, **kw):
        kw.setdefault("block_markers", 1000)
        kw.setdefault("block_samples", 97)
        return pcrelate_blocks(lambda s, c: gm[:, s:s + c].T, n_markers or M, Xd, train, cutoff,
                               log=lambda s: None, keep_all=True, **kw)

    pairs, f, n_used = run_est(g, X, unrel)
    kin = {k: v / 2 for k, v in pairs.items()}

    def key(i, j):
        return (i, j) if i > j else (j, i)

    # positive control: every planted pair recovered
    by_kind = {}
    for (i, j), e in expected.items():
        by_kind.setdefault(e, []).append(kin[key(i, j)] - e)
    for e, devs in sorted(by_kind.items()):
        devs = np.array(devs)
        check(abs(devs.mean()) < 0.01 and np.abs(devs).max() < 0.05,
              f"planted kinship {e}: n {len(devs)}, mean error {devs.mean():+.4f}, max |error| {np.abs(devs).max():.4f}")
    cross = [kin[key(i, j)] for (i, j), e in expected.items()
             if labels[i] == "cross_child" or labels[j] == "cross_child"]
    check(abs(np.mean(cross) - 0.25) < 0.01, f"cross-population parent-offspring at 0.25: mean {np.mean(cross):.4f}")
    # negative control: unrelated founders, within and across populations
    pop = np.array([{"pop1": 0, "pop2": 1}.get(l, 2) for l in labels])
    unrel_pairs = [(i, j) for i in range(nf) for j in range(i)]
    uk = np.array([kin[(i, j)] for i, j in unrel_pairs])
    across = np.array([pop[i] != pop[j] and pop[i] < 2 and pop[j] < 2 for i, j in unrel_pairs])
    within = np.array([pop[i] == pop[j] and pop[i] < 2 for i, j in unrel_pairs])
    check(abs(uk[across].mean()) < 0.005 and abs(uk[within].mean()) < 0.005,
          f"unrelated pairs at 0: across-population mean {uk[across].mean():+.4f}, within {uk[within].mean():+.4f}, sd {uk.std():.4f}")
    # the null stated first: the estimator's noise is Gaussian with the sd above, so at or above
    # the 0.125 cutoff (second degree, ~12 sd out here) no unrelated pair may appear. At 0.05 the
    # Gaussian expectation is ~0.03 pairs at 16,000 markers, but PC-Relate carries a second,
    # non-Gaussian term for ADMIXED pairs: the PC scores are estimated with noise eps_i, and
    # sum_m R_im R_jm picks up 4 eps_i eps_j sum_m beta_m^2, an individual-level product that
    # recurs across one sample's pairs and shrinks as 1/M (Conomos 2016 recommends many markers
    # for this reason). So: the pairs past 0.05 must all involve an admixed founder, and halving
    # the markers must not reduce their number (the direction the 1/M term predicts).
    from math import erfc, sqrt
    z = cutoff / (2 * uk.std())
    exp_above = len(uk) * 0.5 * erfc(z / sqrt(2))
    above = [(i, j) for (i, j), v in zip(unrel_pairs, uk) if 2 * v >= cutoff]
    n_above = len(above)
    adm = all(labels[i] == "adm" or labels[j] == "adm" for i, j in above)
    check(int((2 * uk >= 0.125).sum()) == 0,
          f"unrelated pairs at or above 0.125 (GRM scale, second degree): {int((2 * uk >= 0.125).sum())} of {len(uk)}")
    check(adm and n_above <= 10,
          f"unrelated pairs at or above the {cutoff} cutoff: {n_above} of {len(uk)} (Gaussian expectation "
          f"{exp_above:.3f}; every one involves an admixed founder: {adm}) -- the finite-marker PC-noise term")
    pairs_h, _, _ = run_est(g[:, : M // 2], X, unrel, n_markers=M // 2)
    n_above_h = sum(1 for (i, j) in unrel_pairs if pairs_h[(i, j)] >= cutoff)
    check(n_above_h >= n_above,
          f"the PC-noise term shrinks with markers: {n_above_h} unrelated pairs past the cutoff at {M // 2} markers, "
          f"{n_above} at {M}")
    # the PCs are load-bearing: intercept-only mu (a global allele frequency) must misread ancestry
    pairs0, _, _ = run_est(g, np.ones((N, 1)), unrel)
    uk0 = np.array([pairs0[(i, j)] / 2 for i, j in unrel_pairs])
    n_above0 = int((2 * uk0 >= cutoff).sum())
    check(uk0[across].mean() < -0.02 and uk0[within].mean() > 0.01 and n_above0 > 100,
          f"without PCs: across {uk0[across].mean():+.4f}, within {uk0[within].mean():+.4f}, {n_above0} unrelated pairs above the cutoff")
    # inbreeding: founders ~0, the child of full sibs ~0.25
    check(abs(np.median(f[:nf])) < 0.01 and abs(f[inbred] - 0.25) < 0.05,
          f"inbreeding: founders median {np.median(f[:nf]):+.4f}, child of full sibs {f[inbred]:.4f} (expected 0.25)")
    # two ways: one block against the small blocks above
    pairs1, f1, n1 = run_est(g, X, unrel, block_markers=M, block_samples=N)
    d = max(abs(pairs1[k] - pairs[k]) for k in pairs)
    check(d < 1e-4 and np.abs(f1 - f).max() < 1e-4 and n1 == n_used,
          f"blocking: one block vs 1000 x 97 blocks, max |diff| {d:.2e} (GRM scale), markers used {n_used}")
    # null: permute each marker across samples
    gp = g.copy()
    for m in range(M):
        gp[:, m] = gp[rng.permutation(N), m]
    pcs_p = pcs_from_unrelated(gp, unrel, k=4)
    pairs_p, _, _ = run_est(gp, np.column_stack([np.ones(N), pcs_p]), unrel)
    vp = np.array(list(pairs_p.values())) / 2     # kinship scale
    zp = cutoff / (2 * vp.std())
    exp_p = len(vp) * 0.5 * erfc(zp / sqrt(2))
    check(int((2 * vp >= cutoff).sum()) == 0 and abs(vp.mean()) < 0.002 and exp_p < 0.1,
          f"permuted genotypes: {int((2 * vp >= cutoff).sum())} of {len(vp)} pairs above the cutoff "
          f"(expected {exp_p:.3f}), kinship mean {vp.mean():+.5f}, sd {vp.std():.4f}")
    # missingness: 2 % missing at random changes nothing material
    gm = g.copy()
    gm[rng.uniform(size=g.shape) < 0.02] = -1
    pairs_m, _, _ = run_est(gm, X, unrel)
    dm = np.array([pairs_m[key(i, j)] / 2 - kin[key(i, j)] for (i, j) in expected])
    check(np.abs(dm).max() < 0.02, f"2 % missing genotypes: planted pairs move by at most {np.abs(dm).max():.4f}")
    # the .mtx round trip through run() on a .bed
    td_ctx = tempfile.TemporaryDirectory()
    td = a.keep or td_ctx.name
    os.makedirs(td, exist_ok=True)
    if True:
        pre = os.path.join(td, "sim")
        write_bed(pre, g, [f"S{i}" for i in range(N)])
        with open(pre + ".pcs.tsv", "w") as fh:
            fh.write("#FID\tIID\t" + "\t".join(f"PC{k + 1}" for k in range(4)) + "\n")
            for i in range(N):
                fh.write(f"0\tS{i}\t" + "\t".join(f"{v:.8g}" for v in pcs[i]) + "\n")
        with open(pre + ".unrel.txt", "w") as fh:
            fh.write("\n".join(f"S{i}" for i in range(nf)) + "\n")
        ns = argparse.Namespace(bed=pre, pgen=None, pcs=pre + ".pcs.tsv", npcs=None, unrelated=pre + ".unrel.txt",
                                out=pre + "_pcrel", cutoff=cutoff, maf=0.01, diag="inbreeding",
                                block_markers=1000, block_samples=97, threads=a.threads, mask=1)
        run(ns)
        n_out, ent = read_mtx(pre + "_pcrel.mtx")
        off = {k: v for k, v in ent.items() if k[0] != k[1]}
        want = {k: v for k, v in pairs.items() if v >= cutoff}
        same = set(off) == set(want) and all(abs(off[k] - want[k]) < 1e-4 for k in want)
        check(n_out == N and same and abs(ent[(inbred, inbred)] - 1.25) < 0.05,
              f".bed round trip: {len(off)} off-diagonal pairs kept, identical to the in-memory run; "
              f"inbred diagonal {ent[(inbred, inbred)]:.3f}")
        ids = read_ids(pre + "_pcrel.mtx.sampleIDs.txt")
        check(ids == [f"S{i}" for i in range(N)], "sample IDs in .fam order")
        # the PGEN arm: the same cohort through plink2 --make-pgen and pgenlib must give the same
        # pairs (plink2 may swap REF/ALT, which cannot move a kinship: g -> 2 - g flips R's sign for
        # every sample and leaves A alone)
        import shutil, subprocess
        plink2 = shutil.which("plink2")
        try:
            import pgenlib  # noqa: F401
            have_pgenlib = True
        except ImportError:
            have_pgenlib = False
        if plink2 and have_pgenlib:
            subprocess.run([plink2, "--bfile", pre, "--make-pgen", "--out", pre + "_p"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            ns.bed, ns.pgen, ns.out = None, pre + "_p", pre + "_pcrel_pgen"
            run(ns)
            n_p, ent_p = read_mtx(pre + "_pcrel_pgen.mtx")
            dp = max(abs(ent_p[k] - ent[k]) for k in ent) if set(ent_p) == set(ent) else float("inf")
            check(n_p == N and dp < 1e-5 and read_ids(pre + "_pcrel_pgen.mtx.sampleIDs.txt") == ids,
                  f"PGEN arm: plink2 --make-pgen + pgenlib reproduce the .bed run, max |diff| {dp:.2e}")
        else:
            print(f"skip PGEN arm: plink2 {'found' if plink2 else 'not on PATH'}, "
                  f"pgenlib {'importable' if have_pgenlib else 'not importable'}")
        if a.keep:
            # every pair, for an external reconciliation (GENESIS pcrelate on the same bed and PCs)
            with open(pre + ".pairs_all.tsv", "w") as fh:
                fh.write("i\tj\tkinship\n")
                for (i, j), v in sorted(pairs.items()):
                    fh.write(f"S{i}\tS{j}\t{v / 2:.8g}\n")
            with open(pre + ".labels.tsv", "w") as fh:
                fh.write("iid\tlabel\tinbreeding\n")
                for i in range(N):
                    fh.write(f"S{i}\t{labels[i]}\t{f[i]:.6g}\n")
            print(f"fixture kept in {td}")
    td_ctx.cleanup()
    if fails:
        raise SystemExit(f"selftest: {len(fails)} FAILED")
    print(f"selftest OK: {N} samples, {M} markers, {len(expected)} planted pairs")


def write_bed(prefix, g, ids):
    """int8 (N x M) with -1 missing -> PLINK 1 SNP-major .bed/.bim/.fam"""
    import numpy as np
    N, M = g.shape
    code = np.array([0, 2, 3], dtype=np.uint8)   # 0 copies -> 00, 1 -> 10, 2 -> 11; missing -> 01
    gc = np.where(g < 0, 1, code[np.clip(g, 0, 2)]).astype(np.uint8).T   # M x N codes
    bpm = (N + 3) // 4
    pad = np.zeros((M, bpm * 4), dtype=np.uint8)
    pad[:, :N] = gc
    pad = pad.reshape(M, bpm, 4)
    byts = (pad[:, :, 0] | (pad[:, :, 1] << 2) | (pad[:, :, 2] << 4) | (pad[:, :, 3] << 6)).astype(np.uint8)
    with open(prefix + ".bed", "wb") as fh:
        fh.write(b"\x6c\x1b\x01")
        fh.write(byts.tobytes())
    with open(prefix + ".bim", "w") as fh:
        for m in range(M):
            fh.write(f"1\tm{m}\t0\t{m + 1}\tA\tC\n")
    with open(prefix + ".fam", "w") as fh:
        for s in ids:
            fh.write(f"0\t{s}\t0\t0\t0\t-9\n")


# ----------------------------------------------------------------------------- CLI

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    gi = r.add_mutually_exclusive_group(required=True)
    gi.add_argument("--bed", help="PLINK 1 prefix (.bed/.bim/.fam), LD-pruned markers")
    gi.add_argument("--pgen", help="PLINK 2 prefix (.pgen/.pvar/.psam), hard calls, LD-pruned markers")
    r.add_argument("--pcs", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--npcs", type=int, default=None)
    r.add_argument("--unrelated", default=None)
    r.add_argument("--cutoff", type=float, default=0.05, help="GRM scale (2 x kinship) [%(default)s]")
    r.add_argument("--maf", type=float, default=0.01)
    r.add_argument("--diag", choices=["one", "inbreeding"], default="one")
    r.add_argument("--block-markers", type=int, default=2048)
    r.add_argument("--block-samples", type=int, default=None,
                   help="sample tile; default: the largest that keeps the two accumulators within --memory-gb")
    r.add_argument("--memory-gb", type=float, default=16.0,
                   help="budget for the two float32 accumulators (2 x tile x N x 4 B) [%(default)s]")
    r.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    r.add_argument("--mask", type=int, default=20)
    s = sub.add_parser("selftest")
    s.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    s.add_argument("--keep", default=None, help="keep the simulated .bed/.pcs/.unrel and every pair here")
    a = ap.parse_args()
    if a.cmd == "run":
        run(a)
    else:
        selftest(a)


if __name__ == "__main__":
    main()
