#!/usr/bin/env python3
"""A sparse GRM for SAIGE from a KING-robust kinship table, and its comparison with an
existing SAIGE .mtx (LEDGER #183, 2026-10-28).

Why: amr's sparse GRM from createSparseGRM (5,000 random markers, standardised genotypes)
carries the admixture background in its off-diagonals -- 46M pairs above 0.05 -- which forced
the 0.125 cutoff, and a GRM thinned that hard is not positive definite (aou_sigma_pd_check.R).
KING-robust estimates kinship from within-pair genotype counts, so ancestry does not leak into
the off-diagonals; SAIGE accepts any sparse kinship in Matrix Market form, so this is an INPUT
change: no code path moves.

  convert   plink2 --make-king-table output (.kin0: #FID1 ID1 FID2 ID2 NSNP HETHET IBS0 KINSHIP,
            or the IID1/IID2 spelling) -> <out>.mtx on SAIGE's scale (off-diagonal = 2 x kinship,
            diagonal = 1, symmetric coordinate, lower triangle, 1-based) + <out>.mtx.sampleIDs.txt
            in the order of --samples (one IID per line, e.g. the .psam/.fam IIDs), keeping pairs
            at or above --cutoff on the GRM scale. An ID in the table but not in --samples is an
            error (the table was computed on these samples or it was not).
  compare   the KING-derived .mtx against an existing SAIGE .mtx on the SAME samples: (1) close
            relatives agree (pairs the existing GRM puts at or above 0.4: the new value's median
            ratio and spread -- the positive control; read it first); (2) pairs at or above each
            cutoff under each estimator (the admixture background is the difference); (3) the
            0.05..0.125 band: pairs both keep (real third-degree) against pairs only the old one
            keeps (background); (4) the diagonal of the old one. Counts of pairs below --mask are
            written "<mask" (real data); no sample ID is ever printed.
  plot      three panels to <out>.png: (a) every pair of the old matrix against KING's value
            (absent from the KING table = 0) as a 2-D histogram, (b) the pairs both keep against
            the identity line, (c) the old values split by whether KING keeps the pair. Cells
            holding fewer than --mask pairs are not drawn, so no small count is readable.
  selftest  fabricated data: three families planted, a background of small GRM-only entries,
            every claim checked (compare and plot).

Usage:
  king_to_mtx.py convert --king king.kin0 --samples ids.txt --out amr_king [--cutoff 0.05]
  king_to_mtx.py compare --new amr_king.mtx --old old.sparseGRM.mtx --report cmp.tsv [--mask 20]
  king_to_mtx.py plot --new amr_king.mtx --old old.sparseGRM.mtx --out king_vs_grm [--mask 20]
  king_to_mtx.py selftest
The old matrix at 0.05 has ~46M entries on amr: compare and plot read it with pandas into flat
arrays (a few hundred MB), never a per-entry Python object.
"""
import argparse, gzip, io, os, sys, tempfile
import numpy as np


def openr(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_ids(p):
    with openr(p) as fh:
        ids = [l.strip().split()[0] for l in fh if l.strip()]
    if len(set(ids)) != len(ids):
        raise SystemExit("--samples has a repeated ID (SAIGE keys on IID alone)")
    return ids


def read_king(p):
    """(id1, id2, kinship) per row of a plink2 --make-king-table file."""
    with openr(p) as fh:
        head = fh.readline().lstrip("#").split()
        cols = {c.upper(): i for i, c in enumerate(head)}
        i1 = cols.get("ID1", cols.get("IID1")); i2 = cols.get("ID2", cols.get("IID2")); ik = cols.get("KINSHIP")
        if None in (i1, i2, ik):
            raise SystemExit(f"KING table header not recognised: {head}")
        for l in fh:
            f = l.split()
            if not f:
                continue
            yield f[i1], f[i2], float(f[ik])


def read_mtx(p):
    """A SAIGE .mtx: dict (i, j) with i >= j (0-based) -> value, plus n."""
    ent = {}
    with openr(p) as fh:
        head = fh.readline()
        if "symmetric" not in head:
            raise SystemExit("expected a symmetric coordinate .mtx")
        n = int(fh.readline().split()[0])
        for l in fh:
            f = l.split()
            if not f:
                continue
            i, j, x = int(f[0]) - 1, int(f[1]) - 1, float(f[2])
            if i < j:
                i, j = j, i
            ent[(i, j)] = x
    return n, ent


def read_mtx_frame(p):
    """A SAIGE .mtx as flat arrays: n, and a DataFrame (i, j, x) with i >= j, 0-based."""
    import pandas as pd
    with openr(p) as fh:
        head = fh.readline()
        if "symmetric" not in head:
            raise SystemExit("expected a symmetric coordinate .mtx")
        n = int(fh.readline().split()[0])
    df = pd.read_csv(p, sep=r"\s+", skiprows=2, header=None, names=["i", "j", "x"],
                     dtype={"i": np.int64, "j": np.int64, "x": np.float64}, comment="%")
    df["i"] -= 1; df["j"] -= 1
    lo = df["i"] < df["j"]
    df.loc[lo, ["i", "j"]] = df.loc[lo, ["j", "i"]].values
    return n, df


def pair_tables(new_path, old_path):
    """Off-diagonal pairs of both matrices on the old matrix's indices, joined on IID:
    (old off-diagonal frame, new off-diagonal frame re-indexed, old diagonal values, counts)."""
    import pandas as pd
    nN, new = read_mtx_frame(new_path); nO, old = read_mtx_frame(old_path)
    idsN = read_ids(new_path + ".sampleIDs.txt"); idsO = read_ids(old_path + ".sampleIDs.txt")
    if len(idsN) != nN or len(idsO) != nO:
        raise SystemExit("an .mtx and its sampleIDs file disagree on the sample count")
    posO = {s: k for k, s in enumerate(idsO)}
    remap = np.array([posO.get(s, -1) for s in idsN], dtype=np.int64)
    offN = new[new["i"] != new["j"]].copy()
    offN["i"] = remap[offN["i"].values]; offN["j"] = remap[offN["j"].values]
    offN = offN[(offN["i"] >= 0) & (offN["j"] >= 0)]
    lo = offN["i"] < offN["j"]
    offN.loc[lo, ["i", "j"]] = offN.loc[lo, ["j", "i"]].values
    diagO = old[old["i"] == old["j"]]["x"].values
    offO = old[old["i"] != old["j"]]
    common = int((remap >= 0).sum())
    return offO, offN, diagO, dict(samples_new=nN, samples_old=nO, samples_common=common)


def write_mtx(path, n, ent):
    with open(path, "w") as fh:
        fh.write("%%MatrixMarket matrix coordinate real symmetric\n")
        fh.write(f"{n} {n} {len(ent)}\n")
        for (i, j) in sorted(ent):
            fh.write(f"{i + 1} {j + 1} {ent[(i, j)]!r}\n")


def convert(a):
    ids = read_ids(a.samples); pos = {s: k for k, s in enumerate(ids)}; n = len(ids)
    ent = {(k, k): 1.0 for k in range(n)}
    kept = dropped = unknown = 0
    for id1, id2, kin in read_king(a.king):
        if id1 not in pos or id2 not in pos:
            unknown += 1
            continue
        v = 2.0 * kin
        if v < a.cutoff:
            dropped += 1
            continue
        i, j = pos[id1], pos[id2]
        if i == j:
            continue
        if i < j:
            i, j = j, i
        ent[(i, j)] = v
        kept += 1
    if unknown:
        raise SystemExit(f"{unknown} KING pair(s) name a sample absent from --samples")
    write_mtx(a.out + ".mtx", n, ent)
    with open(a.out + ".mtx.sampleIDs.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    print(f"samples {n}; pairs kept at >= {a.cutoff} (GRM scale): {kept}; below: {dropped}")
    print(f"wrote {a.out}.mtx and {a.out}.mtx.sampleIDs.txt")


def mask(k, m):
    return f"<{m}" if 0 < k < m else str(k)


def compare(a):
    import pandas as pd
    offO, offN, diagO, cnt = pair_tables(a.new, a.old)
    rows = []
    def note(k, v):
        rows.append((k, str(v))); print(f"{k:48s} {v}")
    for k in ("samples_new", "samples_old", "samples_common"):
        note(k, cnt[k])
    m = offO.merge(offN, on=["i", "j"], how="left", suffixes=("_old", "_new"))
    m["x_new"] = m["x_new"].fillna(0.0)
    # (1) close relatives agree: the positive control
    close = m[m["x_old"] >= 0.4]
    if len(close):
        present = close[close["x_new"] > 0]
        r = (present["x_new"] / present["x_old"]).values
        note("old_pairs_at_or_above_0.4", mask(len(close), a.mask))
        note("  absent_from_new", mask(len(close) - len(present), a.mask))
        if len(r):
            note("  new_over_old_median", f"{np.median(r):.3f}")
            note("  new_over_old_5th_95th", f"{np.percentile(r, 5):.3f} .. {np.percentile(r, 95):.3f}")
    else:
        note("old_pairs_at_or_above_0.4", "0 (no first-degree pair in the old matrix on the common samples)")
    # (2) how many pairs each keeps at each cutoff
    for c in (0.05, 0.125, 0.25, 0.4):
        note(f"pairs_at_or_above_{c}_old", mask(int((offO["x"] >= c).sum()), a.mask))
        note(f"pairs_at_or_above_{c}_new", mask(int((offN["x"] >= c).sum()), a.mask))
    # (3) the 0.05..0.125 band: real third-degree against background
    band = m[(m["x_old"] >= 0.05) & (m["x_old"] < 0.125)]
    both = int((band["x_new"] >= 0.05).sum())
    note("band_0.05_0.125_old_pairs", mask(len(band), a.mask))
    note("  of which_new_also_at_or_above_0.05", mask(both, a.mask))
    note("  of which_new_below_0.05_(background)", mask(len(band) - both, a.mask))
    # (4) the old diagonal
    if len(diagO):
        note("old_diag_mean_min_max", f"{diagO.mean():.4f} {diagO.min():.4f} {diagO.max():.4f}")
    note("old_diag_entries_missing", mask(cnt["samples_old"] - len(diagO), a.mask))
    with open(a.report, "w") as fh:
        fh.write("key\tvalue\n")
        for k, v in rows:
            fh.write(f"{k.strip()}\t{v}\n")
    print("done")


def plot(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    offO, offN, diagO, cnt = pair_tables(a.new, a.old)
    m = offO.merge(offN, on=["i", "j"], how="outer", suffixes=("_old", "_new")).fillna(0.0)
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.5))
    top = max(float(m["x_old"].max()), float(m["x_new"].max()), 0.6)
    # (a) every pair either keeps: cells below the mask are not drawn
    h = ax[0].hexbin(m["x_old"], m["x_new"], gridsize=60, bins="log", mincnt=a.mask, extent=(0, top, 0, top))
    ax[0].plot([0, top], [0, top], color="0.3", lw=0.8)
    ax[0].set_title(f"every pair either matrix keeps (cells below {a.mask} pairs not drawn)", fontsize=10)
    fig.colorbar(h, ax=ax[0], label="pairs (log10)")
    # (b) pairs both keep
    both = m[(m["x_old"] > 0) & (m["x_new"] > 0)]
    if len(both) >= a.mask:
        h2 = ax[1].hexbin(both["x_old"], both["x_new"], gridsize=60, bins="log", mincnt=a.mask, extent=(0, top, 0, top))
        fig.colorbar(h2, ax=ax[1], label="pairs (log10)")
    ax[1].plot([0, top], [0, top], color="0.3", lw=0.8)
    for c in (0.05, 0.125, 0.25, 0.5):
        ax[1].axvline(c, color="r", lw=0.4, ls=":"); ax[1].axhline(c, color="r", lw=0.4, ls=":")
    ax[1].set_title(f"pairs both keep: {mask(len(both), a.mask)}", fontsize=10)
    # (c) the old values, split by whether KING keeps the pair
    old_only = m[(m["x_old"] > 0) & (m["x_new"] < 0.05)]["x_old"]
    kept = m[(m["x_old"] > 0) & (m["x_new"] >= 0.05)]["x_old"]
    edges = np.linspace(0, top, 61)
    for vals, lab, col in ((old_only, "KING below 0.05 (background)", "#7f7f7f"), (kept, "KING at or above 0.05", "#d62728")):
        cnts, _ = np.histogram(vals, bins=edges)
        cnts = np.where(cnts < a.mask, 0, cnts)   # bins below the mask are not drawn
        ax[2].step(edges[:-1], cnts, where="post", color=col, label=f"{lab}: {mask(len(vals), a.mask)} pairs")
    ax[2].set_yscale("log"); ax[2].legend(fontsize=8); ax[2].set_title("old GRM values by KING's verdict", fontsize=10)
    for x in ax:
        x.set_xlabel("step-0 GRM (2 x kinship)"); x.grid(alpha=0.3)
    ax[0].set_ylabel("KING-robust (2 x kinship)"); ax[1].set_ylabel("KING-robust (2 x kinship)"); ax[2].set_ylabel("pairs")
    fig.suptitle(f"KING-robust against the step-0 GRM: {cnt['samples_common']:,} shared samples", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(a.out + ".png", dpi=110); plt.close(fig)
    print("wrote", a.out + ".png")


def selftest(_):
    rng = np.random.default_rng(5); n = 300
    ids = [f"s{k}" for k in range(n)]
    T = tempfile.mkdtemp()
    with open(f"{T}/ids.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    # truth: three families (a trio-ish star at 0.5, a sib pair at 0.5, two cousins at 0.125)
    fam = {(1, 0): 0.5, (2, 0): 0.5, (11, 10): 0.5, (21, 20): 0.125}
    # the KING table: true pairs plus noise, kinship = GRM/2
    with open(f"{T}/king.kin0", "w") as fh:
        fh.write("#FID1\tID1\tFID2\tID2\tNSNP\tHETHET\tIBS0\tKINSHIP\n")
        for (i, j), v in fam.items():
            fh.write(f"0\t{ids[i]}\t0\t{ids[j]}\t100000\t0\t0\t{v / 2 + rng.normal(0, 0.002):.5f}\n")
        for _ in range(50):   # unrelated pairs KING reports near zero (some slightly negative)
            i, j = rng.choice(n, 2, replace=False)
            fh.write(f"0\t{ids[i]}\t0\t{ids[j]}\t100000\t0\t0\t{rng.normal(0, 0.005):.5f}\n")
    # the old GRM: the families plus an admixture background of 400 entries at 0.05..0.12
    old = {(k, k): 1.0 + rng.normal(0, 0.01) for k in range(n)}
    for (i, j), v in fam.items():
        old[(i, j)] = v + rng.normal(0, 0.01)
    bg = 0
    while bg < 400:
        i, j = sorted(rng.choice(n, 2, replace=False), reverse=True)
        if (i, j) not in old:
            old[(i, j)] = rng.uniform(0.05, 0.12); bg += 1
    write_mtx(f"{T}/old.mtx", n, old)
    with open(f"{T}/old.mtx.sampleIDs.txt", "w") as fh:
        fh.write("\n".join(ids) + "\n")
    convert(argparse.Namespace(king=f"{T}/king.kin0", samples=f"{T}/ids.txt", out=f"{T}/new", cutoff=0.05))
    nN, new = read_mtx(f"{T}/new.mtx")
    off = {k: x for k, x in new.items() if k[0] != k[1]}
    assert nN == n and len(off) == 4, f"expected the 4 planted pairs at >= 0.05, got {len(off)}"
    assert all(abs(off[k] - v) < 0.02 for k, v in fam.items()), "planted pairs not on the GRM scale"
    assert all(new[(k, k)] == 1.0 for k in range(n)), "diagonal must be 1"
    buf = io.StringIO(); so = sys.stdout; sys.stdout = buf
    try:
        compare(argparse.Namespace(new=f"{T}/new.mtx", old=f"{T}/old.mtx", report=f"{T}/cmp.tsv", mask=0))
    finally:
        sys.stdout = so
    rep = dict(l.rstrip("\n").split("\t") for l in open(f"{T}/cmp.tsv") if not l.startswith("key"))
    assert rep["old_pairs_at_or_above_0.4"] == "3", rep
    assert rep["absent_from_new"] == "0"
    assert 0.95 < float(rep["new_over_old_median"]) < 1.05, rep["new_over_old_median"]
    assert rep["pairs_at_or_above_0.05_old"] == str(400 + 4) and rep["pairs_at_or_above_0.05_new"] == "4"
    assert rep["pairs_at_or_above_0.25_old"] == "3" and rep["pairs_at_or_above_0.25_new"] == "3"
    assert rep["band_0.05_0.125_old_pairs"] == "400", rep["band_0.05_0.125_old_pairs"]
    assert rep["of which_new_also_at_or_above_0.05"] == "0" and rep["of which_new_below_0.05_(background)"] == "400"
    assert rep["old_diag_entries_missing"] == "0"
    # masking: a count below the mask is written "<mask"
    buf = io.StringIO(); sys.stdout = buf
    try:
        compare(argparse.Namespace(new=f"{T}/new.mtx", old=f"{T}/old.mtx", report=f"{T}/cmp2.tsv", mask=20))
    finally:
        sys.stdout = so
    rep2 = dict(l.rstrip("\n").split("\t") for l in open(f"{T}/cmp2.tsv") if not l.startswith("key"))
    assert rep2["old_pairs_at_or_above_0.4"] == "<20" and rep2["pairs_at_or_above_0.05_old"] == "404", rep2
    plot(argparse.Namespace(new=f"{T}/new.mtx", old=f"{T}/old.mtx", out=f"{T}/king_vs_grm", mask=20))
    assert os.path.getsize(f"{T}/king_vs_grm.png") > 10000
    print("selftest OK: convert (4 planted pairs on the GRM scale, diagonal 1), compare (positive control "
          "ratio ~1, 400 background pairs only in the old matrix, the band split, masking), plot (PNG written)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("convert"); c.add_argument("--king", required=True); c.add_argument("--samples", required=True)
    c.add_argument("--out", required=True); c.add_argument("--cutoff", type=float, default=0.05); c.set_defaults(fn=convert)
    k = sub.add_parser("compare"); k.add_argument("--new", required=True); k.add_argument("--old", required=True)
    k.add_argument("--report", required=True); k.add_argument("--mask", type=int, default=20); k.set_defaults(fn=compare)
    g = sub.add_parser("plot"); g.add_argument("--new", required=True); g.add_argument("--old", required=True)
    g.add_argument("--out", required=True); g.add_argument("--mask", type=int, default=20); g.set_defaults(fn=plot)
    s = sub.add_parser("selftest"); s.set_defaults(fn=selftest)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
