#!/usr/bin/env python3
"""AlphaMissense rows for paralog copies AlphaMissense does not score (2026-10-22).

AlphaMissense scores each UniProt protein on ONE transcript. Where several genes encode the identical protein (the
histone clusters H3C*/H4C*/H2AC*/H2BC*, CT45A/CT47A, SPDYE, USP17L, GAGE, AMY1A, ...) the other copies have no row at
all, so their missense variants got no score. For an identical protein the same amino-acid change is the same
mutation to AlphaMissense, so this writes, for every protein-coding transcript of a gene with NO AlphaMissense row
on any of its transcripts whose protein is identical to a protein AlphaMissense scores in its canonical file, every
possible missense SNV of that transcript with AlphaMissense's score for that protein change
(AlphaMissense_aa_substitutions.tsv.gz, keyed by UniProt). Output: AlphaMissense isoforms-file layout, one file per
chromosome, read by flexrv_score_from_alphamissense.py --extra-isoforms; a transcript is never a near match (one
residue different: not used).

Controls: (1) the generator on a random sample of transcripts AlphaMissense DOES score must reproduce AlphaMissense's
own canonical rows exactly for >= 95% (2026-10-22: 288 of 297; the rest are UniProt sequences that differ from
GENCODE's protein); (2) a copy is borrowed from only if the generator reproduces ITS AlphaMissense rows exactly -- the (pos, ref, alt) set and every score -- which tests the
coordinate mapping, the translation and the protein-level lookup together.

usage: alphamissense_paralogs.py --gtf gencode.v39.annotation.gtf.gz --translations gencode.v39.pc_translations.fa.gz \
         --transcripts gencode.v39.pc_transcripts.fa.gz --am AlphaMissense_hg38.tsv.gz \
         --isoforms AlphaMissense_isoforms_hg38.tsv.gz --aa AlphaMissense_aa_substitutions.tsv.gz \
         --out-dir resources/am_by_chrom [--check-n 200]
"""
import argparse, gzip, os, random, re, sys
from collections import defaultdict

BASES = "TCAG"
CODE = dict(zip((a + b + c for a in BASES for b in BASES for c in BASES),
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path, key):
    out, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                if name: out[name[0]] = (name[1], "".join(buf))
                h = l[1:].rstrip("\n").split("|"); name = key(h); buf = []
            else:
                buf.append(l.strip())
    if name: out[name[0]] = (name[1], "".join(buf))
    return out


def read_gtf(path):
    """{transcript (versioned): (chrom, strand, gene (stripped), [(start, end)] CDS in transcript order)}"""
    tx = {}
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.split("\t", 8)
            if f[2] != "CDS": continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            if t.endswith("_PAR_Y"): continue
            g = re.search(r'gene_id "([^".]+)', f[8]).group(1)
            e = tx.setdefault(t, [f[0], f[6], g, []]); e[3].append((int(f[3]), int(f[4])))
    for t, e in tx.items():
        e[3].sort(reverse=(e[1] == "-"))
    return tx


def snvs(tx_full, gtf, translations, transcripts):
    """Every missense SNV of a transcript: [(chrom, pos, ref, alt, protein_variant)], or None when the CDS does not
    translate to the GENCODE protein (incomplete CDS, selenocysteine, ...)."""
    if tx_full not in gtf or tx_full not in translations or tx_full not in transcripts: return None
    chrom, strand, _, segs = gtf[tx_full]
    prot = translations[tx_full][1]; meta, seq = transcripts[tx_full]
    m = re.search(r"CDS:(\d+)-(\d+)", meta)
    if not m: return None
    cds = seq[int(m.group(1)) - 1:int(m.group(1)) - 1 + 3 * len(prot)].upper()
    gpos = [p for s, e in segs for p in (range(s, e + 1) if strand == "+" else range(e, s - 1, -1))]
    if len(cds) != 3 * len(prot) or len(gpos) < len(cds):
        return None
    gpos = gpos[:len(cds)]
    if "".join(CODE.get(cds[i:i + 3], "X") for i in range(0, len(cds), 3)) != prot:
        return None
    out = []
    for i, b in enumerate(cds):
        ci = i - i % 3
        if ci == 0: continue   # the start codon: AlphaMissense's genomic files have no row there (start lost)
        codon = cds[ci:ci + 3]; ref_aa = CODE[codon]
        for a in "ACGT":
            if a == b: continue
            alt_aa = CODE[codon[:i % 3] + a + codon[i % 3 + 1:]]
            if alt_aa in (ref_aa, "*"): continue
            ref_g, alt_g = (b, a) if strand == "+" else (b.translate(COMP), a.translate(COMP))
            out.append((chrom, gpos[i], ref_g, alt_g, "%s%d%s" % (ref_aa, ci // 3 + 1, alt_aa)))
    return out


def main(a):
    gtf = read_gtf(a.gtf)
    translations = read_fasta(a.translations, lambda h: (h[1], h[0]))
    translations = {t: (p, s.rstrip("*")) for t, (p, s) in translations.items()}
    transcripts = read_fasta(a.transcripts, lambda h: (h[0], "|".join(h)))
    stripped = {t.split(".")[0]: t for t in gtf}
    # transcripts AlphaMissense scores; canonical ones with their UniProt accession
    scored, canon_up = set(), {}
    for path, it, iu in ((a.am, 6, 5), (a.isoforms, 5, None)):
        seen = None
        with gzip.open(path, "rt") as fh:
            for l in fh:
                if l[0] == "#": continue
                f = l.split("\t", 7); t = f[it].replace("_PAR_Y", "").split(".")[0]
                if t != seen:
                    scored.add(t); seen = t
                    if iu is not None: canon_up.setdefault(t, f[iu].split("-")[0])
    print("# AlphaMissense transcripts: %d (%d canonical)" % (len(scored), len(canon_up)), file=sys.stderr)
    by_seq = defaultdict(set)   # protein sequence -> (UniProt, transcript) of AlphaMissense canonical transcripts encoding it
    for t, up in canon_up.items():
        full = stripped.get(t)
        if full in translations: by_seq[translations[full][1]].add((up, full))
    genes_scored = {gtf[stripped[t]][2] for t in scored if t in stripped}
    targets = {}   # transcript -> UniProt
    for t, (chrom, strand, g, segs) in gtf.items():
        if g in genes_scored or t not in translations: continue
        src = by_seq.get(translations[t][1])
        if src: targets[t] = sorted(src)
    print("# target transcripts (genes with no AlphaMissense row, protein identical to a scored one): %d in %d genes"
          % (len(targets), len({gtf[t][2] for t in targets})), file=sys.stderr)
    # control: regenerate a random sample of AlphaMissense canonical transcripts
    rng = random.Random(20261022)
    ctl = rng.sample(sorted(t for t in canon_up if stripped.get(t) in translations), a.check_n)
    need = {up for src in targets.values() for up, _ in src} | {canon_up[t] for t in ctl}
    aa = defaultdict(dict)
    with gzip.open(a.aa, "rt") as fh:
        for l in fh:
            if l[0] == "#" or l.startswith("uniprot_id"): continue
            f = l.rstrip("\n").split("\t")
            if f[0] in need: aa[f[0]][f[1]] = (f[2], f[3])
    ctl_tx = {stripped[t]: canon_up[t] for t in ctl}
    src_tx = {t: up for src in targets.values() for up, t in src}
    want = defaultdict(dict)
    with gzip.open(a.am, "rt") as fh:
        ctl_s = {t.split(".")[0] for t in list(ctl_tx) + list(src_tx)}
        for l in fh:
            if l[0] == "#": continue
            f = l.rstrip("\n").split("\t")
            if f[6].split(".")[0] in ctl_s:
                want[f[6].split(".")[0]][(f[0], int(f[1]), f[2], f[3])] = f[8]
    def reproduces(t, up):
        """The generator on an AlphaMissense canonical transcript gives exactly AlphaMissense's rows."""
        rows = snvs(t, gtf, translations, transcripts)
        if rows is None: return None, 0, 0, 0, 0
        got = {(c, p, r, al): aa[up].get(pv, ("NA",))[0] for c, p, r, al, pv in rows}
        w = want[t.split(".")[0]]; sh = set(got) & set(w)
        return set(got) == set(w) and all(got[k] == w[k] for k in w), len(got), len(w), len(sh), sum(got[k] == w[k] for k in sh)
    ok = bad = untranslatable = 0
    for t, up in ctl_tx.items():
        res, ng, nw, nsh, neq = reproduces(t, up)
        if res is None: untranslatable += 1; continue
        if res: ok += 1
        else:
            bad += 1
            if bad <= 5:
                print("# control mismatch %s: generated %d rows, AlphaMissense %d, shared %d, equal scores %d" % (
                    t, ng, nw, nsh, neq), file=sys.stderr)
    print("# control (AlphaMissense canonical transcripts regenerated): %d exact, %d differ, %d not translatable"
          % (ok, bad, untranslatable), file=sys.stderr)
    if a.check_only: return 0 if bad == 0 else 1
    if bad > 0.05 * max(1, ok + bad):
        sys.exit("control failed: more than 5% of regenerated transcripts differ from AlphaMissense; nothing written")
    # per target: borrow only from a scored copy the generator reproduces EXACTLY (its UniProt sequence lines up with
    # the GENCODE protein); the few that do not (UniProt / GENCODE sequence offsets) are never used
    good_src = {t: reproduces(t, up)[0] is True for t, up in src_tx.items()}
    print("# source copies reproduced exactly: %d of %d" % (sum(good_src.values()), len(good_src)), file=sys.stderr)
    by_chrom = defaultdict(list); report = []
    for t, src in sorted(targets.items()):
        ok_src = [(u, st) for u, st in src if good_src.get(st)]
        if not ok_src:
            report.append((t, gtf[t][2], ",".join(u for u, _ in src), "source not reproduced", 0)); continue
        up, st = ok_src[0]
        rows = snvs(t, gtf, translations, transcripts)
        if rows is None: report.append((t, gtf[t][2], up, "untranslatable", 0)); continue
        n = 0
        for c, p, r, al, pv in rows:
            if pv in aa[up]:
                s, cl = aa[up][pv]; by_chrom[c].append((p, r, al, t, pv, s, cl)); n += 1
        report.append((t, gtf[t][2], up, "ok from " + st, n))
    os.makedirs(a.out_dir, exist_ok=True)
    for c, rows in by_chrom.items():
        with gzip.open(os.path.join(a.out_dir, "AlphaMissense_paralogs_hg38.%s.tsv.gz" % c), "wt") as fh:
            fh.write("# AlphaMissense protein-level scores (aa_substitutions) placed on paralog copies AlphaMissense "
                     "does not score; alphamissense_paralogs.py\n")
            fh.write("#CHROM\tPOS\tREF\tALT\tgenome\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n")
            for p, r, al, t, pv, s, cl in sorted(rows):
                fh.write("%s\t%d\t%s\t%s\thg38\t%s\t%s\t%s\t%s\n" % (c, p, r, al, t, pv, s, cl))
    with open(os.path.join(a.out_dir, "AlphaMissense_paralogs_report.tsv"), "w") as fh:
        fh.write("transcript\tgene\tuniprot\tstatus\trows\n")
        for r in report: fh.write("%s\t%s\t%s\t%s\t%d\n" % r)
    print("# wrote %d rows for %d transcripts on %d chromosomes" % (
        sum(len(v) for v in by_chrom.values()), sum(1 for r in report if r[3].startswith("ok")), len(by_chrom)), file=sys.stderr)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for k in ("gtf", "translations", "transcripts", "am", "isoforms", "aa", "out-dir"):
        ap.add_argument("--" + k, required=(k != "out-dir"))
    ap.add_argument("--check-n", type=int, default=300)
    ap.add_argument("--check-only", action="store_true")
    a = ap.parse_args()
    a.out_dir = a.out_dir or "resources/am_by_chrom"
    sys.exit(main(a))
