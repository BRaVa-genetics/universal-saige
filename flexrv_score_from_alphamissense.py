#!/usr/bin/env python3
"""Add an AlphaMissense `score:NAME` line to a SAIGE group file, for FlexRV.

FlexRV reads ONE per-marker pathogenicity score in [0, 1] from a group file
line `GENE score:NAME v1 v2 ...`, placed IMMEDIATELY AFTER the `anno` line and
BEFORE any `weight` line (the position is fixed, not merely conventional:
`R/SAIGE_SPATest_Region_Func.R` reads region 1's third line to decide whether
the file carries a score at all). It is NOT a `weight` line -- `--weights.beta`
and `weight:NAME` sets are inert on FlexRV's transform sets.

`checkGroupFile()` enforces four things this builder must satisfy, so they are
also what its self-test asserts:
  * every region carries a score line, named identically across regions;
  * its length equals that region's `var` line length -- EVERY variant needs a
    value, including the non-coding and synonymous ones that
    `--annotation_in_groupTest` will drop later;
  * no NA and nothing non-numeric;
  * every value in [0, 1].

Score construction, following the paper (Schwartzentruber et al., bioRxiv
2025.11.04.686320) and `docs/design/flexrv-design.md` A2:
  * LoF (annotation in --lof-anno)        -> 1.0 by construction, NOT from
                                             AlphaMissense; the paper's LoF
                                             mask is annotation-based.
  * missense (annotation in --missense-anno)
                                          -> am_pathogenicity, joined on
                                             (CHROM, POS, REF, ALT).
  * missense with no AlphaMissense row    -> --missing policy. Default is the
                                             gene mean of that gene's scored
                                             missense variants (the design
                                             doc's suggested default; 0 =
                                             "treat as benign" is the other,
                                             and the paper does not say).
  * everything else (synonymous, non_coding, ...)
                                          -> --other (default 0.0). These are
                                             never tested -- markers matching
                                             no requested annotation are
                                             dropped in R before any decode --
                                             but the parser still demands a
                                             value, so this is a placeholder,
                                             not a claim about the variant.

TRANSCRIPTS. `AlphaMissense_hg38.tsv.gz` covers 19k genes on ONE transcript
each -- GENCODE v32 canonical, as defined in the AlphaMissense publication.
It is NOT MANE Select, and `AlphaMissense_isoforms_hg38.tsv.gz` is the
non-canonical COMPLEMENT (60k isoforms), which the release explicitly warns
"were not thoroughly evaluated and should be used with caution". So
"MANE Select, canonical where no MANE Select exists" is not a selection the
release offers directly. Two honest routes:
  * default: take the canonical file as shipped, and pass --mane <MANE summary>
    to REPORT how many of the joined transcripts are MANE Select. The report is
    the number that says whether the distinction mattered here.
  * --isoforms <file> --mane <file>: prefer the MANE Select transcript's row,
    drawing it from the isoforms file where it is not the canonical one, and
    falling back to canonical otherwise. This is the literal request, and it
    inherits DeepMind's caution on every row it takes from the isoforms file.
    The run reports how many rows came from each source.
  * --isoforms --anno-table <BRaVa table> --prefer-anno-transcript (what
    04_flexrv_groupfile.sh does with --annoTable): the row scored on the
    transcript the annotation used. Where AlphaMissense has none, the
    canonical row, else a same-gene isoform row, else --missing; never a row
    on a transcript another gene uses (pick() in build() has the rule and the
    measurement behind it).

Join keys are built from the group file's own variant IDs; the format is
auto-detected from the first ID and can be forced with --id-format. Only the
keys the group file actually needs are retained while streaming AlphaMissense,
so memory is O(missense variants in the group file), not O(71M).

Self-test (no downloads, positive and negative controls):
    python3 harness/tools/flexrv_score_from_alphamissense.py --selftest
"""

import argparse
import gzip
import io
import os
import random
import sys
import tempfile

# ---------------------------------------------------------------- id parsing

# Group file variant IDs seen in this repo and in the AoU exome group files:
#   20:881:A:T      chr20:881:A:T      20:881_A/T      chr20:881_A/T
#   20-881-A-T      chr20-881-A-T
# Each yields (chrom, pos, ref, alt) with the chrom's "chr" prefix stripped.
# colon_underscore and underscore: BRaVa biobank group files, e.g. GEL/Qatar "chrX:100629994_TAAAG_T",
# Genes & Health "chr21_41801194_C_G", PMBB "1_169853716_C_A" (2026-10-10 audit).
ID_FORMATS = ("colon", "underscore_slash", "dash", "colon_underscore", "underscore")


def detect_id_format(vid):
    """Return one of ID_FORMATS, or None if the ID carries no coordinates."""
    if vid.count(":") >= 3:
        return "colon"
    if "_" in vid and "/" in vid and ":" in vid:
        return "underscore_slash"
    if vid.count("-") >= 3:
        return "dash"
    if vid.count(":") == 1 and vid.count("_") >= 2 and "/" not in vid:
        return "colon_underscore"
    if ":" not in vid and vid.count("_") >= 3:
        return "underscore"
    return None


def parse_variant_id(vid, fmt):
    """(chrom, pos, ref, alt) or None if this ID does not parse in `fmt`."""
    try:
        if fmt == "colon":
            chrom, pos, ref, alt = vid.split(":", 3)
        elif fmt == "underscore_slash":
            chrom, rest = vid.split(":", 1)
            pos, alleles = rest.split("_", 1)
            ref, alt = alleles.split("/", 1)
        elif fmt == "dash":
            chrom, pos, ref, alt = vid.split("-", 3)
        elif fmt == "colon_underscore":
            chrom, rest = vid.split(":", 1)
            pos, ref, alt = rest.split("_", 2)
        elif fmt == "underscore":
            chrom, pos, ref, alt = vid.split("_", 3)
        else:
            return None
    except ValueError:
        return None
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if not pos.isdigit():
        return None
    return (chrom, pos, ref.upper(), alt.upper())


# --------------------------------------------------------------- file access

def smart_open(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


# GRCh38 pseudoautosomal regions. AlphaMissense files a PAR gene's rows under chrY (transcript IDs '..._PAR_Y'),
# sometimes with no chrX copy at all (CRLF2, PLCXD1, ZBED1, CSF2RA, WASH6P), so chrX variants there found no score.
# PAR1 has the same coordinates on X and Y; PAR2 is X 155,701,383-156,030,895 = Y 56,887,903-57,217,415.
PAR1_END, PAR2_Y_START, PAR2_Y_END, PAR2_SHIFT = 2781479, 56887903, 57217415, 155701383 - 56887903
PAR_Y = {}   # AlphaMissense path -> the AlphaMissense file whose chrY '_PAR_Y' rows are read with it, as chrX (build())
EXTRA_ROWS = {}   # AlphaMissense path -> further files in its layout read after it (--extra-isoforms: paralog copies)


def par_y_lines(path):
    """The '_PAR_Y' rows of an AlphaMissense file, rewritten as chrX rows: chrom chrX, PAR2 positions moved to X, the
    '_PAR_Y' suffix dropped from the transcript ID (so it matches GENCODE's chrX transcript)."""
    with smart_open(path) as fh:
        it = None
        for raw in fh:
            if raw.startswith("#"):
                if raw.startswith("#CHROM"):
                    it = raw.rstrip("\n").split("\t").index("transcript_id")
                continue
            f = raw.rstrip("\n").split("\t")
            if it is None or not f[it].endswith("_PAR_Y"):
                continue
            pos = int(f[1])
            if pos > PAR1_END:
                if not PAR2_Y_START <= pos <= PAR2_Y_END:
                    continue
                pos += PAR2_SHIFT
            f[0], f[1], f[it] = "chrX", str(pos), f[it][:-len("_PAR_Y")]
            yield "\t".join(f) + "\n"


def smart_lines(path):
    with smart_open(path) as fh:
        for raw in fh:
            yield raw


class am_open(object):
    """smart_open for AlphaMissense files: the file's lines, then -- when PAR_Y names a source for it -- that
    source's PAR rows as chrX rows (after the header has been read, so the readers parse them like any other)."""
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.fh = smart_open(self.path)
        src = PAR_Y.get(self.path)
        import itertools
        parts = [self.fh] + ([par_y_lines(src)] if src else [])
        for x in EXTRA_ROWS.get(self.path, ()):
            parts.append(smart_lines(x))
        return itertools.chain(*parts) if len(parts) > 1 else self.fh

    def __exit__(self, *a):
        self.fh.close()


def read_group_regions(path):
    """Yield (gene, kind, name, values, raw) per line, preserving file order.

    `kind` is one of var / anno / score / weight / other; `name` is the set
    name for score and weight lines. Whitespace runs are the separator, which
    is what checkGroupFile uses (`[ \\t]+`).
    """
    with smart_open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                yield (None, "other", None, [], line)
                continue
            gene, typ = parts[0], parts[1]
            head = typ.split(":")[0].split(",")[0]
            name = None
            if head in ("score", "weight"):
                bits = typ.replace(",", ":").split(":", 1)
                name = bits[1] if len(bits) > 1 and bits[1] else head.upper()
                kind = head
            elif head in ("var", "anno"):
                kind = head
            else:
                kind = "other"
            yield (gene, kind, name, parts[2:], line)


# ---------------------------------------------------------- alphamissense io

AM_CANONICAL_COLS = ("#CHROM", "POS", "REF", "ALT", "genome", "uniprot_id",
                     "transcript_id", "protein_variant", "am_pathogenicity",
                     "am_class")


def stream_alphamissense(path, wanted, chrom_filter=None, source="canonical",
                         collect_all=False):
    """Return {key: (score, transcript_id)} for keys in `wanted`.

    Only keys the group file needs are kept, so memory is bounded by the group
    file, not by the 71M-row release. When two rows share a key -- which the
    canonical file can do where genes overlap -- the HIGHER am_pathogenicity is
    kept and the collision is counted, because the alternative (first wins)
    would depend on file order.
    """
    hits, collisions, rows = {}, 0, 0
    allrows = {} if collect_all else None
    with am_open(path) as fh:
        header = None
        for raw in fh:
            if raw.startswith("#"):
                if raw.startswith("#CHROM"):
                    header = raw.rstrip("\n").split("\t")
                continue
            if header is None:
                raise SystemExit(
                    "%s: no '#CHROM...' header line found; is this an "
                    "AlphaMissense release file?" % path)
            rows += 1
            f = raw.rstrip("\n").split("\t")
            chrom = f[0][3:] if f[0].lower().startswith("chr") else f[0]
            if chrom_filter is not None and chrom != chrom_filter:
                continue
            key = (chrom, f[1], f[2].upper(), f[3].upper())
            if key not in wanted:
                continue
            idx_t = header.index("transcript_id")
            idx_p = header.index("am_pathogenicity")
            score = float(f[idx_p])
            tx = f[idx_t]
            if collect_all:
                # UniProt accession, isoform suffix dropped; the isoforms file
                # has no such column. It is what ties two canonical-file
                # transcripts to one protein (NEB: ENST00000409198 and
                # ENST00000172853 are both P20929).
                up = (f[header.index("uniprot_id")].split("-")[0]
                      if "uniprot_id" in header else None)
                allrows.setdefault(key, []).append((score, tx, source, up))
            prev = hits.get(key)
            if prev is None:
                hits[key] = (score, tx)
            else:
                collisions += 1
                if score > prev[0]:
                    hits[key] = (score, tx)
    if collect_all:
        return hits, collisions, rows, allrows
    return hits, collisions, rows


def tx_version(tx):
    """Numeric version of an Ensembl transcript ID ('ENST...354772.10' -> 10), -1 if none."""
    v = tx.rsplit(".", 1)[1] if "." in tx else ""
    return int(v) if v.isdigit() else -1


def read_reference_scores(sources, transcripts, canonical_genes, chrom_filter=None, gtf_versions=None):
    """All AlphaMissense scores on the transcripts that define each gene's percentile scale.

    `sources`: [(path, "canonical" | "isoforms")]. `transcripts`: version-stripped IDs to collect
    from either file (the annotation's own transcripts). `canonical_genes`: {version-stripped
    transcript: gene} of canonical-file transcripts to collect as well, the fallback scale for a gene
    whose annotation transcript AlphaMissense does not score (pass {} to skip). Returns
    {version-stripped transcript: (source, full ID, sorted scores)}, one version per transcript: the
    canonical file's if it has one, else the highest version (the isoform file lists two versions of
    some transcripts with different scores, e.g. ENST00000354772.8 / .10): the version the annotation's GENCODE
    holds (`gtf_versions`), else the canonical file's, else the highest -- the same choice pick() makes for the score
    itself, so a matched variant is always ranked against its own row (2026-10-21: neither "first in the file" nor
    "highest" alone is the annotation's version; on AoU chr1 the two disagreed for 5,780 entries, v39 siding with the
    highest for 5,133 and with the first for 647).
    """
    got = {}  # stripped -> {(source, full): [scores]}
    for path, source in sources:
        with am_open(path) as fh:
            header = None
            for raw in fh:
                if raw.startswith("#"):
                    if raw.startswith("#CHROM"):
                        header = raw.rstrip("\n").split("\t")
                        i_t, i_p = header.index("transcript_id"), header.index("am_pathogenicity")
                    continue
                f = raw.rstrip("\n").split("\t")
                if chrom_filter is not None:
                    c = f[0][3:] if f[0].lower().startswith("chr") else f[0]
                    if c != chrom_filter:
                        continue
                t = f[i_t].split(".")[0]
                if t in transcripts or (source == "canonical" and t in canonical_genes):
                    got.setdefault(t, {}).setdefault((source, f[i_t]), []).append(float(f[i_p]))
    out = {}
    gv = gtf_versions or {}
    for t, by in got.items():
        (src, full), scores = max(by.items(), key=lambda kv: (gv.get(t) == kv[0][1], kv[0][0] == "canonical",
                                                             tx_version(kv[0][1])))
        out[t] = (src, full, sorted(scores))
    return out


def percentile_in(x, ref):
    """Percentile of score x among the sorted reference scores `ref`: the average 1-based rank of x
    (ties share their mean rank) divided by N, i.e. (#below + (#equal + 1) / 2) / N. For a score in
    `ref` this is exactly a within-gene percentile over all possible missense variants (FlexRV,
    Methods); a fallback score from another transcript is placed at its position on the same scale, and one
    above every score on it is at the top, 1.0 (unclamped it would be (N + 0.5) / N; 64 files of the first
    genome-wide build, 2026-10-21)."""
    import bisect
    lo, hi = bisect.bisect_left(ref, x), bisect.bisect_right(ref, x)
    return min(1.0, (lo + (hi - lo + 1) / 2.0) / len(ref))


def read_mane_transcripts(path):
    """Ensembl transcript IDs (version-stripped) of MANE Select from a MANE
    summary file (`MANE.GRCh38.vX.summary.txt.gz`). Falls back to scanning any
    column that looks like an ENST accession if the header is unexpected."""
    mane = set()
    with smart_open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ens = header.index("Ensembl_nuc")
            kind = header.index("MANE_status")
        except ValueError:
            ens, kind = None, None
        for raw in fh:
            f = raw.rstrip("\n").split("\t")
            if ens is not None:
                if kind is not None and f[kind] != "MANE Select":
                    continue
                mane.add(f[ens].split(".")[0])
            else:
                for cell in f:
                    if cell.startswith("ENST"):
                        mane.add(cell.split(".")[0])
    return mane


def read_mane_genes(path):
    """{Ensembl transcript (version-stripped): {Ensembl gene, symbol}} for
    EVERY transcript in a MANE summary (Select and Plus Clinical), so a row's
    transcript can be checked against the gene it is being used for."""
    genes = {}
    with smart_open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_t, i_g = header.index("Ensembl_nuc"), header.index("Ensembl_Gene")
            i_s = header.index("symbol")
        except ValueError as e:
            raise SystemExit("%s: expected Ensembl_nuc/Ensembl_Gene/symbol columns "
                             "of a MANE summary (%s)" % (path, e))
        for raw in fh:
            f = raw.rstrip("\n").split("\t")
            genes[f[i_t].split(".")[0]] = {f[i_g].split(".")[0], f[i_s]}
    return genes


def read_gtf_genes(path, chrom_filter=None):
    """{Ensembl transcript: Ensembl gene}, both version-stripped, from the
    `transcript` lines of a GENCODE GTF. It must be the release the annotation
    table was built on (GENCODE v39 for AoU v8: the earliest release holding
    every one of its 19,543 gene and transcript IDs; v38 lacks 22 genes),
    so that gene IDs compare without translation. An AlphaMissense transcript
    (GENCODE v32) retired by then is simply absent: no evidence either way."""
    genes = {}
    with smart_open(path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            f = raw.split("\t", 8)
            if len(f) < 9 or f[2] != "transcript":
                continue
            if chrom_filter is not None:
                c = f[0][3:] if f[0].lower().startswith("chr") else f[0]
                if c != chrom_filter:
                    continue
            g = t = None
            for kv in f[8].split(";"):
                kv = kv.strip()
                if kv.startswith("gene_id "):
                    g = kv[8:].strip('"').split(".")[0]
                elif kv.startswith("transcript_id "):
                    t = kv[14:].strip('"').split(".")[0]
            if g and t:
                genes[t] = g
    return genes


def read_gtf_versions(path, chrom_filter=None):
    """{Ensembl transcript (version-stripped): its full versioned ID} from the GTF's `transcript` lines: which version of
    a transcript the annotation's GENCODE release holds."""
    out = {}
    with smart_open(path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            f = raw.split("\t", 8)
            if len(f) < 9 or f[2] != "transcript":
                continue
            if chrom_filter is not None:
                c = f[0][3:] if f[0].lower().startswith("chr") else f[0]
                if c != chrom_filter:
                    continue
            i = f[8].find('transcript_id "')
            if i >= 0:
                t = f[8][i + 15:f[8].index('"', i + 15)]
                out[t.split(".")[0]] = t
    return out


def read_gtf_spans(path, genes, chrom_filter=None):
    """{Ensembl gene (version-stripped): (start, end)} over the `transcript` lines of a GENCODE GTF, for `genes`."""
    span = {}
    with smart_open(path) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            f = raw.split("\t", 8)
            if len(f) < 9 or f[2] != "transcript":
                continue
            if chrom_filter is not None:
                c = f[0][3:] if f[0].lower().startswith("chr") else f[0]
                if c != chrom_filter:
                    continue
            i = f[8].find('gene_id "')
            g = f[8][i + 9:f[8].index('"', i + 9)].split(".")[0] if i >= 0 else None
            if g in genes:
                lo, hi = span.get(g, (int(f[3]), int(f[4])))
                span[g] = (min(lo, int(f[3])), max(hi, int(f[4])))
    return span


def read_unmapped_in_spans(sources, spans, mapped, chrom_filter=None):
    """{gene: {version-stripped transcript: (source, full ID, sorted scores)}} for AlphaMissense transcripts that
    GENCODE does not hold (`mapped`: its transcripts) with rows inside a gene's span -- the scale candidates for a gene
    whose own transcripts AlphaMissense does not score, chosen from AlphaMissense and GENCODE alone."""
    import bisect
    order = sorted((lo, hi, g) for g, (lo, hi) in spans.items()); starts = [o[0] for o in order]
    got = {}
    for path, source in sources:
        with am_open(path) as fh:
            for raw in fh:
                if raw.startswith("#"):
                    if raw.startswith("#CHROM"):
                        h = raw.rstrip("\n").split("\t"); i_t, i_p = h.index("transcript_id"), h.index("am_pathogenicity")
                    continue
                f = raw.rstrip("\n").split("\t")
                if chrom_filter is not None:
                    c = f[0][3:] if f[0].lower().startswith("chr") else f[0]
                    if c != chrom_filter:
                        continue
                t = f[i_t].split(".")[0]
                if t in mapped:
                    continue
                pos = int(f[1]); j = bisect.bisect_right(starts, pos)
                for lo, hi, g in order[max(0, j - 50):j]:
                    if lo <= pos <= hi:
                        got.setdefault(g, {}).setdefault(t, {}).setdefault((source, f[i_t]), []).append(float(f[i_p]))
    out = {}
    for g, by_t in got.items():
        for t, by in by_t.items():
            (src, full), sc = max(by.items(), key=lambda kv: (kv[0][0] == "canonical", tx_version(kv[0][1])))
            out.setdefault(g, {})[t] = (src, full, sorted(sc))
    return out


def read_anno_table(path, id_col="ID", gene_col="GENE", tx_col="TRANSCRIPT",
                    mane_col="MANE_SELECT", csq_col="CSQ", score_col=None,
                    zero_is_missing=True, keep_ids=None):
    """Per-variant transcript context from a BRaVa long-form annotation table.

    That table is the SOURCE of the group file -- one row per (gene, variant),
    carrying the transcript each consequence was called against, the MANE
    Select accession when the transcript has one, and the raw VEP consequence.
    It is what makes "the MANE Select transcript, or canonical where there is
    none" a checkable claim rather than an assumption: AlphaMissense names the
    transcript it scored, and this names the transcript the annotation used.

    Returns {(gene, variant_id): (transcript, mane_accession_or_None, csq)}.
    `keep_ids`: keep only rows at these variant IDs (every gene, every
    consequence). Every lookup is at a missense-bin variant, so passing that
    set changes no output; it drops the non_coding bulk, which on AoU is most
    of the table and made 4 large chromosomes at once swap a 16 GB machine.
    """
    ctx = {}
    with smart_open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_id, i_g = header.index(id_col), header.index(gene_col)
            i_t, i_m = header.index(tx_col), header.index(mane_col)
            i_c = header.index(csq_col)
            i_x = header.index(score_col) if score_col else None
        except ValueError as e:
            raise SystemExit("%s: expected columns %s; got %s (%s)" % (
                path, [id_col, gene_col, tx_col, mane_col, csq_col], header, e))
        for raw in fh:
            if keep_ids is not None and raw.split("\t", i_id + 1)[i_id] not in keep_ids:
                continue
            f = raw.rstrip("\n").split("\t")
            mane = f[i_m]
            extra = None
            if i_x is not None:
                raw = f[i_x]
                try:
                    extra = float(raw)
                except ValueError:
                    extra = None
                if extra is not None and zero_is_missing and extra == 0.0:
                    extra = None
            ctx[(f[i_g], f[i_id])] = (
                f[i_t], None if mane in ("NA", "", ".") else mane, f[i_c], extra)
    return ctx


def is_snv(vid, fmt):
    key = parse_variant_id(vid, fmt)
    return key is not None and len(key[2]) == 1 and len(key[3]) == 1


# Consequences that change the protein but are not a missense substitution: under --nonmissense
# gene-mean these stay in the FlexRV file at the gene mean (in-frame indels, LOFTEE-LC LoF types,
# start/stop lost). A missense-bin entry with none of them (a SpliceAI-route intronic, splice-region,
# synonymous or UTR variant) has no protein change for AlphaMissense to speak to and is dropped, as in
# the FlexRV paper and the group files sent to BRaVa (2026-10-20), which hold no such variants.
PROTEIN_ALTERING_NONMISSENSE = {
    "inframe_insertion", "inframe_deletion", "protein_altering_variant", "start_lost", "stop_lost",
    "frameshift_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant",
    "transcript_ablation"}


def is_protein_altering_csq(csq):
    return csq is not None and bool(PROTEIN_ALTERING_NONMISSENSE.intersection(csq.split("&")))


def is_missense_csq(csq):
    """A CSQ is missense when missense_variant is one of its '&'-joined VEP terms: an annotation that writes
    every term (e.g. the BRaVa 2026-10 re-annotation: missense_variant&splice_region_variant) is as missense as
    one that writes only the most severe (AoU v8: missense_variant). Before 2026-10-15 the test was equality,
    and --nonmissense drop removed every missense variant in a splice region from such annotations."""
    return csq is not None and "missense_variant" in csq.split("&")


SCORE_COLUMNS = {
    # column -> (default transform, whether 0 means "absent", one-line note)
    "REVEL_SCORE": ("none", True,
                    "REVEL, already a [0,1] ensemble probability; missense only "
                    "by construction, so it is 0 for every pLoF row"),
    "CADD_PHRED":  ("percentile-in-gene", True,
                    "CADD PHRED is UNBOUNDED (0.001-63 here), so it needs a "
                    "transform; it is genome-wide per variant, NOT "
                    "transcript-specific (measured: 49,641 multi-gene pairs, 0 "
                    "disagreements)"),
    # 0 is a REAL spliceAI output ("no predicted splice effect"), not a gap:
    # 89.8% of the nonzero values on chr20 sit BELOW the 0.20 threshold the
    # file is named for, so the column is populated across the whole range
    # rather than only where it drove a call. REVEL and CADD are the opposite:
    # REVEL is 0 for 100% of pLoF/synonymous/non_coding rows and CADD for 100%
    # of synonymous, which is a pipeline gap, not a measurement.
    "DS_MAX":      ("none", False,
                    "spliceAI max delta, already [0,1]; scored against ONE gene "
                    "per variant (measured: never nonzero in two genes at once)"),
}


def apply_transform(vals, transform):
    """vals is a list of (index, raw) for ONE gene. Returns {index: [0,1] score}.

    percentile-in-gene is the paper's own construction for PrimateAI-3D -- a
    within-gene rank, which gives the transforms real spread. cadd-phred maps a
    PHRED to the genome-wide percentile it encodes (1 - 10^(-p/10)); it is
    principled but compresses everything above PHRED 10 into [0.9, 1], which is
    most of what a damaging-missense mask contains, so it is NOT the default.
    """
    if not vals:
        return {}
    if transform == "none":
        out = {}
        for i, v in vals:
            if not (0.0 <= v <= 1.0):
                raise SystemExit(
                    "score %.6g is outside [0, 1] and --score-transform is "
                    "'none'; pick a transform" % v)
            out[i] = v
        return out
    if transform == "cadd-phred":
        return {i: 1.0 - 10.0 ** (-v / 10.0) for i, v in vals}
    if transform == "percentile-in-gene":
        # average ranks for ties, mapped onto (0, 1]; a gene with one scored
        # variant gets 1.0, which is what a percentile of a single value is.
        order = sorted(vals, key=lambda t: t[1])
        n, out, i = len(order), {}, 0
        while i < n:
            j = i
            while j + 1 < n and order[j + 1][1] == order[i][1]:
                j += 1
            rank = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k][0]] = rank / n
            i = j + 1
        return out
    if transform.startswith("minmax:"):
        lo, hi = (float(x) for x in transform.split(":", 1)[1].split(","))
        if hi <= lo:
            raise SystemExit("minmax needs LO < HI")
        return {i: min(1.0, max(0.0, (v - lo) / (hi - lo))) for i, v in vals}
    raise SystemExit("unknown --score-transform %r" % transform)


# ------------------------------------------------------------------- builder

def build(group_in, group_out, am_path, lof_anno, missense_anno,
          score_name="AM", missing="gene-mean", other=0.0, id_format=None,
          chrom_filter=None, isoforms_path=None, mane_path=None,
          replace_existing=False, anno_table_path=None,
          prefer_anno_transcript=False, drop_anno=(), score_column=None,
          score_transform=None, zero_is_missing=None, fallback="safe",
          provenance=None, gtf_path=None, nonmissense="impute", am_scale="raw",
          par_y_am=None, par_y_isoforms=None, empty_gene_fill=None, extra_isoforms=(), log=sys.stderr):
    """Write `group_out` = `group_in` with a score line per region.

    `am_scale`: "raw" writes am_pathogenicity itself. "percentile" (2026-10-21) writes the FlexRV
    paper's score, the variant's percentile among ALL possible missense variants of its gene --
    here, every AlphaMissense row on the gene's annotation transcript (MANE Select, else canonical:
    the TRANSCRIPT most of the gene's missense-bin rows carry in --anno-table), from whichever file
    holds it; where AlphaMissense does not score that transcript, the gene's canonical-file
    transcript (mapped to the gene by --gtf). The row a variant's score comes from is chosen exactly
    as in raw mode (pick()); only the scale changes. Needs --anno-table, --prefer-anno-transcript,
    --isoforms and --gtf. --missing / --nonmissense gene-mean fills are taken on the percentile scale,
    after ranking, so a fill never moves a real variant's rank.

    `fallback` is what --prefer-anno-transcript does when AlphaMissense has no
    row on the annotation's transcript: "safe" (see pick()), or "max", the rule
    before 2026-10, kept only so files built with it can be reproduced.
    `gtf_path`: GENCODE GTF of the annotation's release; a fallback row whose
    transcript it puts in another gene is never used (see pick()).
    `nonmissense`: what to do with a missense-bin variant whose annotated consequence
    (--anno-table CSQ) is not missense_variant -- in-frame indels, LOFTEE-LC
    frameshifts/stop-gains, start/stop-lost, SpliceAI-route splice/intronic/
    synonymous variants. "impute" (default, as before): AlphaMissense if a row
    exists on some transcript, else the --missing policy. "drop": remove it from
    the output, because AlphaMissense scores missense substitutions only, so any
    value it got would be another transcript's score or a gene-mean placeholder
    (AoU v8 chr20: 27,911 such entries, 27,421 at the gene mean). "gene-mean"
    (2026-10-21): keep it, with AlphaMissense only from a row on its OWN annotation
    transcript (a start-lost substitution can have one), else the gene mean of the
    gene's scored missense entries -- the FlexRV group files sent to BRaVa
    (2026-10-20) do exactly this for in-frame indels, LOFTEE-LC LoF and
    start/stop-lost (82 of 82 chr21 genes reproduced), and it keeps the FlexRV
    run's variants those of the SKAT-O masks, LC LoF included without calling it
    pLoF. Entries whose consequence changes no protein (SpliceAI-route intronic,
    splice-region, synonymous) are dropped, as under "drop".
    `provenance`, a path, gets one line per missense (gene, variant): the
    source of its score, the transcript(s) and the value before imputation.

    Returns a dict of counts -- the report is part of the output, because
    "how many missense variants got a real AlphaMissense score" is the number
    that decides whether the run means anything.
    """
    PAR_Y.clear(); EXTRA_ROWS.clear()
    extra_isoforms = [x for x in (extra_isoforms or ()) if x]
    if extra_isoforms:
        if not isoforms_path:
            raise SystemExit("--extra-isoforms needs --isoforms: its rows are read as isoform rows")
        EXTRA_ROWS[isoforms_path] = extra_isoforms
    if chrom_filter in (None, "X"):
        if par_y_am and am_path:
            PAR_Y[am_path] = par_y_am
        if par_y_isoforms and isoforms_path:
            PAR_Y[isoforms_path] = par_y_isoforms
    empty_fill = (0.5 if am_scale == "percentile" else 0.0) if empty_gene_fill is None else float(empty_gene_fill)
    lof = set(lof_anno)
    mis = set(missense_anno)

    regions, order = {}, []
    for gene, kind, name, values, raw in read_group_regions(group_in):
        if gene is None:
            continue
        if gene not in regions:
            regions[gene] = {"var": None, "anno": None, "lines": []}
            order.append(gene)
        r = regions[gene]
        if kind == "var":
            r["var"] = values
        elif kind == "anno":
            r["anno"] = values
        elif kind == "score":
            if not replace_existing:
                raise SystemExit(
                    "%s already carries a score line ('score:%s' in region %s). "
                    "Pass --replace-existing to overwrite it." % (group_in, name, gene))
            continue  # dropped; the new one is written in its place
        r["lines"].append((kind, name, values, raw))

    for gene in order:
        r = regions[gene]
        if r["var"] is None or r["anno"] is None:
            raise SystemExit("region %s has no %s line" %
                             (gene, "var" if r["var"] is None else "anno"))
        if len(r["var"]) != len(r["anno"]):
            raise SystemExit("region %s: var line has %d entries, anno line %d"
                             % (gene, len(r["var"]), len(r["anno"])))

    # --- which annotation labels exist, so a typo in --lof-anno is caught here
    seen_anno = set()
    for gene in order:
        seen_anno.update(regions[gene]["anno"])
    for label, flag in [(l, "--lof-anno") for l in lof] + \
                       [(m, "--missense-anno") for m in mis]:
        if label not in seen_anno:
            raise SystemExit(
                "%s=%s: no variant in %s carries that annotation. Labels present: %s"
                % (flag, label, group_in, ", ".join(sorted(seen_anno))))

    counts_base = {"genes": len(order)}
    zim = SCORE_COLUMNS.get(score_column, (None, True))[1]
    if zero_is_missing is not None:
        zim = zero_is_missing
    mis_ids = set(v for g in order
                  for v, a in zip(regions[g]["var"], regions[g]["anno"]) if a in mis)
    anno_ctx = (read_anno_table(anno_table_path, score_col=score_column,
                                zero_is_missing=zim, keep_ids=mis_ids)
                if anno_table_path else None)

    # --- the keys we need: missense variants only
    if id_format is None:
        probe = next((regions[g]["var"][0] for g in order if regions[g]["var"]), None)
        id_format = detect_id_format(probe) if probe else None
        if id_format is None:
            raise SystemExit(
                "could not detect a coordinate format from the first variant ID "
                "(%r). AlphaMissense joins on (CHROM, POS, REF, ALT), so IDs "
                "like rsIDs cannot be used; pass --id-format if the format is "
                "one of %s." % (probe, "/".join(ID_FORMATS)))
        print("# id format auto-detected: %s (from %r)" % (id_format, probe), file=log)

    if score_column and anno_ctx is None:
        raise SystemExit("--score-column needs --anno-table: the column lives there")

    wanted, unparsed = set(), 0
    for gene in order:
        r = regions[gene]
        for vid, anno in zip(r["var"], r["anno"]):
            if anno not in mis:
                continue
            key = parse_variant_id(vid, id_format)
            if key is None:
                unparsed += 1
                continue
            wanted.add(key)
    if not wanted:
        raise SystemExit(
            "no missense variant IDs parsed into coordinates; nothing to join. "
            "Check --missense-anno (%s) and --id-format (%s)."
            % (",".join(sorted(mis)), id_format))

    if score_column:
        # Column mode: the score already sits in the annotation table, one
        # value per (gene, variant), so there is nothing to join and no
        # transcript to reconcile -- whatever transcript-specificity the score
        # has is already baked into the row.
        return _build_from_column(
            regions, order, lof, mis, drop_anno, group_out, score_name,
            anno_ctx, score_column, score_transform, missing, other, counts_base,
            zim, log, nonmissense=nonmissense)

    # --- join. With an annotation table we keep EVERY AlphaMissense row per
    # coordinate, because the row to use is decided per (gene, variant) by the
    # transcript the annotation itself chose -- not by a global max.
    prefer_tx = anno_ctx is not None and prefer_anno_transcript
    if prefer_tx:
        hits, collisions, rows, allrows = stream_alphamissense(
            am_path, wanted, chrom_filter, collect_all=True)
    else:
        hits, collisions, rows = stream_alphamissense(am_path, wanted, chrom_filter)
        allrows = None
    print("# AlphaMissense canonical: %d of %d wanted keys matched (%d rows scanned, "
          "%d key collisions resolved by max)" % (len(hits), len(wanted), rows, collisions),
          file=log)

    from_isoform = 0
    mane = read_mane_transcripts(mane_path) if mane_path else None
    if isoforms_path:
        if mane is None and not prefer_tx:
            raise SystemExit(
                "--isoforms needs a rule for which isoform to prefer: either "
                "--mane (the MANE summary) or --prefer-anno-transcript with "
                "--anno-table (the transcript the annotation itself used).")
        if prefer_tx:
            iso_hits, iso_coll, iso_rows, iso_all = stream_alphamissense(
                isoforms_path, wanted, chrom_filter, source="isoforms",
                collect_all=True)
            for k, v in iso_all.items():
                allrows.setdefault(k, []).extend(v)
        else:
            iso_hits, iso_coll, iso_rows = stream_alphamissense(
                isoforms_path, wanted, chrom_filter, source="isoforms")
        for key, (score, tx) in (iso_hits.items() if not prefer_tx else []):
            if tx.split(".")[0] not in mane:
                continue
            cur = hits.get(key)
            if cur is not None and cur[1].split(".")[0] in mane:
                continue  # canonical already IS the MANE Select transcript
            hits[key] = (score, tx)
            from_isoform += 1
        if prefer_tx:
            print("# AlphaMissense isoforms: %d rows scanned and merged; the row "
                  "per (gene, variant) is chosen by transcript below"
                  % iso_rows, file=log)
        else:
            print("# AlphaMissense isoforms: %d keys replaced by their MANE Select "
                  "transcript (%d rows scanned)" % (from_isoform, iso_rows), file=log)

    mane_hits = None
    if mane is not None:
        mane_hits = sum(1 for _, tx in hits.values() if tx.split(".")[0] in mane)
        print("# MANE Select: %d of %d joined transcripts are MANE Select (%.1f%%)"
              % (mane_hits, len(hits), 100.0 * mane_hits / max(1, len(hits))), file=log)

    # --- the row to use, chosen per (gene, variant) ------------------------
    # "the MANE Select transcript, or the canonical one where no MANE Select
    # exists" is a choice the ANNOTATION already made: its TRANSCRIPT column is
    # that transcript, and its MANE_SELECT column says which of the two cases
    # applies. So the correct AlphaMissense row is the one scored against THAT
    # transcript -- wherever it lives, canonical file or isoforms file. Keying
    # on coordinates alone would let a neighbouring gene's MANE transcript win
    # at a shared position, which is real here: 606 coordinates carry more than
    # one row.
    #
    # When AlphaMissense has no row on that transcript, the fallback used to be
    # the max over EVERY row at the coordinate, any transcript of any gene.
    # That is biased upward (a max over correlated scores) and borrows an
    # overlapping gene's score: on AoU chr20, against the shipped-canonical
    # build, changed rows moved up 63.5% / down 36.5%, 0.906 was crossed 1,007
    # times upward and 411 downward, and 392 of 4,172 variants in two genes
    # with different transcripts got the identical score in both
    # (genebass-burden controls/AM_ACMG_CUTOFF.md, 2026-10-08). Now, in order:
    #   a. the canonical-file row (the transcript DeepMind evaluated). That file
    #      is one transcript per UniProt entry, not per gene, so a coordinate
    #      can carry several: rows sharing ONE accession are one protein and
    #      are averaged (NEB, 5,904 chr2 variants); rows of different
    #      accessions are different genes, ambiguous, and fall through to b;
    #   b. the eligible isoform rows, restricted to those the MANE summary puts
    #      in this gene when there are any; several are AVERAGED, a choice that
    #      does not look at the scores (a max would reintroduce the bias);
    #   c. nothing, so the --missing policy applies.
    # A row is ineligible when its transcript is the one the annotation assigns
    # to a DIFFERENT gene at this coordinate, or (with --mane) the MANE summary,
    # or (with --gtf) the annotation's own GENCODE release, puts it in a
    # different gene. The GTF is what catches a partner gene's NON-annotated
    # transcript: on AoU chr20, 1,559 fallback rows took one -- readthrough and
    # overlapping genes (AL031681.2 <- L3MBTL1/SRSF6, AL121845.3 <- LIME1/ZGPAT)
    # -- which neither the annotation nor MANE marks. Genes are compared by
    # Ensembl ID, not name, deliberately: in v39 ENSG00000105501 and
    # ENSG00000268500 are both "SIGLEC5" and both tested; AlphaMissense only
    # has the latter's transcripts, so the former is left unscored rather than
    # duplicating it (decided 2026-10-08). A gene left with no
    # scored missense variant is named in the log. No rule here picks a row by
    # its score.
    pick_stats = dict.fromkeys(
        ("by_transcript", "fallback_canonical", "fallback_isoform",
         "fallback_canonical_averaged", "fallback_isoform_averaged",
         "missing_no_row", "missing_excluded",
         "fallback_max"), 0)
    mane_genes = read_mane_genes(mane_path) if (mane_path and prefer_tx) else None
    gtf_genes = read_gtf_genes(gtf_path, chrom_filter) if (gtf_path and prefer_tx) else None
    gtf_versions = read_gtf_versions(gtf_path, chrom_filter) if (gtf_path and prefer_tx and am_scale == "percentile") else {}
    if gtf_genes is not None:
        print("# GTF: %d transcripts mapped to genes%s" % (
            len(gtf_genes), "" if chrom_filter is None else " on chr%s" % chrom_filter),
            file=log)
    excluded_by_gene = {}
    tx_at = {}  # variant -> {gene: annotation transcript}, version-stripped
    if prefer_tx:
        for (g, vid), ctx in anno_ctx.items():
            if ctx[0] and ctx[0] not in ("NA", "."):
                tx_at.setdefault(vid, {})[g] = ctx[0].split(".")[0]

    def pick(gene, vid, key):
        """(score, transcript, matched_annotation_transcript, source) or None."""
        if key is None:
            return None
        if not prefer_tx:
            h = hits.get(key)
            return None if h is None else (h[0], h[1], None, "canonical")
        want = anno_ctx.get((gene, vid), (None, None, None, None))[0]
        cands = allrows.get(key)
        if not cands:
            pick_stats["missing_no_row"] += 1
            return None
        if want:
            w = want.split(".")[0]
            if am_scale == "percentile":
                # one version per transcript, the one read_reference_scores() ranks against
                same = [c for c in cands if c[1].split(".")[0] == w]
                if same:
                    score, tx = max(((c[0], c[1]) for c in same),
                                    key=lambda st: (gtf_versions.get(w) == st[1],
                                                    any(c[2] == "canonical" and c[1] == st[1] for c in same),
                                                    tx_version(st[1])))
                    pick_stats["by_transcript"] += 1
                    return (score, tx, True, "matched")
            for score, tx, src, _ in cands:
                if tx.split(".")[0] == w:
                    pick_stats["by_transcript"] += 1
                    return (score, tx, True, "matched")
        if fallback == "max":
            best = max(cands, key=lambda c: c[0])
            pick_stats["fallback_max"] += 1
            return (best[0], best[1], False, "max")
        others = set(t for g, t in tx_at.get(vid, {}).items() if g != gene)

        ok = []
        for score, tx, src, up in cands:
            t = tx.split(".")[0]
            mg = mane_genes.get(t) if mane_genes is not None else None
            gg = gtf_genes.get(t) if gtf_genes is not None else None
            if (t in others or (mg is not None and gene not in mg)
                    or (gg is not None and gg != gene)):
                continue
            ok.append((score, tx, src, mg is not None or gg is not None, up))
        canon = [c for c in ok if c[2] == "canonical"]
        if canon and len(set(c[4] for c in canon)) == 1:
            pick_stats["fallback_canonical"] += 1
            if len(canon) > 1:
                pick_stats["fallback_canonical_averaged"] += 1
            return (sum(c[0] for c in canon) / len(canon),
                    ",".join(sorted(c[1] for c in canon)), False, "canonical")
        iso = [c for c in ok if c[2] == "isoforms"]
        if any(c[3] for c in iso):
            iso = [c for c in iso if c[3]]
        if iso:
            pick_stats["fallback_isoform"] += 1
            if len(iso) > 1:
                pick_stats["fallback_isoform_averaged"] += 1
            return (sum(c[0] for c in iso) / len(iso),
                    ",".join(sorted(c[1] for c in iso)), False, "isoform")
        pick_stats["missing_excluded"] += 1
        excluded_by_gene[gene] = excluded_by_gene.get(gene, 0) + 1
        return None

    # --- percentile scale: each gene's reference transcript and its scores --
    ref_of, ref_stats = {}, {"annotation_transcript": 0, "canonical_fallback": 0, "isoform_fallback": 0,
                             "used_transcript_fallback": 0, "none": 0}
    if am_scale == "percentile":
        if not (prefer_tx and isoforms_path and gtf_genes is not None):
            raise SystemExit("--am-scale percentile needs --anno-table with --prefer-anno-transcript, "
                             "--isoforms and --gtf: the scale is the gene's annotation transcript, which "
                             "may be in either AlphaMissense file, and its fallback is mapped by the GTF")
        from collections import Counter
        anno_tx = {}
        for gene in order:
            r = regions[gene]
            c = Counter((anno_ctx.get((gene, v)) or ("NA",))[0].split(".")[0]
                        for v, a in zip(r["var"], r["anno"]) if a in mis)
            c.pop("NA", None); c.pop(".", None); c.pop("", None)
            if c:  # the transcript most of the gene's rows were called on; ties -> smallest ID
                anno_tx[gene] = min(c, key=lambda t: (-c[t], t))
        genes_here = set(order)
        canon_genes = {t: g for t, g in gtf_genes.items() if g in genes_here}
        ref = read_reference_scores([(am_path, "canonical"), (isoforms_path, "isoforms")],
                                    set(anno_tx.values()), canon_genes, chrom_filter, gtf_versions)
        canon_by_gene = {}
        for t, (src, full, scores) in ref.items():
            if src == "canonical" and t in canon_genes:
                canon_by_gene.setdefault(canon_genes[t], []).append((len(scores), t))
        # third, as in pick(): a gene whose annotation transcript and canonical-file transcripts AlphaMissense does
        # not score, but whose OTHER isoforms it does (its scores then come from pick()'s isoform fallback), is
        # ranked on the isoform transcript scoring the most variants; a second pass reads only those genes'
        # transcripts (2026-10-21: without it 275,626 entries with an isoform score lost it genome-wide)
        need = {g for g in order if anno_tx.get(g) not in ref and not canon_by_gene.get(g)}
        iso_by_gene = {}
        if need:
            iso_tx = {t: g for t, g in gtf_genes.items() if g in need}
            iso_ref = read_reference_scores([(isoforms_path, "isoforms")], set(iso_tx), {}, chrom_filter, gtf_versions)
            for t, v in iso_ref.items():
                iso_by_gene.setdefault(iso_tx[t], []).append((len(v[2]), t, v))
        # fourth: a gene still without a scale whose scores come from AlphaMissense transcripts GENCODE v39 no longer
        # holds (AlphaMissense is on v32; pick()'s fallbacks accept them) is ranked on the largest such transcript with
        # rows inside the gene's v39 span -- chosen from AlphaMissense and GENCODE alone, never from a cohort's variants,
        # so one gene has one scale in every cohort (2026-10-21: choosing by the cohort's own picks made 3 genes /
        # 88 (gene, variant) differ between biobanks; without this fallback 45,532 entries in ~100 genes were filled
        # with 0.0 where the raw build gave them real scores)
        left = [g for g in order if anno_tx.get(g) not in ref and not canon_by_gene.get(g) and not iso_by_gene.get(g)]
        unmapped_by_gene = {}
        if left:
            spans = read_gtf_spans(gtf_path, set(left), chrom_filter)
            unmapped_by_gene = read_unmapped_in_spans([(am_path, "canonical"), (isoforms_path, "isoforms")],
                                                      spans, set(gtf_genes), chrom_filter)
        for gene in order:
            t = anno_tx.get(gene)
            if t in ref:
                ref_of[gene] = (t,) + ref[t][1:]
                ref_stats["annotation_transcript"] += 1
            elif canon_by_gene.get(gene):
                n, t2 = max(canon_by_gene[gene], key=lambda x: (x[0], x[1]))  # the one scoring the most variants
                ref_of[gene] = (t2,) + ref[t2][1:]
                ref_stats["canonical_fallback"] += 1
            elif iso_by_gene.get(gene):
                n, t2, v = max(iso_by_gene[gene], key=lambda x: (x[0], x[1]))
                ref_of[gene] = (t2,) + v[1:]
                ref_stats["isoform_fallback"] += 1
            elif unmapped_by_gene.get(gene):
                cand = unmapped_by_gene[gene]
                t2 = max(cand, key=lambda t: (len(cand[t][2]), t))
                ref_of[gene] = (t2,) + cand[t2][1:]
                ref_stats["used_transcript_fallback"] += 1
            else:
                ref_stats["none"] += 1
        print("# percentile scale: %d genes on their annotation transcript, %d on the canonical-file "
              "transcript (AlphaMissense does not score the annotation's), %d on another isoform (neither), "
              "%d on the largest AlphaMissense transcript in its span that GENCODE no longer holds, %d with no AlphaMissense "
              "transcript (their missense take --missing)" % (
                  ref_stats["annotation_transcript"], ref_stats["canonical_fallback"],
                  ref_stats["isoform_fallback"], ref_stats["used_transcript_fallback"], ref_stats["none"]), file=log)

    # --- transcript diagnostics, and an EXPLAINED coverage gap ------------
    # A bare coverage percentage is not interpretable: AlphaMissense scores
    # single-nucleotide MISSENSE substitutions only, while a BRaVa
    # "*_missense" bin also holds in-frame indels and variants whose most
    # severe consequence is something else entirely. So the unjoined set is
    # broken down by cause, and the number to read is the coverage of the
    # AM-ELIGIBLE subset, not of the bin.
    diag = {"eligible": 0, "eligible_joined": 0,
            "unjoined_indel": 0, "unjoined_not_missense_csq": 0,
            "unjoined_eligible": 0,
            "tx_agree": 0, "tx_disagree": 0, "tx_unknown": 0,
            "joined_mane": 0, "joined_no_mane": 0}
    for gene in order:
        r = regions[gene]
        for vid, anno in zip(r["var"], r["anno"]):
            if anno not in mis:
                continue
            key = parse_variant_id(vid, id_format)
            hit = pick(gene, vid, key)
            csq = anno_ctx.get((gene, vid), (None, None, None, None))[2] if anno_ctx else None
            snv = is_snv(vid, id_format)
            eligible = snv and (csq is None or is_missense_csq(csq))
            if eligible:
                diag["eligible"] += 1
                if hit is not None:
                    diag["eligible_joined"] += 1
                else:
                    diag["unjoined_eligible"] += 1
            elif hit is None:
                if not snv:
                    diag["unjoined_indel"] += 1
                else:
                    diag["unjoined_not_missense_csq"] += 1
            if hit is None or anno_ctx is None:
                continue
            tx_aou, mane_acc, _, _ = anno_ctx.get((gene, vid), (None, None, None, None))
            if tx_aou is None:
                diag["tx_unknown"] += 1
            elif tx_aou.split(".")[0] == hit[1].split(".")[0]:
                diag["tx_agree"] += 1
            else:
                diag["tx_disagree"] += 1
            diag["joined_mane" if mane_acc else "joined_no_mane"] += 1

    if diag["eligible"]:
        print("# AM-eligible missense (SNV, CSQ missense_variant): %d; joined %d "
              "(%.2f%%), unjoined %d"
              % (diag["eligible"], diag["eligible_joined"],
                 100.0 * diag["eligible_joined"] / diag["eligible"],
                 diag["unjoined_eligible"]), file=log)
        print("# unjoined for a KNOWN reason: %d in-frame indel / non-SNV, "
              "%d whose most severe CSQ is not missense_variant"
              % (diag["unjoined_indel"], diag["unjoined_not_missense_csq"]), file=log)
    if anno_ctx is not None and (diag["tx_agree"] + diag["tx_disagree"]):
        tot = diag["tx_agree"] + diag["tx_disagree"] + diag["tx_unknown"]
        print("# transcript concordance: AlphaMissense scored the SAME transcript "
              "the annotation used in %d of %d joined rows (%.2f%%); %d disagree"
              % (diag["tx_agree"], tot, 100.0 * diag["tx_agree"] / tot,
                 diag["tx_disagree"]), file=log)
        print("# of joined rows, %d sit on a transcript with a MANE Select "
              "accession and %d do not (canonical fallback)"
              % (diag["joined_mane"], diag["joined_no_mane"]), file=log)
    counts_diag = diag

    for k in pick_stats:
        pick_stats[k] = 0
    excluded_by_gene.clear()
    no_scored = []  # (gene, missense count) with no scored missense variant

    # --- per-gene score vectors
    drop = set(drop_anno)
    counts = {"genes": len(order), "variants": 0, "lof": 0, "missense": 0,
              "dropped_anno": 0, "emptied_regions": 0,
              "missense_scored": 0, "missense_imputed": 0, "missense_dropped": 0,
              "other": 0, "nonmissense_dropped": 0,
              "unparsed_ids": unparsed, "am_keys_matched": len(hits),
              "am_keys_wanted": len(wanted), "am_collisions": collisions,
              "from_isoform": from_isoform, "mane_hits": mane_hits}
    counts.update({"diag_" + k: v for k, v in counts_diag.items()})

    if nonmissense in ("drop", "gene-mean") and anno_ctx is None:
        raise SystemExit("--nonmissense %s needs --anno-table: the consequence lives there" % nonmissense)
    nonmis_csq, nonmis_kept = {}, {}   # dropped / kept (gene-mean) non-missense entries, by consequence
    pct_mode = am_scale == "percentile"
    counts.update({"nonmissense_matched": 0, "nonmissense_filled": 0, "missense_no_scale": 0})
    prov = open(provenance, "w") if provenance else None
    if prov is not None:
        # ANNOTATION/CSQ/ANNO_TRANSCRIPT: the bin, the most severe VEP consequence
        # and the transcript it was called on -- whether a missense-bin variant
        # is missense at all; SOURCE/TRANSCRIPT/AM: where its score came from;
        # with --am-scale percentile also REF_TRANSCRIPT (the gene's scale) and SCORE (the value written)
        prov.write("GENE\tID\tANNOTATION\tCSQ\tANNO_TRANSCRIPT\tSOURCE\tTRANSCRIPT\tAM"
                   + ("\tREF_TRANSCRIPT\tSCORE" if pct_mode else "") + "\n")
    out_lines = []
    for gene in order:
        r = regions[gene]
        raw_scores, kinds, pbuf = [], [], []   # pbuf: (index, provenance line without SCORE), written after the fill
        gref = ref_of.get(gene)
        for vid, anno in zip(r["var"], r["anno"]):
            counts["variants"] += 1
            csq = anno_ctx.get((gene, vid), (None, None, None, None))[2] if anno_ctx else None
            nonmis = anno in mis and anno_ctx is not None and not is_missense_csq(csq)
            if anno in lof:
                raw_scores.append(1.0)
                kinds.append("lof")
                counts["lof"] += 1
            elif nonmis and (nonmissense == "drop" or
                             (nonmissense == "gene-mean" and not is_protein_altering_csq(csq))):
                raw_scores.append(None)
                kinds.append("nonmissense")
                nonmis_csq[csq] = nonmis_csq.get(csq, 0) + 1
                if prov is not None:
                    actx = anno_ctx.get((gene, vid)) or ("NA",) * 3
                    pbuf.append((len(raw_scores) - 1, "%s\t%s\t%s\t%s\t%s\tdropped_nonmissense\t.\tNA" % (
                        gene, vid, anno, actx[2], actx[0])))
            elif anno in mis:
                counts["missense"] += 1
                key = parse_variant_id(vid, id_format)
                hit = pick(gene, vid, key)
                if nonmis and nonmissense == "gene-mean":
                    # a non-missense consequence takes AlphaMissense only from its OWN transcript's row
                    if hit is not None and hit[3] != "matched":
                        hit = None
                    counts["nonmissense_matched" if hit is not None else "nonmissense_filled"] += 1
                    nonmis_kept[csq] = nonmis_kept.get(csq, 0) + 1
                if pct_mode and hit is not None and gref is None:
                    counts["missense_no_scale"] += 1   # a row, but no transcript to rank it on: --missing
                    hit = None
                if prov is not None:
                    actx = (anno_ctx.get((gene, vid)) if anno_ctx else None) or ("NA",) * 3
                    pbuf.append((len(raw_scores), "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
                        gene, vid, anno, actx[2], actx[0],
                        "missing" if hit is None else hit[3],
                        "." if hit is None else hit[1],
                        "NA" if hit is None else "%.6f" % hit[0])))
                if hit is None:
                    raw_scores.append(None)
                    kinds.append("missense_missing")
                else:
                    raw_scores.append(percentile_in(hit[0], gref[2]) if pct_mode else hit[0])
                    kinds.append("missense")
                    counts["missense_scored"] += 1
            else:
                raw_scores.append(float(other))
                kinds.append("other")
                counts["other"] += 1

        # resolve missing missense
        scored = [s for s, k in zip(raw_scores, kinds) if k == "missense"]
        n_mis = sum(1 for k in kinds if k in ("missense", "missense_missing"))
        if n_mis and not scored:
            no_scored.append((gene, n_mis))
        if missing == "gene-mean":
            # a gene with no scored missense variant has no mean: --empty-gene-fill (percentile scale: 0.5, the
            # expected within-gene percentile, so FlexRV still tests its missense variants; raw: 0.0, as released)
            fill = sum(scored) / len(scored) if scored else empty_fill
        elif missing == "drop":
            fill = None
        else:
            fill = float(missing)

        keep = [True] * len(raw_scores)
        for i, k in enumerate(kinds):
            if k == "nonmissense":
                keep[i] = False
                counts["nonmissense_dropped"] += 1
                continue
            if k != "missense_missing":
                continue
            if fill is None:
                keep[i] = False
                counts["missense_dropped"] += 1
            else:
                raw_scores[i] = fill
                counts["missense_imputed"] += 1

        if prov is not None:
            for i, line in pbuf:
                if pct_mode:
                    s = raw_scores[i] if keep[i] else None
                    line += "\t%s\t%s" % (gref[1] if gref else ".", "NA" if s is None else "%.6f" % s)
                prov.write(line + "\n")

        if drop:
            for i, a in enumerate(r["anno"]):
                if a in drop and keep[i]:
                    keep[i] = False
                    counts["dropped_anno"] += 1

        final_var = [v for v, k in zip(r["var"], keep) if k]
        final_anno = [a for a, k in zip(r["anno"], keep) if k]
        final_score = [s for s, k in zip(raw_scores, keep) if k]
        if not final_var:
            counts["emptied_regions"] += 1
            print("# region %s has no variants left after dropping %s and is "
                  "omitted (an empty var line is not a valid region)"
                  % (gene, "/".join(sorted(drop))), file=log)
            continue

        for s in final_score:
            if s is None or not (0.0 <= s <= 1.0):
                raise SystemExit("region %s produced an out-of-range score %r; "
                                 "checkGroupFile would reject it" % (gene, s))

        # Emit: var, anno, score, then any weight/other lines, in that order --
        # the score line's position is FIXED immediately after anno.
        for kind, name, values, raw in r["lines"]:
            if kind == "var":
                out_lines.append("%s var %s" % (gene, " ".join(final_var)))
                out_lines.append("%s anno %s" % (gene, " ".join(final_anno)))
                # percentile scores keep 6 decimals: in the largest genes (N > 30,000 possible missense) one
                # rank step is < 1e-4, so 4 decimals would write the bottom variant as 0 and tie neighbours
                # (2026-10-21); raw mode keeps 4, so files built before stay byte-identical
                out_lines.append("%s score:%s %s" % (
                    gene, score_name, " ".join(("%.6f" if pct_mode else "%.4f") % s for s in final_score)))
            elif kind == "anno":
                continue  # written above, so the three stay adjacent and ordered
            elif kind == "weight":
                kept = [v for v, k in zip(values, keep) if k]
                out_lines.append("%s %s %s" % (
                    gene, raw.split()[1], " ".join(kept)))
            else:
                out_lines.append(raw)

    if prov is not None:
        prov.close()
    if nonmissense == "gene-mean":
        print("# --nonmissense gene-mean: %d protein-altering missense-bin variants whose consequence is not "
              "missense_variant (%s): %d scored from their own transcript's AlphaMissense row, %d at the gene mean; "
              "%d with no protein change dropped (%s)"
              % (counts["nonmissense_matched"] + counts["nonmissense_filled"], ", ".join(
                  "%s %d" % (k, v) for k, v in sorted(nonmis_kept.items(), key=lambda x: -x[1])),
                 counts["nonmissense_matched"], counts["nonmissense_filled"], counts["nonmissense_dropped"],
                 ", ".join("%s %d" % (k, v) for k, v in sorted(nonmis_csq.items(), key=lambda x: -x[1]))), file=log)
    if counts["missense_no_scale"]:
        print("# --am-scale percentile: %d missense rows had an AlphaMissense row but their gene has no "
              "AlphaMissense transcript to rank on; they take --missing" % counts["missense_no_scale"], file=log)
    if nonmissense == "drop":  # (gene-mean reports its drops above)
        print("# --nonmissense drop: removed %d missense-bin variants whose consequence is not "
              "missense_variant (%s)" % (counts["nonmissense_dropped"], ", ".join(
                  "%s %d" % (k, v) for k, v in sorted(nonmis_csq.items(), key=lambda x: -x[1]))),
              file=log)
    if prefer_tx:
        ps = pick_stats
        fell = (ps["fallback_canonical"] + ps["fallback_isoform"] + ps["fallback_max"]
                + ps["missing_excluded"])
        tot = ps["by_transcript"] + fell
        print("# transcript-matched selection: %d of %d missense rows with an "
              "AlphaMissense row took the row AlphaMissense computed on the "
              "ANNOTATION's own transcript (%.2f%%); %d had no row on it"
              % (ps["by_transcript"], tot,
                 100.0 * ps["by_transcript"] / max(1, tot), fell), file=log)
        if fallback == "max":
            print("# fallback (LEGACY max over every row at the coordinate): %d"
                  % ps["fallback_max"], file=log)
        else:
            print("# fallback%s: %d fallback-canonical (%d averaged over >1 transcript "
                  "of one UniProt entry), %d fallback-same-gene-isoform "
                  "(%d averaged over >1 isoform), %d missing (only rows another gene "
                  "uses, or canonical rows of >1 UniProt entry -> --missing); plus %d with no "
                  "AlphaMissense row at all"
                  % ("" if chrom_filter is None else " on chr%s" % chrom_filter,
                     ps["fallback_canonical"], ps["fallback_canonical_averaged"],
                     ps["fallback_isoform"],
                     ps["fallback_isoform_averaged"], ps["missing_excluded"],
                     ps["missing_no_row"]), file=log)
            if mane_genes is None and gtf_genes is None:
                print("# (no --mane or --gtf: a row is excluded only when the "
                      "annotation assigns its transcript to another gene at that "
                      "coordinate)", file=log)
        lost = [(g, n) for g, n in no_scored if excluded_by_gene.get(g)]
        if no_scored:
            print("# %d regions have NO scored missense variant (--missing %s fills "
                  "them%s); %d of these only because every AlphaMissense row "
                  "there is another gene's%s"
                  % (len(no_scored), missing,
                     " with %g" % empty_fill if missing == "gene-mean" else "", len(lost),
                     (": " + ", ".join("%s (%d of %d missense)" % (
                         g, excluded_by_gene[g], n) for g, n in lost)) if lost else ""),
                  file=log)

    with open(group_out, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    return counts


def _build_from_column(regions, order, lof, mis, drop_anno, group_out,
                       score_name, anno_ctx, score_column, score_transform,
                       missing, other, counts, zero_missing, log, nonmissense="impute"):
    """Score every region from an annotation-table COLUMN.

    Same construction as the AlphaMissense path -- LoF 1.0 by annotation, the
    score classes from the source, everything else --other -- but the value is
    read rather than joined, and the per-gene transform runs before the
    missing-value policy so that an imputed value is on the TRANSFORMED scale
    (imputing a raw CADD PHRED and then ranking it would let the fill change
    every other variant's rank).

    `nonmissense` as in build(): "impute" (default) treats every missense-bin entry alike; "drop" removes
    those whose CSQ is not missense_variant; "gene-mean" removes only those with no protein change and fills
    the protein-altering ones (in-frame indels, LOFTEE-LC LoF, start/stop lost) with the gene mean unless the
    column has a value for them -- so a second weight set (e.g. scallion, 2026-10-21) tests exactly the
    variants the AlphaMissense percentile file does.
    """
    drop = set(drop_anno)
    transform = score_transform or SCORE_COLUMNS.get(score_column, ("none",))[0]
    counts = dict(counts)
    counts.update({"variants": 0, "lof": 0, "missense": 0, "missense_scored": 0,
                   "missense_imputed": 0, "missense_dropped": 0, "other": 0,
                   "dropped_anno": 0, "emptied_regions": 0, "unparsed_ids": 0,
                   "am_keys_matched": 0, "am_keys_wanted": 0, "am_collisions": 0,
                   "from_isoform": 0, "mane_hits": None, "nonmissense_dropped": 0,
                   "nonmissense_filled": 0, "nonmissense_matched": 0})
    print("# column mode: %s, 0 treated as %s, transform '%s'%s"
          % (score_column, "MISSING" if zero_missing else "a real value", transform,
             ("  [" + SCORE_COLUMNS[score_column][2] + "]")
             if score_column in SCORE_COLUMNS else ""), file=log)

    out_lines = []
    for gene in order:
        r = regions[gene]
        raw, kinds = [None] * len(r["var"]), []
        pending = []
        for i, (vid, anno) in enumerate(zip(r["var"], r["anno"])):
            counts["variants"] += 1
            if anno in lof:
                raw[i] = 1.0
                kinds.append("lof")
                counts["lof"] += 1
            elif anno in mis and nonmissense != "impute" and not is_missense_csq(
                    anno_ctx.get((gene, vid), (None, None, None, None))[2]) and (
                    nonmissense == "drop" or not is_protein_altering_csq(
                        anno_ctx.get((gene, vid), (None, None, None, None))[2])):
                kinds.append("nonmissense")
            elif anno in mis:
                counts["missense"] += 1
                ctx = anno_ctx.get((gene, vid), (None, None, None, None))
                val = ctx[3]
                if nonmissense == "gene-mean" and not is_missense_csq(ctx[2]):
                    counts["nonmissense_matched" if val is not None else "nonmissense_filled"] += 1
                if val is None:
                    kinds.append("missense_missing")
                else:
                    kinds.append("missense")
                    pending.append((i, val))
            else:
                raw[i] = float(other)
                kinds.append("other")
                counts["other"] += 1

        for i, v in apply_transform(pending, transform).items():
            raw[i] = v
            counts["missense_scored"] += 1

        scored = [raw[i] for i, k in enumerate(kinds) if k == "missense"]
        if missing == "gene-mean":
            fill = sum(scored) / len(scored) if scored else 0.0
        elif missing == "drop":
            fill = None
        else:
            fill = float(missing)

        keep = [True] * len(raw)
        for i, k in enumerate(kinds):
            if k == "nonmissense":
                keep[i] = False
                counts["nonmissense_dropped"] += 1
                continue
            if k != "missense_missing":
                continue
            if fill is None:
                keep[i] = False
                counts["missense_dropped"] += 1
            else:
                raw[i] = fill
                counts["missense_imputed"] += 1
        if drop:
            for i, a in enumerate(r["anno"]):
                if a in drop and keep[i]:
                    keep[i] = False
                    counts["dropped_anno"] += 1

        fv = [v for v, k in zip(r["var"], keep) if k]
        fa = [a for a, k in zip(r["anno"], keep) if k]
        fs = [v for v, k in zip(raw, keep) if k]
        if not fv:
            counts["emptied_regions"] += 1
            print("# region %s has no variants left after dropping %s and is "
                  "omitted" % (gene, "/".join(sorted(drop))), file=log)
            continue
        for v in fs:
            if v is None or not (0.0 <= v <= 1.0):
                raise SystemExit("region %s produced an out-of-range score %r" % (gene, v))

        for kind, name, values, rawline in r["lines"]:
            if kind == "var":
                out_lines.append("%s var %s" % (gene, " ".join(fv)))
                out_lines.append("%s anno %s" % (gene, " ".join(fa)))
                out_lines.append("%s score:%s %s" % (     # 6 decimals: column scores are often percentiles too
                    gene, score_name, " ".join("%.6f" % v for v in fs)))
            elif kind == "anno":
                continue
            elif kind == "weight":
                out_lines.append("%s %s %s" % (
                    gene, rawline.split()[1],
                    " ".join(v for v, k in zip(values, keep) if k)))
            else:
                out_lines.append(rawline)

    with open(group_out, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    return counts


# ---------------------------------------------------------------- self-test

def _selftest():
    """Positive and negative controls on a planted AlphaMissense table.

    Positive: scores planted for known variants must come back at exactly the
    right offsets. Negative: an AlphaMissense table whose coordinates match
    NOTHING must produce zero joins and must not invent a score -- every
    missense value must be the imputation, not a pathogenicity. A third arm
    shuffles the planted values across keys and asserts the output changes,
    which is what catches a join that silently ignores its key.
    """
    tmp = tempfile.mkdtemp(prefix="flexrv_am_selftest_")
    rng = random.Random(20260919)
    ok = True

    def check(label, cond):
        nonlocal ok
        print("  %-58s %s" % (label, "PASS" if cond else "FAIL"))
        if not cond:
            ok = False

    # --- a two-gene group file: 10 lof, 20 missense, 10 synonymous per gene
    genes, truth = [], {}
    gf = os.path.join(tmp, "group.txt")
    with open(gf, "w") as fh:
        pos = 1000
        for g in ("GENEA", "GENEB"):
            vids, annos = [], []
            for i in range(40):
                pos += rng.randint(3, 50)
                vids.append("20:%d:A:G" % pos)
                annos.append("pLoF" if i < 10 else
                             ("missense" if i < 30 else "synonymous"))
            genes.append((g, vids, annos))
            fh.write("%s var %s\n" % (g, " ".join(vids)))
            fh.write("%s anno %s\n" % (g, " ".join(annos)))

    # --- plant AlphaMissense rows for 14 of each gene's 20 missense variants
    am = os.path.join(tmp, "am.tsv.gz")
    planted = {}
    with gzip.open(am, "wt") as fh:
        fh.write("# copied header shape from the AlphaMissense release\n")
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        for g, vids, annos in genes:
            mis = [v for v, a in zip(vids, annos) if a == "missense"]
            for v in mis[:14]:
                chrom, p, ref, alt = parse_variant_id(v, "colon")
                s = round(rng.random(), 4)
                planted[v] = s
                fh.write("\t".join(["chr20", p, ref, alt, "hg38", "P00000",
                                    "ENST00000001.1", "A1B", "%.4f" % s,
                                    "ambiguous"]) + "\n")
        # decoys: real-looking rows at coordinates the group file does not hold
        for k in range(50):
            fh.write("\t".join(["chr20", str(9_000_000 + k), "C", "T", "hg38",
                                "P00000", "ENST00000009.1", "A1B", "0.9999",
                                "likely_pathogenic"]) + "\n")

    def run(am_file, out, missing="gene-mean"):
        return build(gf, out, am_file, ["pLoF"], ["missense"],
                     score_name="AM", missing=missing, other=0.0,
                     log=open(os.devnull, "w"))

    # ---- POSITIVE CONTROL
    out1 = os.path.join(tmp, "pos.txt")
    c1 = run(am, out1)
    got = {}
    for gene, kind, name, values, raw in read_group_regions(out1):
        if kind == "var":
            cur_var = values
        elif kind == "score":
            got[gene] = dict(zip(cur_var, [float(v) for v in values]))
    check("positive: 28 of 40 missense joined", c1["missense_scored"] == 28)
    check("positive: 12 missense imputed", c1["missense_imputed"] == 12)
    check("positive: 20 lof scored 1.0",
          all(got[g][v] == 1.0 for g, vids, annos in genes
              for v, a in zip(vids, annos) if a == "pLoF"))
    check("positive: planted scores recovered at the right IDs",
          all(abs(got[g][v] - planted[v]) < 5e-5 for g, vids, annos in genes
              for v in vids if v in planted))
    check("positive: synonymous scored 0.0 (--other)",
          all(got[g][v] == 0.0 for g, vids, annos in genes
              for v, a in zip(vids, annos) if a == "synonymous"))
    per_gene_fill = {}
    for g, vids, annos in genes:
        s = [planted[v] for v, a in zip(vids, annos) if a == "missense" and v in planted]
        per_gene_fill[g] = sum(s) / len(s)
    check("positive: missing missense == that GENE's mean, not the pooled mean",
          all(abs(got[g][v] - per_gene_fill[g]) < 5e-5
              for g, vids, annos in genes
              for v, a in zip(vids, annos)
              if a == "missense" and v not in planted)
          and abs(per_gene_fill["GENEA"] - per_gene_fill["GENEB"]) > 1e-3)
    check("positive: every value in [0, 1]",
          all(0.0 <= v <= 1.0 for d in got.values() for v in d.values()))
    check("positive: score line length == var line length",
          all(len(got[g]) == 40 for g, _, _ in genes))

    # ---- NEGATIVE CONTROL: an AlphaMissense file matching nothing
    am_null = os.path.join(tmp, "am_null.tsv.gz")
    with gzip.open(am_null, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        for k in range(200):
            fh.write("\t".join(["chr20", str(50_000_000 + k), "C", "T", "hg38",
                                "P00000", "ENST00000009.1", "A1B", "0.8888",
                                "likely_pathogenic"]) + "\n")
    out2 = os.path.join(tmp, "neg.txt")
    c2 = run(am_null, out2, missing="0")
    got2 = {}
    for gene, kind, name, values, raw in read_group_regions(out2):
        if kind == "var":
            cur_var = values
        elif kind == "score":
            got2[gene] = dict(zip(cur_var, [float(v) for v in values]))
    check("negative: zero keys matched", c2["missense_scored"] == 0)
    check("negative: no score of 0.8888 appears anywhere (nothing invented)",
          not any(abs(v - 0.8888) < 1e-9 for d in got2.values() for v in d.values()))
    check("negative: all 40 missense take the --missing value",
          c2["missense_imputed"] == 40 and
          all(got2[g][v] == 0.0 for g, vids, annos in genes
              for v, a in zip(vids, annos) if a == "missense"))
    check("negative: lof still 1.0 (annotation-based, not from AlphaMissense)",
          all(got2[g][v] == 1.0 for g, vids, annos in genes
              for v, a in zip(vids, annos) if a == "pLoF"))

    # ---- SHUFFLE CONTROL: same values, wrong keys -> different output
    am_shuf = os.path.join(tmp, "am_shuf.tsv.gz")
    keys = [parse_variant_id(v, "colon") for v in planted]
    vals = list(planted.values())
    rng.shuffle(vals)
    with gzip.open(am_shuf, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        for (c, p, r, a), s in zip(keys, vals):
            fh.write("\t".join(["chr20", p, r, a, "hg38", "P00000",
                                "ENST00000001.1", "A1B", "%.4f" % s,
                                "ambiguous"]) + "\n")
    out3 = os.path.join(tmp, "shuf.txt")
    run(am_shuf, out3)
    check("shuffle: same values at permuted keys change the output "
          "(the join reads its key)",
          open(out1).read() != open(out3).read())

    # ---- DROP policy
    out4 = os.path.join(tmp, "drop.txt")
    c4 = run(am, out4, missing="drop")
    lens = {}
    for gene, kind, name, values, raw in read_group_regions(out4):
        lens.setdefault(gene, {})[kind] = len(values)
    check("drop: 12 unscored missense removed from var/anno/score together",
          c4["missense_dropped"] == 12 and
          all(v["var"] == v["anno"] == v["score"] == 34 for v in lens.values()))

    # ---- structural: score line sits immediately after anno, before weights
    seq = [(g, k) for g, k, _, _, _ in read_group_regions(out1)
           if k in ("var", "anno", "score")]
    check("structure: every region is var, anno, score in that order",
          seq == [(g, k) for g, _, _ in genes for k in ("var", "anno", "score")])

    # ---- TRANSCRIPT-PREFERENCE CONTROL
    # The annotation names transcript T_MANE. AlphaMissense's CANONICAL file
    # scores a DIFFERENT transcript (T_OTHER) at the same coordinates, and the
    # ISOFORMS file holds T_MANE. Whether the tool reads the transcript, rather
    # than just the coordinates, is the whole question -- so the two files are
    # given deliberately different scores and the output must show which was
    # taken. Without --prefer-anno-transcript it must take the canonical one;
    # with it (and --isoforms) it must take the annotation's transcript.
    anno_tsv = os.path.join(tmp, "anno_table.tsv")
    mis_ids = [v for g, vids, annos in genes
               for v, a in zip(vids, annos) if a == "missense"]
    with open(anno_tsv, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED",
                            "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL",
                            "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        for g, vids, annos in genes:
            for v, a in zip(vids, annos):
                csq = {"pLoF": "stop_gained", "missense": "missense_variant",
                       "synonymous": "synonymous_variant"}[a]
                fh.write("\t".join([v, g, "NA", "0", "0", csq, "ENST00000MANE",
                                    "NM_000001.1", "YES", "protein_coding",
                                    "0", a]) + "\n")

    CANON_S, MANE_S = 0.1111, 0.9999
    am_canon_other = os.path.join(tmp, "am_canon_other.tsv.gz")
    with gzip.open(am_canon_other, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        for v in mis_ids:
            c, pp, r, aa = parse_variant_id(v, "colon")
            fh.write("\t".join(["chr20", pp, r, aa, "hg38", "P0", "ENST00000OTHER.2",
                                "A1B", "%.4f" % CANON_S, "ambiguous"]) + "\n")
    am_iso_mane = os.path.join(tmp, "am_iso_mane.tsv.gz")
    with gzip.open(am_iso_mane, "wt") as fh:
        fh.write("\t".join(["#CHROM", "POS", "REF", "ALT", "genome",
                            "transcript_id", "protein_variant",
                            "am_pathogenicity", "am_class"]) + "\n")
        for v in mis_ids:
            c, pp, r, aa = parse_variant_id(v, "colon")
            fh.write("\t".join(["chr20", pp, r, aa, "hg38", "ENST00000MANE.3",
                                "A1B", "%.4f" % MANE_S, "likely_pathogenic"]) + "\n")

    def scores_of(path):
        out = {}
        for gene, kind, name, values, raw in read_group_regions(path):
            if kind == "var":
                cv = values
            elif kind == "score":
                out.update(dict(zip(cv, [float(x) for x in values])))
        return out

    out5 = os.path.join(tmp, "tx_canonical.txt")
    build(gf, out5, am_canon_other, ["pLoF"], ["missense"], score_name="AM",
          anno_table_path=anno_tsv, log=open(os.devnull, "w"))
    s5 = scores_of(out5)
    check("transcript: WITHOUT the flag, the canonical row wins",
          all(abs(s5[v] - CANON_S) < 5e-5 for v in mis_ids))

    out6 = os.path.join(tmp, "tx_mane.txt")
    build(gf, out6, am_canon_other, ["pLoF"], ["missense"], score_name="AM",
          anno_table_path=anno_tsv, isoforms_path=am_iso_mane,
          prefer_anno_transcript=True, log=open(os.devnull, "w"))
    s6 = scores_of(out6)
    check("transcript: WITH the flag, the annotation's transcript wins",
          all(abs(s6[v] - MANE_S) < 5e-5 for v in mis_ids))
    check("transcript: the two modes actually differ (the control is live)",
          abs(CANON_S - MANE_S) > 0.5 and
          any(abs(s5[v] - s6[v]) > 0.5 for v in mis_ids))

    out7 = os.path.join(tmp, "tx_nofall.txt")
    build(gf, out7, am_canon_other, ["pLoF"], ["missense"], score_name="AM",
          anno_table_path=anno_tsv, prefer_anno_transcript=True,
          log=open(os.devnull, "w"))
    s7 = scores_of(out7)
    check("transcript: with the flag but NO isoforms file, it falls back to "
          "the other transcript rather than dropping the score",
          all(abs(s7[v] - CANON_S) < 5e-5 for v in mis_ids))

    # ---- DROP-ANNO CONTROL
    # Removing an annotation class must remove exactly those variants and
    # leave every surviving score bit-identical. The second half is the real
    # assertion: gene-mean imputation is computed over scored MISSENSE
    # variants, so dropping a non-missense class must not perturb it -- but
    # that is a property of the code, not a law, so it is checked rather than
    # argued.
    out8 = os.path.join(tmp, "dropped.txt")
    c8 = build(gf, out8, am, ["pLoF"], ["missense"], score_name="AM",
               drop_anno=["synonymous"], log=open(os.devnull, "w"))
    kept = {}
    for gene, kind, name, values, raw in read_group_regions(out8):
        if kind == "var":
            cv = values
        elif kind == "anno":
            ca = values
        elif kind == "score":
            for v, a, sc in zip(cv, ca, values):
                kept[(gene, v)] = (a, float(sc))
    check("drop-anno: 20 synonymous removed, 60 variants remain",
          c8["dropped_anno"] == 20 and len(kept) == 60)
    check("drop-anno: no synonymous survives",
          not any(a == "synonymous" for a, _ in kept.values()))
    check("drop-anno: every SURVIVING score is bit-identical to the undropped run",
          all(kept[(g, v)][1] == got[g][v] for (g, v) in kept))
    check("drop-anno: no region emptied here", c8["emptied_regions"] == 0)

    # a region made entirely of the dropped class must be omitted, not emitted
    # with an empty var line (which checkGroupFile would reject)
    gf2 = os.path.join(tmp, "group_onegene_allsyn.txt")
    with open(gf2, "w") as fh:
        for g, vids, annos in genes:
            fh.write("%s var %s\n" % (g, " ".join(vids)))
            fh.write("%s anno %s\n" % (g, " ".join(annos)))
        fh.write("GENEC var 20:9999991:A:G 20:9999992:A:G\n")
        fh.write("GENEC anno synonymous synonymous\n")
    out9 = os.path.join(tmp, "dropped_empty.txt")
    c9 = build(gf2, out9, am, ["pLoF"], ["missense"], score_name="AM",
               drop_anno=["synonymous"], log=open(os.devnull, "w"))
    seen = set(g for g, k, _, _, _ in read_group_regions(out9) if k == "var")
    check("drop-anno: an all-dropped region is omitted, not left empty",
          c9["emptied_regions"] == 1 and "GENEC" not in seen and len(seen) == 2)

    # ---- TRANSFORM KERNELS, pinned to closed forms
    pc = apply_transform([(0, 5.0), (1, 1.0), (2, 3.0), (3, 9.0)],
                         "percentile-in-gene")
    check("transform: percentile-in-gene ranks correctly",
          pc == {1: 0.25, 2: 0.5, 0: 0.75, 3: 1.0})
    pt = apply_transform([(0, 2.0), (1, 2.0), (2, 9.0)], "percentile-in-gene")
    check("transform: ties take the average rank (1.5/3, 1.5/3, 3/3)",
          abs(pt[0] - 0.5) < 1e-12 and pt[0] == pt[1] and pt[2] == 1.0)
    cp = apply_transform([(0, 10.0), (1, 20.0), (2, 30.0)], "cadd-phred")
    check("transform: cadd-phred is 1 - 10^(-p/10) (10->0.9, 20->0.99, 30->0.999)",
          all(abs(cp[i] - v) < 1e-12
              for i, v in [(0, 0.9), (1, 0.99), (2, 0.999)]))
    mm = apply_transform([(0, 0.0), (1, 5.0), (2, 10.0), (3, 20.0)], "minmax:0,10")
    check("transform: minmax clips to [0,1]",
          mm == {0: 0.0, 1: 0.5, 2: 1.0, 3: 1.0})
    try:
        apply_transform([(0, 7.0)], "none"); bad = False
    except SystemExit:
        bad = True
    check("transform: 'none' REFUSES a value outside [0,1] rather than emitting it",
          bad)

    # ---- COLUMN MODE, with 0 meaning absent
    col_tsv = os.path.join(tmp, "col_table.tsv")
    with open(col_tsv, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED",
                            "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL",
                            "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        planted_col = {}
        for gi, (g, vids, annos) in enumerate(genes):
            k = 0
            for v, a in zip(vids, annos):
                rev = "0"
                if a == "missense":
                    k += 1
                    # every 5th missense has NO revel (stored as 0 = absent)
                    rev = "0" if k % 5 == 0 else "%.4f" % (0.05 * k + 0.01 * gi)
                    if rev != "0":
                        planted_col[(g, v)] = float(rev)
                fh.write("\t".join([v, g, "NA", rev, "0", "missense_variant",
                                    "ENST1", "NM_1", "YES", "protein_coding",
                                    "0", a]) + "\n")
    out10 = os.path.join(tmp, "col_revel.txt")
    c10 = build(gf, out10, None, ["pLoF"], ["missense"], score_name="REVEL",
                anno_table_path=col_tsv, score_column="REVEL_SCORE",
                log=open(os.devnull, "w"))
    g10 = {}
    for gene, kind, name, values, raw in read_group_regions(out10):
        if kind == "var":
            cv = values
        elif kind == "score":
            g10[gene] = dict(zip(cv, [float(x) for x in values]))
            sname = name
    check("column: the score line is named from --name", sname == "REVEL")
    check("column: 32 of 40 missense scored, 8 absent (0) imputed",
          c10["missense_scored"] == 32 and c10["missense_imputed"] == 8)
    check("column: planted REVEL values land verbatim (transform 'none')",
          all(abs(g10[g][v] - planted_col[(g, v)]) < 5e-5
              for (g, v) in planted_col))
    check("column: a 0 in the column is NOT emitted as a 0 score",
          all(g10[g][v] != 0.0 for g, vids, annos in genes
              for v, a in zip(vids, annos)
              if a == "missense" and (g, v) not in planted_col))
    check("column: lof still 1.0, synonymous still --other",
          all(g10[g][v] == 1.0 for g, vids, annos in genes
              for v, a in zip(vids, annos) if a == "pLoF") and
          all(g10[g][v] == 0.0 for g, vids, annos in genes
              for v, a in zip(vids, annos) if a == "synonymous"))

    # percentile-in-gene must rank WITHIN a gene, not across the file
    out11 = os.path.join(tmp, "col_pct.txt")
    build(gf, out11, None, ["pLoF"], ["missense"], score_name="REVEL",
          anno_table_path=col_tsv, score_column="REVEL_SCORE",
          score_transform="percentile-in-gene", log=open(os.devnull, "w"))
    g11 = {}
    for gene, kind, name, values, raw in read_group_regions(out11):
        if kind == "var":
            cv = values
        elif kind == "score":
            g11[gene] = dict(zip(cv, [float(x) for x in values]))
    tops = [max(g11[g][v] for g, vids, annos in [(gg, vv, aa)]
                for v, a in zip(vids, annos)
                if a == "missense" and (gg, v) in planted_col)
            for gg, vv, aa in genes]
    check("column: percentile-in-gene gives EACH gene its own top rank of 1.0",
          len(tops) == 2 and all(abs(t - 1.0) < 1e-9 for t in tops))

    # ---- FALLBACK CONTROL: no row on the annotation's transcript
    # Two genes share coordinates. GENEP's annotation transcript TP is absent
    # from AlphaMissense; GENEQ's TQ is present and scores HIGH. The legacy max
    # fallback hands GENEP GENEQ's score; the fallback must not.
    #   S1 shared: canonical TPC 0.30, TQ 0.95  -> P: 0.30 (canonical)
    #   S2 shared: only TQ 0.90                 -> P: missing (gene mean)
    #   S3 P only: TP 0.20 (iso), TPC 0.70      -> P: 0.20 (matched beats canonical)
    #   S4 shared: TP1 0.40, TP2 0.60 (iso), TQ 0.99 -> P: 0.50 (isoform mean)
    #   S5 P only: TR 0.80 (iso), MANE puts TR in GENER -> with --mane: missing
    #   S6 P only: canonical TX 0.10 (UniProt P1), TY 0.85 (P2) -> ambiguous: missing
    #   S7 P only: canonical TZ1 0.20, TZ2 0.60, both P3 -> one protein: 0.40
    #   M1..M3 P only: matched TP 0.11/0.22/0.33
    #   S8 GENEZ only (annotated on TZZ, absent): the one row is TQC 0.90, a
    #      canonical transcript of GENEQ that no annotation uses -- the
    #      readthrough case. Only --gtf (TQC -> GENEQ) can exclude it.
    fx = {"S1": 5001, "S2": 5002, "S3": 5003, "S4": 5004, "S5": 5005,
          "S6": 5006, "S7": 5007, "S8": 5008, "M1": 5011, "M2": 5012, "M3": 5013}
    vid = dict((k, "20:%d:A:G" % p) for k, p in fx.items())
    p_ids = [vid[k] for k in ("S1", "S2", "S3", "S4", "S5", "S6", "S7",
                              "M1", "M2", "M3")]
    q_ids = [vid[k] for k in ("S1", "S2", "S4")]
    gf3 = os.path.join(tmp, "group_shared.txt")
    with open(gf3, "w") as fh:
        for g, ids in (("GENEP", p_ids), ("GENEQ", q_ids), ("GENEZ", [vid["S8"]])):
            lof_id = "20:%d:A:T" % {"GENEP": 5090, "GENEQ": 5091, "GENEZ": 5092}[g]
            fh.write("%s var %s %s\n" % (g, " ".join(ids), lof_id))
            fh.write("%s anno %s pLoF\n" % (g, " ".join(["missense"] * len(ids))))
    anno3 = os.path.join(tmp, "anno_shared.tsv")
    with open(anno3, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED",
                            "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL",
                            "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        for g, ids, tx in (("GENEP", p_ids, "ENST0000000TP"),
                           ("GENEQ", q_ids, "ENST0000000TQ"),
                           ("GENEZ", [vid["S8"]], "ENST000000TZZ")):
            for v in ids:
                fh.write("\t".join([v, g, "NA", "0", "0", "missense_variant", tx,
                                    "NM_1", "YES", "protein_coding", "0",
                                    "missense"]) + "\n")
    canon_rows = [("S1", "ENST000000TPC.1", 0.30, "P0"), ("S3", "ENST000000TPC.1", 0.70, "P0"),
                  ("S6", "ENST0000000TX.1", 0.10, "P1"), ("S6", "ENST0000000TY.1", 0.85, "P2"),
                  ("S7", "ENST000000TZ1.1", 0.20, "P3"), ("S7", "ENST000000TZ2.1", 0.60, "P3-2"),
                  ("S8", "ENST000000TQC.1", 0.90, "P9")]
    iso_rows = [("S1", "ENST0000000TQ.2", 0.95), ("S2", "ENST0000000TQ.2", 0.90),
                ("S3", "ENST0000000TP.2", 0.20), ("S4", "ENST000000TP1.1", 0.40),
                ("S4", "ENST000000TP2.1", 0.60), ("S4", "ENST0000000TQ.2", 0.99),
                ("S5", "ENST0000000TR.1", 0.80), ("M1", "ENST0000000TP.2", 0.11),
                ("M2", "ENST0000000TP.2", 0.22), ("M3", "ENST0000000TP.2", 0.33)]
    am3c = os.path.join(tmp, "am_shared_canon.tsv.gz")
    with gzip.open(am3c, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        for k, tx, sc, up in canon_rows:
            fh.write("\t".join(["chr20", str(fx[k]), "A", "G", "hg38", up, tx,
                                "A1B", "%.4f" % sc, "ambiguous"]) + "\n")
    am3i = os.path.join(tmp, "am_shared_iso.tsv.gz")
    with gzip.open(am3i, "wt") as fh:
        fh.write("\t".join(["#CHROM", "POS", "REF", "ALT", "genome",
                            "transcript_id", "protein_variant",
                            "am_pathogenicity", "am_class"]) + "\n")
        for k, tx, sc in iso_rows:
            fh.write("\t".join(["chr20", str(fx[k]), "A", "G", "hg38", tx,
                                "A1B", "%.4f" % sc, "ambiguous"]) + "\n")
    mane3 = os.path.join(tmp, "mane_shared.tsv")
    with open(mane3, "w") as fh:
        fh.write("\t".join(["#NCBI_GeneID", "Ensembl_Gene", "HGNC_ID", "symbol",
                            "name", "RefSeq_nuc", "RefSeq_prot", "Ensembl_nuc",
                            "Ensembl_prot", "MANE_status"]) + "\n")
        fh.write("\t".join(["GeneID:9", "ENSG0000000R.3", "HGNC:9", "GENER", "r",
                            "NM_9.1", "NP_9.1", "ENST0000000TR.1", "ENSP9.1",
                            "MANE Select"]) + "\n")

    # GENCODE-style GTF of the annotation's release: TPC is GENEP's, TQC is
    # GENEQ's; TZ1/TZ2 are absent (retired by that release: no evidence).
    gtf3 = os.path.join(tmp, "gencode_shared.gtf")
    with open(gtf3, "w") as fh:
        fh.write("##description: selftest\n")
        for tx, g in (("ENST000000TPC.1", "GENEP.4"), ("ENST0000000TP.2", "GENEP.4"),
                      ("ENST000000TQC.1", "GENEQ.2"), ("ENST0000000TQ.2", "GENEQ.2"),
                      ("ENST000000TQX.1", "GENEQ.2")):
            fh.write("chr20\tHAVANA\tgene\t1\t9\t.\t+\t.\tgene_id \"%s\";\n" % g)
            fh.write("chr20\tHAVANA\ttranscript\t1\t9\t.\t+\t.\tgene_id \"%s\"; "
                     "transcript_id \"%s\"; gene_name \"X\";\n" % (g, tx))

    def per_gene(path):
        out = {}
        for gene, kind, name, values, raw in read_group_regions(path):
            if kind == "var":
                cv = values
            elif kind == "score":
                out[gene] = dict(zip(cv, [float(x) for x in values]))
        return out

    def run3(tag, **kw):
        o = os.path.join(tmp, "shared_%s.txt" % tag)
        c = build(gf3, o, am3c, ["pLoF"], ["missense"], score_name="AM",
                  anno_table_path=anno3, isoforms_path=am3i,
                  prefer_anno_transcript=True, log=open(os.devnull, "w"), **kw)
        return per_gene(o), c

    old3, _ = run3("max", fallback="max")
    new3, _ = run3("safe_nomane")
    newm3, _ = run3("safe_mane", mane_path=mane3)
    gtf_log = os.path.join(tmp, "shared_gtf.log")
    o = os.path.join(tmp, "shared_gtf.txt")
    with open(gtf_log, "w") as lg:
        build(gf3, o, am3c, ["pLoF"], ["missense"], score_name="AM",
              anno_table_path=anno3, isoforms_path=am3i, prefer_anno_transcript=True,
              mane_path=mane3, gtf_path=gtf3, log=lg)
    newg3 = per_gene(o)
    prov3 = os.path.join(tmp, "shared_gtf.prov.tsv")
    build(gf3, os.path.join(tmp, "shared_prov.txt"), am3c, ["pLoF"], ["missense"],
          score_name="AM", anno_table_path=anno3, isoforms_path=am3i,
          prefer_anno_transcript=True, mane_path=mane3, gtf_path=gtf3,
          provenance=prov3, log=open(os.devnull, "w"))
    prows = [l.rstrip("\n").split("\t") for l in open(prov3)]
    pz = [r for r in prows if r[0] == "GENEZ"]
    P, Q = "GENEP", "GENEQ"
    check("fallback: LEGACY rule takes GENEQ's 0.95 for GENEP at a shared site",
          abs(old3[P][vid["S1"]] - 0.95) < 5e-5 and abs(old3[P][vid["S2"]] - 0.90) < 5e-5
          and abs(old3[P][vid["S4"]] - 0.99) < 5e-5)
    check("fallback: new rule takes the canonical row (0.30), not the max",
          abs(new3[P][vid["S1"]] - 0.30) < 5e-5)
    check("fallback: matched transcript still beats canonical (0.20, not 0.70)",
          abs(new3[P][vid["S3"]] - 0.20) < 5e-5 and abs(old3[P][vid["S3"]] - 0.20) < 5e-5)
    check("fallback: same-gene isoforms averaged (0.50), other gene's 0.99 excluded",
          abs(new3[P][vid["S4"]] - 0.50) < 5e-5)
    # with --mane: scored P rows are S1 .30, S3 .20, S4 .50, S7 .40, M .11 .22 .33
    gm = (0.30 + 0.20 + 0.50 + 0.40 + 0.11 + 0.22 + 0.33) / 7
    check("fallback: only another gene's row -> --missing (P's gene mean)",
          abs(newm3[P][vid["S2"]] - gm) < 5e-5)
    check("fallback: canonical rows of ONE UniProt entry averaged (0.40), not max",
          abs(newm3[P][vid["S7"]] - 0.40) < 5e-5 and abs(old3[P][vid["S7"]] - 0.60) < 5e-5)
    check("fallback: canonical rows of two UniProt entries -> missing, not max",
          abs(newm3[P][vid["S6"]] - gm) < 5e-5)
    check("fallback: --mane excludes an isoform MANE puts in another gene",
          abs(newm3[P][vid["S5"]] - gm) < 5e-5 and abs(new3[P][vid["S5"]] - 0.80) < 5e-5)
    check("fallback: GENEQ's matched scores identical under both rules",
          all(old3[Q][v] == new3[Q][v] == newm3[Q][v] for v in q_ids))
    check("fallback: matched rows (S3, M1-M3) bit-identical to the legacy rule",
          all(old3[P][vid[k]] == new3[P][vid[k]] == newm3[P][vid[k]]
              for k in ("S3", "M1", "M2", "M3")))
    check("gtf: WITHOUT it, GENEZ takes GENEQ's unannotated canonical (0.90)",
          abs(newm3["GENEZ"][vid["S8"]] - 0.90) < 5e-5)
    check("gtf: WITH it, that row is excluded; GENEZ has no score -> 0.0 fill",
          newg3["GENEZ"][vid["S8"]] == 0.0)
    check("gtf: the log names GENEZ as left with no scored missense",
          "GENEZ (1 of 1 missense)" in open(gtf_log).read())
    check("gtf: GENEP/GENEQ unchanged (own rows kept; retired TZ1/TZ2 still averaged)",
          all(newg3[g][v] == newm3[g][v] for g in (P, Q) for v in newm3[g]))
    check("provenance: header carries the annotation context and the AM source",
          prows[0] == ["GENE", "ID", "ANNOTATION", "CSQ", "ANNO_TRANSCRIPT",
                       "SOURCE", "TRANSCRIPT", "AM"])
    check("provenance: one row per missense variant, bin/CSQ/transcript from the table",
          len(prows) - 1 == len(p_ids) + len(q_ids) + 1 and
          pz == [["GENEZ", vid["S8"], "missense", "missense_variant",
                  "ENST000000TZZ", "missing", ".", "NA"]])
    # ---- NON-MISSENSE CONTROL (--nonmissense drop)
    # One region, all labelled "missense": n1/n2 missense SNVs matched at 0.20/0.60, n6 a
    # missense SNV with no AlphaMissense row (gene mean), n3 an in-frame deletion, n5 an LC
    # frameshift, n4 an intronic SNV whose only row (0.99) is another transcript's canonical.
    nm = {"n1": "20:7001:A:G", "n2": "20:7002:A:G", "n6": "20:7006:A:G",
          "n3": "20:7003:AGT:A", "n5": "20:7005:AG:A", "n4": "20:7004:A:G"}
    nm_csq = {"n1": "missense_variant", "n2": "missense_variant", "n6": "missense_variant",
              "n3": "inframe_deletion", "n5": "frameshift_variant", "n4": "intron_variant"}
    order_nm = ["n1", "n2", "n3", "n4", "n5", "n6"]
    gf4 = os.path.join(tmp, "group_nonmis.txt")
    with open(gf4, "w") as fh:
        fh.write("GENEM var %s 20:7010:A:T\n" % " ".join(nm[k] for k in order_nm))
        fh.write("GENEM anno %s pLoF\n" % " ".join(["missense"] * len(order_nm)))
        # n7: a missense SNV whose CSQ lists every VEP term (missense_variant&splice_region_variant), in its
        # own gene so GENEM's numbers stay as they were; it is missense and keeps its AlphaMissense 0.80
        fh.write("GENEN var 20:7007:A:G 20:7011:A:T\nGENEN anno missense pLoF\n")
    anno4 = os.path.join(tmp, "anno_nonmis.tsv")
    with open(anno4, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ",
                            "TRANSCRIPT", "MANE_SELECT", "CANONICAL", "BIOTYPE", "DS_MAX",
                            "ANNOTATION"]) + "\n")
        for k in order_nm:
            fh.write("\t".join([nm[k], "GENEM", "LC" if k == "n5" else "NA", "0", "0",
                                nm_csq[k], "ENST00000000TM", "NM_1", "YES", "protein_coding",
                                "0", "missense"]) + "\n")
        fh.write("\t".join(["20:7007:A:G", "GENEN", "NA", "0", "0", "missense_variant&splice_region_variant",
                            "ENST00000000TN", "NM_2", "YES", "protein_coding", "0", "missense"]) + "\n")
    am4c = os.path.join(tmp, "am_nonmis_canon.tsv.gz")
    with gzip.open(am4c, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        fh.write("\t".join(["chr20", "7004", "A", "G", "hg38", "P7", "ENST00000000TX.1",
                            "A1B", "0.9900", "pathogenic"]) + "\n")
    am4i = os.path.join(tmp, "am_nonmis_iso.tsv.gz")
    with gzip.open(am4i, "wt") as fh:
        fh.write("\t".join(["#CHROM", "POS", "REF", "ALT", "genome", "transcript_id",
                            "protein_variant", "am_pathogenicity", "am_class"]) + "\n")
        for k, sc in (("n1", 0.20), ("n2", 0.60)):
            fh.write("\t".join(["chr20", nm[k].split(":")[1], "A", "G", "hg38",
                                "ENST00000000TM.1", "A1B", "%.4f" % sc, "x"]) + "\n")
        fh.write("\t".join(["chr20", "7007", "A", "G", "hg38", "ENST00000000TN.1", "A1B", "0.8000", "x"]) + "\n")

    def run4(mode):
        o = os.path.join(tmp, "nonmis_%s.txt" % mode)
        c = build(gf4, o, am4c, ["pLoF"], ["missense"], score_name="AM",
                  anno_table_path=anno4, isoforms_path=am4i, prefer_anno_transcript=True,
                  nonmissense=mode, log=open(os.devnull, "w"))
        g = per_gene(o)
        return g["GENEM"], c, list(read_group_regions(o)), g.get("GENEN", {})

    imp, c_imp, _, _ = run4("impute")
    drp, c_drp, rg, drp_n = run4("drop")
    check("nonmissense: impute (old) gives the intronic SNV another transcript's 0.99",
          abs(imp[nm["n4"]] - 0.99) < 5e-5)
    check("nonmissense: impute fills indel and LC frameshift with the gene mean",
          abs(imp[nm["n3"]] - (0.2 + 0.6 + 0.99) / 3) < 5e-5 and imp[nm["n5"]] == imp[nm["n3"]])
    check("nonmissense: drop removes the 3 non-missense variants, keeps the rest",
          c_drp["nonmissense_dropped"] == 3 and
          set(drp) == {nm["n1"], nm["n2"], nm["n6"], "20:7010:A:T"})
    check("nonmissense: drop leaves matched missense scores bit-identical",
          drp[nm["n1"]] == imp[nm["n1"]] and drp[nm["n2"]] == imp[nm["n2"]])
    check("nonmissense: the gene mean is taken over missense only (0.40, not 0.597)",
          abs(drp[nm["n6"]] - 0.40) < 5e-5 and abs(imp[nm["n6"]] - (0.2 + 0.6 + 0.99) / 3) < 5e-5)
    check("nonmissense: drop keeps a missense SNV whose CSQ lists several terms "
          "(missense_variant&splice_region_variant) with its own AlphaMissense 0.80",
          abs(drp_n.get("20:7007:A:G", -1) - 0.80) < 5e-5)
    check("nonmissense: var, anno and score lines stay aligned after the drop",
          all(len({len(v) for g2, k, n, v, raw in rg if g2 == g and k in ("var", "anno", "score")}) == 1
              for g in {r[0] for r in rg}) and
          all(len(v) == 4 for g, k, n, v, raw in rg if g == "GENEM" and k in ("var", "anno", "score")))
    check("fallback: no non-imputed score moves UP against the legacy max",
          all(new3[P][vid[k]] <= old3[P][vid[k]] + 5e-5 for k in ("S1", "S3", "S4", "S5")))

    # ---- PERCENTILE SCALE (--am-scale percentile, --nonmissense gene-mean; 2026-10-21)
    # G1: annotation (MANE) transcript T1 is in the ISOFORM file in two versions -- .1 (a decoy whose
    # p2 is 0.95) and .2 (the scale: p1..p5 and p9 = 0.1 0.3 0.3 0.5 0.9 0.7) -- and AlphaMissense's
    # canonical file has another G1 transcript, TC1, scoring p6 (0.4; no T1 row at p6).
    # Group file: p2, p5 missense; p6 missense via the canonical fallback; p7 in-frame deletion and p8
    # LOFTEE-LC frameshift (no row: gene mean); p9 start_lost SNV WITH a T1.2 row; a pLoF; a synonymous.
    # By hand on T1.2's sorted scale [0.1 0.3 0.3 0.5 0.7 0.9], N = 6: p2 = (2+3)/2/6, p5 = 6/6,
    # p9 = 5/6, p6 = (3 + 0.5)/6 (between 0.3 and 0.5); gene mean of those four = p7 = p8.
    # G2: annotation transcript T2 is not in AlphaMissense; its canonical-file transcript TC2 (0.1 0.2
    # 0.2 0.9) is the scale: q1 (0.9) = 4/4, q2 (0.2) = (1 + 1.5)/4. G3: no AlphaMissense at all.
    pos = {"p%d" % i: "20:%d:A:G" % (8000 + i) for i in range(1, 10)}
    pos.update({"q1": "20:8101:A:G", "q2": "20:8102:A:G", "r1": "20:8201:A:G"})
    pos["p7"] = "20:8007:AGT:A"; pos["p8"] = "20:8008:AG:A"
    pos["s1"] = "20:8011:C:T"   # a SpliceAI-route intronic variant in the missense bin: dropped
    gfp = os.path.join(tmp, "group_pct.txt")
    g1 = ["p2", "p5", "p6", "p7", "p8", "p9", "s1"]
    with open(gfp, "w") as fh:
        fh.write("G1 var %s 20:8050:C:T 20:8060:C:T\nG1 anno %s pLoF synonymous\n" % (
            " ".join(pos[k] for k in g1), " ".join(["missense"] * len(g1))))
        fh.write("G2 var %s %s\nG2 anno missense missense\n" % (pos["q1"], pos["q2"]))
        fh.write("G3 var %s\nG3 anno missense\n" % pos["r1"])
        fh.write("G4 var 20:8301:A:G\nG4 anno missense\n")
        fh.write("G5 var 20:8401:A:G\nG5 anno missense\n")
        fh.write("G8 var 20:8501:A:G\nG8 anno missense\n")
    csq_p = {"p7": "inframe_deletion", "p8": "frameshift_variant", "p9": "start_lost", "s1": "intron_variant"}
    annop = os.path.join(tmp, "anno_pct.tsv")
    with open(annop, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ", "TRANSCRIPT",
                            "MANE_SELECT", "CANONICAL", "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        for k in g1:
            fh.write("\t".join([pos[k], "G1", "LC" if k == "p8" else "NA", "0", "0", csq_p.get(k, "missense_variant"),
                                "T1", "NM_1", "YES", "protein_coding", "0", "missense"]) + "\n")
        for k in ("q1", "q2"):
            fh.write("\t".join([pos[k], "G2", "NA", "0", "0", "missense_variant", "T2", "NA", "YES",
                                "protein_coding", "0", "missense"]) + "\n")
        fh.write("\t".join([pos["r1"], "G3", "NA", "0", "0", "missense_variant", "T3", "NA", "YES",
                            "protein_coding", "0", "missense"]) + "\n")
        fh.write("\t".join(["20:8301:A:G", "G4", "NA", "0", "0", "missense_variant", "T4", "NA", "YES",
                            "protein_coding", "0", "missense"]) + "\n")
        fh.write("\t".join(["20:8401:A:G", "G5", "NA", "0", "0", "missense_variant", "T5", "NA", "YES",
                            "protein_coding", "0", "missense"]) + "\n")
        fh.write("\t".join(["20:8501:A:G", "G8", "NA", "0", "0", "missense_variant", "T8", "NM_8", "YES",
                            "protein_coding", "0", "missense"]) + "\n")
    amc_p = os.path.join(tmp, "am_pct_canon.tsv.gz")
    with gzip.open(amc_p, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        fh.write("\t".join(["chr20", "8006", "A", "G", "hg38", "P1", "TC1.3", "A1B", "0.4000", "x"]) + "\n")
        for p_, sc in (("8101", 0.9), ("8102", 0.2), ("8103", 0.2), ("8104", 0.1)):
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "P2", "TC2.1", "A1B", "%.4f" % sc, "x"]) + "\n")
        for p_, sc in (("8401", 0.5), ("8402", 0.1), ("8403", 0.9)):
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "P5", "TR5.1", "A1B", "%.4f" % sc, "x"]) + "\n")
        for p_, sc in (("8410", 0.2), ("8411", 0.3), ("8412", 0.6), ("8413", 0.7), ("8414", 0.8)):   # larger, no group variant
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "P5", "TR5B.1", "A1B", "%.4f" % sc, "x"]) + "\n")
    ami_p = os.path.join(tmp, "am_pct_iso.tsv.gz")
    with gzip.open(ami_p, "wt") as fh:
        fh.write("\t".join(["#CHROM", "POS", "REF", "ALT", "genome", "transcript_id",
                            "protein_variant", "am_pathogenicity", "am_class"]) + "\n")
        for k, sc in (("p1", 0.1), ("p2", 0.3), ("p3", 0.3), ("p4", 0.5), ("p5", 0.9), ("p9", 0.7)):
            fh.write("\t".join(["chr20", pos[k].split(":")[1], "A", "G", "hg38", "T1.2", "A1B", "%.4f" % sc, "x"]) + "\n")
        fh.write("\t".join(["chr20", "8002", "A", "G", "hg38", "T1.1", "A1B", "0.9500", "x"]) + "\n")
        for p_, sc in (("8301", 0.6), ("8302", 0.2), ("8303", 0.4), ("8304", 0.8)):
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "TI4.1", "A1B", "%.4f" % sc, "x"]) + "\n")
        # G8: T8.1 (the GTF's version: 8501 = 0.30 on 0.1 0.3 0.5 0.7) and T8.3 (higher, listed FIRST: 8501 = 0.90 on
        # 0.2 0.9): the annotation's version must win over both "first in the file" and "highest"
        for p_, sc in (("8501", 0.9), ("8502", 0.2)):
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "T8.3", "A1B", "%.4f" % sc, "x"]) + "\n")
        for p_, sc in (("8501", 0.3), ("8502", 0.1), ("8503", 0.5), ("8504", 0.7)):
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "T8.1", "A1B", "%.4f" % sc, "x"]) + "\n")
    gtf_p = os.path.join(tmp, "pct.gtf")
    with open(gtf_p, "w") as fh:
        for t, g in (("T1", "G1"), ("TC1", "G1"), ("T2", "G2"), ("TC2", "G2"), ("T3", "G3"), ("T4", "G4"), ("TI4", "G4")):
            fh.write('chr20\tHAVANA\ttranscript\t1\t2\t.\t+\t.\tgene_id "%s.1"; transcript_id "%s.%s";\n'
                     % (g, t, "2" if t == "T1" else "1"))
        fh.write('chr20\tHAVANA\ttranscript\t1\t2\t.\t+\t.\tgene_id "G8.1"; transcript_id "T8.1";\n')
        fh.write('chr20\tHAVANA\ttranscript\t8400\t8500\t.\t+\t.\tgene_id "G5.1"; transcript_id "T5.1";\n')
    prov_p = os.path.join(tmp, "prov_pct.tsv")

    def run_pct(**kw):
        o = os.path.join(tmp, "pct_%d.txt" % len(os.listdir(tmp)))
        c = build(gfp, o, amc_p, ["pLoF"], ["missense"], score_name="AM", anno_table_path=annop,
                  isoforms_path=ami_p, prefer_anno_transcript=True, gtf_path=gtf_p, log=open(os.devnull, "w"), **kw)
        return per_gene(o), c

    pg, cp = run_pct(am_scale="percentile", nonmissense="gene-mean", provenance=prov_p)
    want = {"p2": 2.5 / 6, "p5": 1.0, "p9": 5 / 6, "p6": 3.5 / 6}
    gm = sum(want.values()) / 4
    check("percentile: matched missense ranked on the gene's annotation transcript, ties at mean rank",
          abs(pg["G1"][pos["p2"]] - 2.5 / 6) < 5e-5 and pg["G1"][pos["p5"]] == 1.0)
    check("percentile: the HIGHER transcript version gives score and scale (decoy .1 ignored)",
          abs(pg["G1"][pos["p2"]] - 2.5 / 6) < 5e-5)
    check("percentile: a canonical-fallback score is placed on the annotation scale",
          abs(pg["G1"][pos["p6"]] - 3.5 / 6) < 5e-5)
    check("percentile: start_lost SNV with a row on its own transcript keeps its percentile",
          abs(pg["G1"][pos["p9"]] - 5 / 6) < 5e-5)
    check("percentile: in-frame del and LC frameshift take the gene mean ON THE PERCENTILE SCALE",
          abs(pg["G1"][pos["p7"]] - gm) < 5e-5 and pg["G1"][pos["p8"]] == pg["G1"][pos["p7"]])
    check("percentile: a SpliceAI-route intronic variant (no protein change) is dropped, not imputed",
          pos["s1"] not in pg["G1"])
    check("percentile: pLoF 1.0, synonymous 0.0 placeholder",
          pg["G1"]["20:8050:C:T"] == 1.0 and pg["G1"]["20:8060:C:T"] == 0.0)
    check("percentile: annotation transcript unscored -> canonical-file transcript's scale",
          pg["G2"][pos["q1"]] == 1.0 and abs(pg["G2"][pos["q2"]] - 2.5 / 4) < 5e-5)
    check("percentile: only another isoform scored -> ranked on that isoform (0.6 = 3/4), not filled",
          abs(pg["G4"]["20:8301:A:G"] - 0.75) < 5e-5)
    check("percentile: scores only from transcripts GENCODE lacks -> the LARGEST such in the gene's span, "
          "even if this cohort's variants never touch it (0.5 on TR5B 0.2 0.3 0.6 0.7 0.8 = 2.5/5; on TR5 2/3)",
          abs(pg["G5"]["20:8401:A:G"] - 0.5) < 5e-5)
    check("percentile: two versions of the annotation transcript -> the GTF's (lower, listed second) for score "
          "AND scale (0.30 on 0.1 0.3 0.5 0.7 = 2/4; first-in-file or highest would give 1.0)",
          abs(pg["G8"]["20:8501:A:G"] - 0.5) < 5e-5)
    check("percentile: a gene with no AlphaMissense transcript is filled at 0.5 (the expected within-gene percentile)",
          pg["G3"][pos["r1"]] == 0.5)
    check("percentile: --empty-gene-fill overrides that (0.25)",
          run_pct(am_scale="percentile", nonmissense="gene-mean", empty_gene_fill=0.25)[0]["G3"][pos["r1"]] == 0.25)
    # paralog rows (alphamissense_paralogs.py): G3's own transcript T3 (no AlphaMissense row) gets rows from a further
    # file in the isoforms layout -- 8201 = 0.9 on 0.2 0.6 0.9 -> (2 + 1) / 3 = 1.0 (not the 0.5 fill)
    exi = os.path.join(tmp, "am_paralogs.tsv.gz")
    with gzip.open(exi, "wt") as fh:
        fh.write("#CHROM\tPOS\tREF\tALT\tgenome\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n")
        for p_, sc in (("8201", 0.9), ("8202", 0.2), ("8203", 0.6)):
            fh.write("\t".join(["chr20", p_, "A", "G", "hg38", "T3.1", "A1B", "%.4f" % sc, "x"]) + "\n")
    pe, _ = run_pct(am_scale="percentile", nonmissense="gene-mean", extra_isoforms=[exi])
    check("paralogs: --extra-isoforms rows score a gene AlphaMissense lacks, on its own transcript's scale (1.0, not 0.5)",
          abs(pe["G3"][pos["r1"]] - 1.0) < 5e-5 and pe["G1"] == pg["G1"] and pe["G2"] == pg["G2"])
    check("raw: a gene with no AlphaMissense row stays 0.0 (released raw files reproduced)",
          run_pct(am_scale="raw", nonmissense="gene-mean")[0]["G3"][pos["r1"]] == 0.0)

    # pseudoautosomal genes: AlphaMissense files them under chrY ('_PAR_Y'), with no chrX row. GX (chrX, TX): variants
    # at PAR1 X:1001 / 1002 (same coordinates on Y) and PAR2 X:155701483 (= Y:56888003); the chrY file also holds a
    # Y-specific transcript (TY) at Y:1001 and a PAR row of another transcript -- neither may reach GX's scale.
    gx = os.path.join(tmp, "group_par.txt"); ax = os.path.join(tmp, "anno_par.tsv")
    xids = ["X:1001:A:G", "X:1002:A:G", "X:155701483:A:G"]
    with open(gx, "w") as fh:
        fh.write("GX var %s X:2000:A:T\nGX anno missense missense missense pLoF\n" % " ".join(xids))
    with open(ax, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ", "TRANSCRIPT",
                            "MANE_SELECT", "CANONICAL", "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        for v in xids:
            fh.write("\t".join([v, "GX", "NA", "0", "0", "missense_variant", "TX", "NM_X", "YES", "protein_coding", "0",
                                "missense"]) + "\n")
    amx, amy = os.path.join(tmp, "am_par_x.tsv.gz"), os.path.join(tmp, "am_par_y.tsv.gz")
    with gzip.open(amx, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        fh.write("\t".join(["chrX", "5000000", "A", "G", "hg38", "PO", "TO.1", "A1B", "0.5000", "x"]) + "\n")
    with gzip.open(amy, "wt") as fh:
        fh.write("\t".join(AM_CANONICAL_COLS) + "\n")
        for p_, sc in (("1001", 0.2), ("1002", 0.8), ("1003", 0.4), ("56888003", 0.6)):
            fh.write("\t".join(["chrY", p_, "A", "G", "hg38", "PX", "TX.1_PAR_Y", "A1B", "%.4f" % sc, "x"]) + "\n")
        fh.write("\t".join(["chrY", "1001", "A", "G", "hg38", "PY", "TY.1", "A1B", "0.9900", "x"]) + "\n")
        fh.write("\t".join(["chrY", "20000000", "A", "G", "hg38", "PX", "TX.1_PAR_Y", "A1B", "0.9900", "x"]) + "\n")
    isx = os.path.join(tmp, "am_par_iso.tsv.gz")
    with gzip.open(isx, "wt") as fh:
        fh.write("\t".join(["#CHROM", "POS", "REF", "ALT", "genome", "transcript_id",
                            "protein_variant", "am_pathogenicity", "am_class"]) + "\n")
    gtfx = os.path.join(tmp, "par.gtf")
    with open(gtfx, "w") as fh:
        fh.write('chrX\tHAVANA\ttranscript\t1000\t155701500\t.\t+\t.\tgene_id "GX.1"; transcript_id "TX.1";\n')

    def run_par(**kw):
        o = os.path.join(tmp, "par_%d.txt" % len(os.listdir(tmp)))
        build(gx, o, amx, ["pLoF"], ["missense"], score_name="AM", anno_table_path=ax, prefer_anno_transcript=True,
              isoforms_path=isx,
              gtf_path=gtfx, chrom_filter="X", am_scale="percentile", nonmissense="gene-mean",
              log=open(os.devnull, "w"), **kw)
        return per_gene(o)["GX"]
    px = run_par(par_y_am=amy)
    # scale 0.2 0.4 0.6 0.8 (PAR1 three, PAR2 one; the Y-specific TY row and the out-of-PAR row excluded)
    check("PAR: chrX variants scored from AlphaMissense's chrY '_PAR_Y' rows (PAR1 same coordinates, PAR2 shifted)",
          abs(px["X:1001:A:G"] - 0.25) < 5e-5 and abs(px["X:1002:A:G"] - 1.0) < 5e-5
          and abs(px["X:155701483:A:G"] - 0.75) < 5e-5)
    check("PAR: without --par-y-am the same gene has no scale (0.5 fill) -- the gap this closes",
          {v for k, v in run_par().items() if k in xids} == {0.5})
    check("percentile: counts -- 2 non-missense filled, 1 matched (start_lost)",
          cp["nonmissense_filled"] == 2 and cp["nonmissense_matched"] == 1)
    pr = [l.rstrip("\n").split("\t") for l in open(prov_p)]
    ixp = {k: i for i, k in enumerate(pr[0])}
    check("percentile: provenance REF_TRANSCRIPT/SCORE; SCORE equals the file's value",
          pr[0][-2:] == ["REF_TRANSCRIPT", "SCORE"] and
          all(abs(float(f[ixp["SCORE"]]) - pg[f[0]][f[1]]) < 5e-5 for f in pr[1:] if f[ixp["SCORE"]] != "NA") and
          {f[ixp["REF_TRANSCRIPT"]] for f in pr[1:] if f[0] == "G1"} == {"T1.2"} and
          {f[ixp["REF_TRANSCRIPT"]] for f in pr[1:] if f[0] == "G2"} == {"TC2.1"})
    scale = [0.1, 0.3, 0.3, 0.5, 0.7, 0.9]
    by_rank = apply_transform(list(enumerate(scale)), "percentile-in-gene")
    check("percentile: percentile_in == average-rank/N transform on the scale (two routes)",
          all(abs(percentile_in(x, sorted(scale)) - by_rank[i]) < 1e-12 for i, x in enumerate(scale)))
    rg_raw, c_raw = run_pct(am_scale="raw", nonmissense="gene-mean")
    check("percentile: raw mode on the same input writes am_pathogenicity (p5 0.90)",
          abs(rg_raw["G1"][pos["p5"]] - 0.9) < 5e-5)
    try:
        build(gfp, os.path.join(tmp, "nogtf.txt"), amc_p, ["pLoF"], ["missense"], anno_table_path=annop,
              isoforms_path=ami_p, prefer_anno_transcript=True, am_scale="percentile", log=open(os.devnull, "w"))
        refused = False
    except SystemExit:
        refused = True
    check("percentile: refused without --gtf (fallback scale cannot be mapped)", refused)
    # column mode with --nonmissense gene-mean (a second weight set, e.g. scallion): the same variants as the
    # AlphaMissense file -- missense from the column, protein-altering non-missense at the gene mean of the column
    # values, the intronic s1 dropped -- and the CLI path (which once read a count column mode never set)
    annoc = os.path.join(tmp, "anno_col.tsv")
    colv = {"p2": 0.2, "p5": 0.8, "p6": 0.5, "p9": 0.6}
    with open(annop) as fi, open(annoc, "w") as fo:
        fo.write(fi.readline().rstrip("\n") + "\tXSCORE\n")
        for l in fi:
            f = l.rstrip("\n").split("\t"); k = next((k for k, v in pos.items() if v == f[0]), None)
            fo.write("\t".join(f + ["%.4f" % colv[k] if k in colv else "NA"]) + "\n")
    oc = os.path.join(tmp, "col_gm.txt")
    cc = build(gfp, oc, None, ["pLoF"], ["missense"], score_name="X", anno_table_path=annoc, score_column="XSCORE",
               score_transform="none", zero_is_missing=False, nonmissense="gene-mean", log=open(os.devnull, "w"))
    gc = per_gene(oc)
    check("column gene-mean: missense from the column, in-frame del / LC frameshift at their mean, intronic dropped",
          abs(gc["G1"][pos["p2"]] - 0.2) < 5e-5 and abs(gc["G1"][pos["p7"]] - 0.525) < 5e-5 and
          gc["G1"][pos["p8"]] == gc["G1"][pos["p7"]] and pos["s1"] not in gc["G1"] and cc["nonmissense_dropped"] == 1)
    import subprocess
    rc = subprocess.run([sys.executable, os.path.abspath(__file__), "--group", gfp, "--out", os.path.join(tmp, "col_cli.txt"),
                         "--anno-table", annoc, "--score-column", "XSCORE", "--score-transform", "none", "--zero-is-value",
                         "--lof-anno", "pLoF", "--missense-anno", "missense", "--name", "X", "--nonmissense", "gene-mean"],
                        capture_output=True, text=True)
    check("column gene-mean: the command line runs to the end (exit 0) and writes the same file",
          rc.returncode == 0 and open(os.path.join(tmp, "col_cli.txt")).read() == open(oc).read())
    big = os.path.join(tmp, "big_gene.txt")   # one gene with 40,000 scale variants: its bottom one is 1/40,000
    with open(big, "w") as fh:
        fh.write("G6 var 20:9000:A:G\nG6 anno missense\n")
    annob = os.path.join(tmp, "anno_big.tsv")
    with open(annob, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED", "CSQ", "TRANSCRIPT", "MANE_SELECT",
                            "CANONICAL", "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        fh.write("\t".join(["20:9000:A:G", "G6", "NA", "0", "0", "missense_variant", "T6", "NM_6", "YES",
                            "protein_coding", "0", "missense"]) + "\n")
    amib = os.path.join(tmp, "am_big_iso.tsv.gz")
    with gzip.open(amib, "wt") as fh:
        fh.write("\t".join(["#CHROM", "POS", "REF", "ALT", "genome", "transcript_id", "protein_variant",
                            "am_pathogenicity", "am_class"]) + "\n")
        fh.write("chr20\t9000\tA\tG\thg38\tT6.1\tA1B\t0.0001\tx\n")
        for i in range(39999):
            fh.write("chr20\t%d\tA\tG\thg38\tT6.1\tA1B\t%.4f\tx\n" % (10000 + i, 0.5 + i / 1e6))
    gtfb = os.path.join(tmp, "big.gtf")
    open(gtfb, "w").write('chr20\tH\ttranscript\t1\t2\t.\t+\t.\tgene_id "G6.1"; transcript_id "T6.1";\n')
    ob = os.path.join(tmp, "big_out.txt")
    build(big, ob, amc_p, [], ["missense"], anno_table_path=annob, isoforms_path=amib,
          prefer_anno_transcript=True, gtf_path=gtfb, am_scale="percentile", log=open(os.devnull, "w"))
    check("percentile: a 40,000-variant gene's bottom variant is written as 1/40000 (6 decimals), not 0.0000",
          per_gene(ob)["G6"]["20:9000:A:G"] == 0.000025)
    check("percentile: a fallback score above the whole scale is 1.0, below it 0.5/N (never outside [0, 1])",
          percentile_in(0.99, sorted(scale)) == 1.0 and abs(percentile_in(0.01, sorted(scale)) - 0.5 / 6) < 1e-12)

    # ---- ID FORMATS of the BRaVa biobank group files
    for vid, fmt, want in (("chr20:881:A:T", "colon", ("20", "881", "A", "T")),
                           ("chrX:100629994_TAAAG_T", "colon_underscore", ("X", "100629994", "TAAAG", "T")),
                           ("chr21_41801194_C_G", "underscore", ("21", "41801194", "C", "G")),
                           ("1_169853716_C_A", "underscore", ("1", "169853716", "C", "A"))):
        check("id format: %-24s detected %s and parsed" % (vid, fmt),
              detect_id_format(vid) == fmt and parse_variant_id(vid, fmt) == want)
    check("id format: an alt-contig name does not parse as a coordinate",
          parse_variant_id("chr1_KI270706v1_random_100_A", "underscore") is None)

    print("\n%s" % ("selftest PASSED" if ok else "selftest FAILED"))
    return 0 if ok else 1


# --------------------------------------------------------------------- cli

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", help="input SAIGE group file (.txt or .txt.gz)")
    ap.add_argument("--out", help="output group file (plain text)")
    ap.add_argument("--am", help="AlphaMissense_hg38.tsv.gz (canonical transcripts)")
    ap.add_argument("--isoforms", help="AlphaMissense_isoforms_hg38.tsv.gz; with "
                                       "--mane, prefer the MANE Select transcript")
    ap.add_argument("--mane", help="MANE summary file (MANE.GRCh38.vX.summary.txt.gz). "
                                   "With --prefer-anno-transcript, a fallback row "
                                   "whose transcript MANE puts in another gene is "
                                   "never used.")
    ap.add_argument("--lof-anno", default="pLoF",
                    help="comma-separated annotation labels scored 1.0 by "
                         "construction (default: pLoF). Pass the same labels to "
                         "--flexRV_lofAnno.")
    ap.add_argument("--missense-anno", default="damaging_missense,other_missense",
                    help="comma-separated annotation labels joined to "
                         "AlphaMissense (default: damaging_missense,other_missense)")
    ap.add_argument("--name", default="AM",
                    help="score set name; the run needs --flexRV_score=<name> "
                         "(default: AM)")
    ap.add_argument("--missing", default="gene-mean",
                    help="missense with no AlphaMissense row: 'gene-mean' "
                         "(default), 'drop', or a literal value in [0, 1]")
    ap.add_argument("--empty-gene-fill", type=float, default=None,
                    help="with --missing gene-mean, the value for missense variants of a gene with NO scored missense "
                         "variant (no mean to take): default 0.5 with --am-scale percentile (the expected within-gene "
                         "percentile), else 0.0")
    ap.add_argument("--par-y-am", help="with --chrom X: an AlphaMissense canonical file holding chrY; its '_PAR_Y' "
                    "rows (pseudoautosomal genes, which AlphaMissense files under chrY) are read as chrX rows")
    ap.add_argument("--par-y-isoforms", help="the same for the isoforms file")
    ap.add_argument("--extra-isoforms", action="append", default=[],
                    help="a further file in the AlphaMissense isoforms layout read as isoform rows (repeatable): "
                         "alphamissense_paralogs.py output, AlphaMissense's protein-level scores placed on paralog "
                         "copies it does not score")
    ap.add_argument("--other", type=float, default=0.0,
                    help="score for every other annotation (synonymous, "
                         "non_coding, ...); these are dropped by "
                         "--annotation_in_groupTest before any decode, so this "
                         "is a placeholder (default: 0.0)")
    ap.add_argument("--id-format", choices=ID_FORMATS,
                    help="force the variant-ID coordinate format (default: "
                         "auto-detect from the first ID)")
    ap.add_argument("--anno-table", help="BRaVa long-form annotation table "
                    "(ID/GENE/TRANSCRIPT/MANE_SELECT/CSQ ...). Supplies the "
                    "transcript each annotation was called against, so the "
                    "transcript AlphaMissense scored can be CHECKED against it "
                    "and the coverage gap explained by consequence class.")
    ap.add_argument("--score-column",
                    help="read the score from this column of --anno-table "
                         "instead of joining AlphaMissense (e.g. REVEL_SCORE, "
                         "CADD_PHRED, DS_MAX). A 0 in these columns means "
                         "ABSENT, not 'benign', and is treated as missing.")
    ap.add_argument("--score-transform",
                    help="none | percentile-in-gene | cadd-phred | minmax:LO,HI. "
                         "Default depends on the column: CADD_PHRED is unbounded "
                         "so it defaults to percentile-in-gene; REVEL_SCORE and "
                         "DS_MAX are already [0,1] and default to none.")
    ap.add_argument("--zero-is-missing", dest="zim", action="store_true",
                    default=None, help="force a 0 in --score-column to count as "
                                       "absent (per-column default otherwise)")
    ap.add_argument("--zero-is-value", dest="zim", action="store_false",
                    help="force a 0 in --score-column to count as a real score")
    ap.add_argument("--drop-anno", default="",
                    help="comma-separated annotation labels to REMOVE from the "
                         "output entirely (e.g. non_coding). These are markers "
                         "--annotation_in_groupTest would discard in R before "
                         "any decode, so dropping them changes no result while "
                         "shrinking the file; a region left with nothing is "
                         "omitted and named.")
    ap.add_argument("--prefer-anno-transcript", action="store_true",
                    help="choose the AlphaMissense row scored on the transcript "
                         "the annotation table used for that (gene, variant) -- "
                         "i.e. MANE Select where the annotation has one, its "
                         "canonical transcript otherwise. Requires --anno-table; "
                         "pair with --isoforms so the MANE transcript can be "
                         "found when it is not AlphaMissense's canonical one.")
    ap.add_argument("--chrom", help="only read AlphaMissense rows on this "
                                    "chromosome (large speedup on the 71M-row file)")
    ap.add_argument("--gtf", help="GENCODE GTF of the release the --anno-table was "
                    "built on (v39 for AoU v8). With --prefer-anno-transcript, a "
                    "fallback row whose transcript it puts in another gene is never "
                    "used -- e.g. a readthrough gene taking its partner's score.")
    ap.add_argument("--nonmissense", choices=("impute", "drop", "gene-mean"), default="impute",
                    help="missense-bin variants whose --anno-table CSQ is not "
                         "missense_variant (in-frame indels, LOFTEE-LC LoF, start/stop "
                         "lost, SpliceAI-route variants): 'impute' (default; AlphaMissense "
                         "from any transcript, else --missing), 'drop' (remove them: "
                         "AlphaMissense scores missense substitutions only) or 'gene-mean' "
                         "(AlphaMissense only from a row on their own transcript, else the "
                         "gene mean of the scored missense entries, on the --am-scale scale)")
    ap.add_argument("--am-scale", choices=("raw", "percentile"), default="raw",
                    help="'raw' (default): am_pathogenicity. 'percentile': the FlexRV paper's "
                         "score, the percentile among all possible missense variants of the "
                         "gene -- every AlphaMissense row on the gene's annotation transcript "
                         "(MANE Select, else canonical), the canonical-file transcript where "
                         "AlphaMissense does not score it; needs --anno-table, "
                         "--prefer-anno-transcript, --isoforms and --gtf")
    ap.add_argument("--provenance", help="write one line per missense-bin variant: "
                    "GENE, ID, ANNOTATION, CSQ, ANNO_TRANSCRIPT (from --anno-table) "
                    "and SOURCE, TRANSCRIPT, AM (which AlphaMissense row its score "
                    "came from, before imputation)")
    # Reproduces group files built before 2026-10 (fallback = max over every
    # row at the coordinate); biased upward, so not offered in --help.
    ap.add_argument("--legacy-max-fallback", action="store_true",
                    help=argparse.SUPPRESS)
    ap.add_argument("--replace-existing", action="store_true",
                    help="overwrite an existing score line instead of refusing")
    ap.add_argument("--selftest", action="store_true",
                    help="run the positive/negative controls and exit")
    a = ap.parse_args()

    if a.selftest:
        return _selftest()
    need = ("group", "out") if a.score_column else ("group", "out", "am")
    for req in need:
        if getattr(a, req) is None:
            ap.error("--%s is required (or pass --selftest)" % req)

    counts = build(a.group, a.out, a.am,
                   [s for s in a.lof_anno.split(",") if s],
                   [s for s in a.missense_anno.split(",") if s],
                   score_name=a.name, missing=a.missing, other=a.other,
                   id_format=a.id_format, chrom_filter=a.chrom,
                   par_y_am=a.par_y_am, par_y_isoforms=a.par_y_isoforms, empty_gene_fill=a.empty_gene_fill,
                   extra_isoforms=a.extra_isoforms,
                   isoforms_path=a.isoforms, mane_path=a.mane,
                   replace_existing=a.replace_existing,
                   anno_table_path=a.anno_table,
                   prefer_anno_transcript=a.prefer_anno_transcript,
                   drop_anno=[x for x in a.drop_anno.split(",") if x],
                   score_column=a.score_column,
                   score_transform=a.score_transform, zero_is_missing=a.zim,
                   fallback="max" if a.legacy_max_fallback else "safe",
                   provenance=a.provenance, gtf_path=a.gtf, nonmissense=a.nonmissense,
                   am_scale=a.am_scale)

    print("# wrote %s" % a.out)
    print("# %(genes)d regions, %(variants)d variants: %(lof)d lof -> 1.0, "
          "%(missense)d missense (%(missense_scored)d scored, "
          "%(missense_imputed)d imputed, %(missense_dropped)d dropped), "
          "%(other)d other -> --other" % counts)
    if counts["nonmissense_dropped"]:
        print("# dropped %d non-missense variants from the missense bins (--nonmissense drop)"
              % counts["nonmissense_dropped"])
    if counts["dropped_anno"]:
        print("# dropped %d variants by --drop-anno (%d regions emptied and omitted)"
              % (counts["dropped_anno"], counts["emptied_regions"]))
    if counts["unparsed_ids"]:
        print("# WARNING: %d missense IDs did not parse into coordinates and "
              "could not be joined" % counts["unparsed_ids"])
    if counts["missense"]:
        frac = 100.0 * counts["missense_scored"] / counts["missense"]
        print("# %s coverage of missense variants: %.1f%%"
              % (a.score_column or "AlphaMissense", frac))
        if frac < 50.0:
            print("# WARNING: under half the missense variants got a real "
                  "AlphaMissense score. Check the build (AlphaMissense is hg38) "
                  "and the REF/ALT orientation of the group file's IDs before "
                  "reading anything into a FlexRV result.")
    print("#\n# Run it with, e.g.:")
    print("#   --groupFile=%s --flexRV_score=%s --flexRV_maxMAF=0.001 \\\n"
          "#   --flexRV_lofAnno=%s --annotation_in_groupTest=%s --r.corr=1"
          % (a.out, a.name, a.lof_anno,
             ":".join([s for s in (a.lof_anno + "," + a.missense_anno).split(",") if s])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
