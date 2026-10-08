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
ID_FORMATS = ("colon", "underscore_slash", "dash")


def detect_id_format(vid):
    """Return one of ID_FORMATS, or None if the ID carries no coordinates."""
    if vid.count(":") >= 3:
        return "colon"
    if "_" in vid and "/" in vid and ":" in vid:
        return "underscore_slash"
    if vid.count("-") >= 3:
        return "dash"
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
    with smart_open(path) as fh:
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


def read_anno_table(path, id_col="ID", gene_col="GENE", tx_col="TRANSCRIPT",
                    mane_col="MANE_SELECT", csq_col="CSQ", score_col=None,
                    zero_is_missing=True):
    """Per-variant transcript context from a BRaVa long-form annotation table.

    That table is the SOURCE of the group file -- one row per (gene, variant),
    carrying the transcript each consequence was called against, the MANE
    Select accession when the transcript has one, and the raw VEP consequence.
    It is what makes "the MANE Select transcript, or canonical where there is
    none" a checkable claim rather than an assumption: AlphaMissense names the
    transcript it scored, and this names the transcript the annotation used.

    Returns {(gene, variant_id): (transcript, mane_accession_or_None, csq)}.
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
          provenance=None, log=sys.stderr):
    """Write `group_out` = `group_in` with a score line per region.

    `fallback` is what --prefer-anno-transcript does when AlphaMissense has no
    row on the annotation's transcript: "safe" (see pick()), or "max", the rule
    before 2026-10, kept only so files built with it can be reproduced.
    `provenance`, a path, gets one line per missense (gene, variant): the
    source of its score, the transcript(s) and the value before imputation.

    Returns a dict of counts -- the report is part of the output, because
    "how many missense variants got a real AlphaMissense score" is the number
    that decides whether the run means anything.
    """
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
    anno_ctx = (read_anno_table(anno_table_path, score_col=score_column,
                                zero_is_missing=zim)
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
            zim, log)

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
    # to a DIFFERENT gene at this coordinate, or (with --mane) the MANE summary
    # puts it in a different gene. No rule here picks a row by its score.
    pick_stats = dict.fromkeys(
        ("by_transcript", "fallback_canonical", "fallback_isoform",
         "fallback_canonical_averaged", "fallback_isoform_averaged",
         "missing_no_row", "missing_excluded",
         "fallback_max"), 0)
    mane_genes = read_mane_genes(mane_path) if (mane_path and prefer_tx) else None
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
            for score, tx, src, _ in cands:
                if tx.split(".")[0] == w:
                    pick_stats["by_transcript"] += 1
                    return (score, tx, True, "matched")
        if fallback == "max":
            best = max(cands, key=lambda c: c[0])
            pick_stats["fallback_max"] += 1
            return (best[0], best[1], False, "max")
        others = set(t for g, t in tx_at.get(vid, {}).items() if g != gene)

        def gene_of(tx):
            return mane_genes.get(tx) if mane_genes is not None else None

        ok = []
        for score, tx, src, up in cands:
            t = tx.split(".")[0]
            mg = gene_of(t)
            if t in others or (mg is not None and gene not in mg):
                continue
            ok.append((score, tx, src, mg is not None, up))
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
        return None

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
            eligible = snv and (csq is None or csq == "missense_variant")
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

    # --- per-gene score vectors
    drop = set(drop_anno)
    counts = {"genes": len(order), "variants": 0, "lof": 0, "missense": 0,
              "dropped_anno": 0, "emptied_regions": 0,
              "missense_scored": 0, "missense_imputed": 0, "missense_dropped": 0,
              "other": 0, "unparsed_ids": unparsed, "am_keys_matched": len(hits),
              "am_keys_wanted": len(wanted), "am_collisions": collisions,
              "from_isoform": from_isoform, "mane_hits": mane_hits}
    counts.update({"diag_" + k: v for k, v in counts_diag.items()})

    prov = open(provenance, "w") if provenance else None
    if prov is not None:
        prov.write("GENE\tID\tSOURCE\tTRANSCRIPT\tAM\n")
    out_lines = []
    for gene in order:
        r = regions[gene]
        raw_scores, kinds = [], []
        for vid, anno in zip(r["var"], r["anno"]):
            counts["variants"] += 1
            if anno in lof:
                raw_scores.append(1.0)
                kinds.append("lof")
                counts["lof"] += 1
            elif anno in mis:
                counts["missense"] += 1
                key = parse_variant_id(vid, id_format)
                hit = pick(gene, vid, key)
                if prov is not None:
                    prov.write("%s\t%s\t%s\t%s\t%s\n" % (
                        gene, vid, "missing" if hit is None else hit[3],
                        "." if hit is None else hit[1],
                        "NA" if hit is None else "%.6f" % hit[0]))
                if hit is None:
                    raw_scores.append(None)
                    kinds.append("missense_missing")
                else:
                    raw_scores.append(hit[0])
                    kinds.append("missense")
                    counts["missense_scored"] += 1
            else:
                raw_scores.append(float(other))
                kinds.append("other")
                counts["other"] += 1

        # resolve missing missense
        scored = [s for s, k in zip(raw_scores, kinds) if k == "missense"]
        if missing == "gene-mean":
            fill = sum(scored) / len(scored) if scored else 0.0
        elif missing == "drop":
            fill = None
        else:
            fill = float(missing)

        keep = [True] * len(raw_scores)
        for i, k in enumerate(kinds):
            if k != "missense_missing":
                continue
            if fill is None:
                keep[i] = False
                counts["missense_dropped"] += 1
            else:
                raw_scores[i] = fill
                counts["missense_imputed"] += 1

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
                out_lines.append("%s score:%s %s" % (
                    gene, score_name, " ".join("%.4f" % s for s in final_score)))
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
            if mane_genes is None:
                print("# (no --mane: a row is excluded only when the annotation "
                      "assigns its transcript to another gene at that coordinate)",
                      file=log)

    with open(group_out, "w") as fh:
        fh.write("\n".join(out_lines) + "\n")
    return counts


def _build_from_column(regions, order, lof, mis, drop_anno, group_out,
                       score_name, anno_ctx, score_column, score_transform,
                       missing, other, counts, zero_missing, log):
    """Score every region from an annotation-table COLUMN.

    Same construction as the AlphaMissense path -- LoF 1.0 by annotation, the
    score classes from the source, everything else --other -- but the value is
    read rather than joined, and the per-gene transform runs before the
    missing-value policy so that an imputed value is on the TRANSFORMED scale
    (imputing a raw CADD PHRED and then ranking it would let the fill change
    every other variant's rank).
    """
    drop = set(drop_anno)
    transform = score_transform or SCORE_COLUMNS.get(score_column, ("none",))[0]
    counts = dict(counts)
    counts.update({"variants": 0, "lof": 0, "missense": 0, "missense_scored": 0,
                   "missense_imputed": 0, "missense_dropped": 0, "other": 0,
                   "dropped_anno": 0, "emptied_regions": 0, "unparsed_ids": 0,
                   "am_keys_matched": 0, "am_keys_wanted": 0, "am_collisions": 0,
                   "from_isoform": 0, "mane_hits": None})
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
            elif anno in mis:
                counts["missense"] += 1
                val = anno_ctx.get((gene, vid), (None, None, None, None))[3]
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
                out_lines.append("%s score:%s %s" % (
                    gene, score_name, " ".join("%.4f" % v for v in fs)))
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
    fx = {"S1": 5001, "S2": 5002, "S3": 5003, "S4": 5004, "S5": 5005,
          "S6": 5006, "S7": 5007, "M1": 5011, "M2": 5012, "M3": 5013}
    vid = dict((k, "20:%d:A:G" % p) for k, p in fx.items())
    p_ids = [vid[k] for k in ("S1", "S2", "S3", "S4", "S5", "S6", "S7",
                              "M1", "M2", "M3")]
    q_ids = [vid[k] for k in ("S1", "S2", "S4")]
    gf3 = os.path.join(tmp, "group_shared.txt")
    with open(gf3, "w") as fh:
        for g, ids in (("GENEP", p_ids), ("GENEQ", q_ids)):
            lof_id = "20:%d:A:T" % (5090 if g == "GENEP" else 5091)
            fh.write("%s var %s %s\n" % (g, " ".join(ids), lof_id))
            fh.write("%s anno %s pLoF\n" % (g, " ".join(["missense"] * len(ids))))
    anno3 = os.path.join(tmp, "anno_shared.tsv")
    with open(anno3, "w") as fh:
        fh.write("\t".join(["ID", "GENE", "LOF", "REVEL_SCORE", "CADD_PHRED",
                            "CSQ", "TRANSCRIPT", "MANE_SELECT", "CANONICAL",
                            "BIOTYPE", "DS_MAX", "ANNOTATION"]) + "\n")
        for g, ids, tx in (("GENEP", p_ids, "ENST0000000TP"),
                           ("GENEQ", q_ids, "ENST0000000TQ")):
            for v in ids:
                fh.write("\t".join([v, g, "NA", "0", "0", "missense_variant", tx,
                                    "NM_1", "YES", "protein_coding", "0",
                                    "missense"]) + "\n")
    canon_rows = [("S1", "ENST000000TPC.1", 0.30, "P0"), ("S3", "ENST000000TPC.1", 0.70, "P0"),
                  ("S6", "ENST0000000TX.1", 0.10, "P1"), ("S6", "ENST0000000TY.1", 0.85, "P2"),
                  ("S7", "ENST000000TZ1.1", 0.20, "P3"), ("S7", "ENST000000TZ2.1", 0.60, "P3-2")]
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
    check("fallback: no non-imputed score moves UP against the legacy max",
          all(new3[P][vid[k]] <= old3[P][vid[k]] + 5e-5 for k in ("S1", "S3", "S4", "S5")))

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
    ap.add_argument("--provenance", help="write GENE/ID/SOURCE/TRANSCRIPT/AM per "
                    "missense variant: which row each score came from")
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
                   isoforms_path=a.isoforms, mane_path=a.mane,
                   replace_existing=a.replace_existing,
                   anno_table_path=a.anno_table,
                   prefer_anno_transcript=a.prefer_anno_transcript,
                   drop_anno=[x for x in a.drop_anno.split(",") if x],
                   score_column=a.score_column,
                   score_transform=a.score_transform, zero_is_missing=a.zim,
                   fallback="max" if a.legacy_max_fallback else "safe",
                   provenance=a.provenance)

    print("# wrote %s" % a.out)
    print("# %(genes)d regions, %(variants)d variants: %(lof)d lof -> 1.0, "
          "%(missense)d missense (%(missense_scored)d scored, "
          "%(missense_imputed)d imputed, %(missense_dropped)d dropped), "
          "%(other)d other -> --other" % counts)
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
