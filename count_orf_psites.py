#!/usr/bin/env python3
"""Count strand-specific P-sites in the three frames of spliced ORFs.

Python 3.9+; standard library only. Run with --help for arguments and an example.
"""

import argparse
from array import array
from collections import Counter, defaultdict
import csv
import gzip
import math
from pathlib import Path
import re
import sys

METADATA = ('orf_id', 'orf_biotype', 'transcript_id', 'gene_id', 'gene_name')
BIN_SIZE = 16384
ATTRIBUTE = re.compile(r'(\w+)\s+"([^"\r\n]*)"')


def open_text(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    return opener(path, 'rt', encoding='utf-8-sig')


def read_gtf(path):
    """Group CDS by (orf_id, transcript_id) and index exon-aware frame offsets."""
    orfs = {}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            where = f'{path}:{line_number}'
            if len(fields) != 9:
                raise ValueError(f'{where}: expected nine tab-separated GTF fields')
            if fields[2] != 'CDS':
                continue
            attrs = dict(ATTRIBUTE.findall(fields[8]))
            if not attrs.get('orf_id') or not attrs.get('transcript_id'):
                raise ValueError(f'{where}: CDS needs orf_id and transcript_id')
            start, end = int(fields[3]) - 1, int(fields[4])
            strand = fields[6]
            if start < 0 or end <= start or strand not in ('+', '-'):
                raise ValueError(f'{where}: invalid CDS coordinates or strand')
            key = (attrs['orf_id'], attrs['transcript_id'])
            metadata = tuple(attrs.get(name, '.') for name in METADATA)
            if key not in orfs:
                orfs[key] = (metadata, fields[0], strand, set())
            old_metadata, chrom, old_strand, segments = orfs[key]
            if (metadata, fields[0], strand) != (old_metadata, chrom, old_strand):
                raise ValueError(f'{where}: inconsistent annotation for {key}')
            segments.add((start, end))  # Exact duplicate records count once.
    if not orfs:
        raise ValueError(f'{path}: no CDS records with ORF annotations found')
    rows = []
    index = defaultdict(lambda: defaultdict(list))
    for row_number, (key, (metadata, chrom, strand, segments)) in enumerate(orfs.items()):
        ordered = sorted(segments)
        if any(left[1] > right[0] for left, right in zip(ordered, ordered[1:])):
            raise ValueError(f'{path}: overlapping CDS segments within {key}')
        if strand == '-':
            ordered.reverse()
        offset = 0
        for start, end in ordered:
            block = (start, end, offset % 3, row_number)
            for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
                index[(chrom, strand)][bin_id].append(block)
            offset += end - start
        rows.append(metadata)
    repeated = sum(n > 1 for n in Counter(row[0] for row in rows).values())
    if repeated:
        print(f'Note: {repeated} ORF IDs occur in multiple transcripts; '
              'keeping separate (orf_id, transcript_id) rows.', file=sys.stderr)
    return rows, index


def read_bedgraph_list(path):
    """Read whitespace-separated: bedgraph_path tissue plus|minus."""
    samples, seen = {}, set()
    aliases = {'plus': '+', '+': '+', 'minus': '-', '-': '-'}
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split()
            where = f'{path}:{line_number}'
            if len(fields) != 3 or fields[2].lower() not in aliases:
                raise ValueError(f'{where}: expected bedgraph_path tissue plus|minus')
            filename, tissue, strand_text = fields
            if tissue in METADATA:
                raise ValueError(f'{where}: tissue name conflicts with metadata column')
            bedgraph = Path(filename).expanduser()
            if not bedgraph.is_absolute():
                bedgraph = Path(path).resolve().parent / bedgraph
            if not bedgraph.is_file():
                raise ValueError(f'{where}: bedGraph not found: {bedgraph}')
            strand = aliases[strand_text.lower()]
            identity = (bedgraph.resolve(), tissue, strand)
            if identity in seen:
                raise ValueError(f'{where}: duplicate bedGraph/tissue/strand entry')
            seen.add(identity)
            samples.setdefault(tissue, []).append((bedgraph, strand))
    if not samples:
        raise ValueError(f'{path}: empty bedGraph list')
    for tissue, entries in samples.items():
        if {strand for _, strand in entries} != {'+', '-'}:
            print(f'Warning: {tissue} has only one strand; the other gets zero counts.',
                  file=sys.stderr)
    return samples


def count_bedgraph(path, strand, index, counts):
    """Count residues modulo three without expanding bedGraph intervals to bases."""
    matched, chromosomes = False, set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields or fields[0].startswith('#') or fields[0] in ('track', 'browser'):
                continue
            where = f'{path}:{line_number}'
            if len(fields) != 4:
                raise ValueError(f'{where}: expected four bedGraph fields')
            chrom, start_text, end_text, value_text = fields
            start, end, value = int(start_text), int(end_text), float(value_text)
            if start < 0 or end <= start or not math.isfinite(value) or value < 0:
                raise ValueError(f'{where}: invalid interval or P-site count (must be nonnegative)')
            chromosomes.add(chrom)
            bins = index.get((chrom, strand))
            if bins is None or value == 0:
                continue
            for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
                for exon_start, exon_end, offset, row in bins.get(bin_id, ()):
                    lo, hi = max(start, exon_start), min(end, exon_end)
                    # A segment spanning bins contributes only once per interval.
                    if lo >= hi or bin_id != lo // BIN_SIZE:
                        continue
                    matched = True
                    first = (offset + (lo - exon_start if strand == '+' else exon_end - hi)) % 3
                    length = hi - lo
                    for frame in range(3):
                        bases = (length + 2 - (frame - first) % 3) // 3
                        counts[3 * row + frame] += bases * value
    missing = chromosomes - {chrom for chrom, s in index if s == strand}
    if missing:
        print(f'Note: {path.name}: {len(missing)} chromosome(s) have no {strand} CDS '
              f"in the GTF (e.g. {', '.join(sorted(missing)[:5])}).", file=sys.stderr)
    if not matched:
        print(f'Warning: {path}: no positive P-sites overlapped matching-strand CDS; '
              'check chromosome names and coordinates.', file=sys.stderr)


def report_cds_frames(tissue, rows, counts):
    """Report pooled raw frame counts for rows with orf_biotype exactly CDS."""
    cds_rows = [i for i, row in enumerate(rows) if row[1] == 'CDS']
    totals = [math.fsum(counts[3 * i + frame] for i in cds_rows)
              for frame in range(3)]
    total = math.fsum(totals)
    percentages = ', '.join(
        f'p{frame}={100 * value / total:.2f}%' if total else f'p{frame}=NA'
        for frame, value in enumerate(totals))
    raw = ', '.join(f'p{frame}={value:.12g}' for frame, value in enumerate(totals))
    note = '' if total else '; no CDS counts to calculate percentages'
    print(f'{tissue}: orf_biotype=CDS frame percentages: {percentages} '
          f'(raw counts: {raw}; total={total:.12g}{note})',
          file=sys.stderr, flush=True)


def write_tables(prefix, rows, samples, counts, target_total):
    """Use one normalization factor per tissue, shared across all three frames."""
    factors = {}
    for tissue in samples:
        total = math.fsum(counts[tissue])
        factors[tissue] = target_total / total if total else 0.0
        if not total:
            print(f'Warning: {tissue}: zero total; normalized column stays zero.', file=sys.stderr)
        print(f'{tissue}: assigned total={total:.12g}; normalization factor='
              f'{factors[tissue]:.12g}', file=sys.stderr)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    for frame in range(3):
        for normalized in (False, True):
            suffix = '.normalized' if normalized else ''
            path = f'{prefix}.p{frame}{suffix}.tsv'
            with open(path, 'w', newline='', encoding='utf-8') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow((*METADATA, *samples))
                for row_number, metadata in enumerate(rows):
                    values = []
                    for tissue in samples:
                        value = counts[tissue][3 * row_number + frame]
                        if normalized:
                            value *= factors[tissue]
                        values.append(format(value, '.15g'))
                    writer.writerow((*metadata, *values))
            print(f'Wrote {path}', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=r'''Example (from this script's directory):
  python count_orf_psites.py \
    --gtf ../human_nicos_20231005_pooled/collapsed/nicos_20231005.orfs.gtf \
    --bedgraphs checked_tissues.txt --out-prefix nicos_psites

Input list: one whitespace-separated 'bedgraph_path tissue plus|minus' per line.
The third column also accepts + or -. Blank lines and # comments are ignored.
Paths without spaces are required; relative bedGraph paths resolve relative to
this list file. Multiple distinct files for one tissue/strand are added together.
GTF, bedGraph and list inputs may be gzip-compressed (.gz). Sorting is not needed.

Coordinates and frames:
  Only GTF CDS features are counted. GTF [start,end] is 1-based inclusive;
  bedGraph [start,end) is 0-based, end-exclusive. Chromosome names must match
  exactly (e.g. '1' differs from 'chr1'). Each bedGraph value applies to EVERY
  base in its interval. Counts must be nonnegative; overlapping bedGraph records
  contribute additively. Each file is matched only to its listed strand.
  The first CDS base in the ORF's 5'-to-3' direction is p0, followed by p1,p2.
  Frames continue across CDS junctions, excluding introns. On minus strand,
  segments and bases are traversed in descending genomic order. GTF phase is
  not used: supplied CDS segments must start at the ORF's first nucleotide.

Rows and normalization:
  One row per (orf_id, transcript_id), in GTF encounter order; shared orf_id
  values can have distinct transcript-specific CDS structures. Exact duplicate
  CDS segments within a row count once. Overlaps between different ORFs count
  toward each ORF independently. Missing optional annotations are written as '.'.
  Tissue columns follow first appearance in the list; unobserved ORFs get zero.
  Each tissue uses ONE factor: target_total / sum(all ORFs, all three frames).
  Thus p0+p1+p2 together sum to target_total per nonzero tissue, preserving frame
  ratios. The denominator counts ORF assignments (including assignments to
  overlapping ORFs), not all P-sites in the input library. Zero-total tissues
  stay zero with a warning. No length normalization is applied.

Progress report (stderr):
  After all files for each tissue are counted, report p0/p1/p2 percentages from
  summed raw counts over rows with orf_biotype exactly CDS, combining strands.
  Each percentage is 100 * CDS frame total / (CDS p0 + CDS p1 + CDS p2).
  These are pooled count percentages, not averages of per-ORF percentages.
  Raw frame totals are included; percentages are NA when the CDS total is zero.

Outputs (existing files with these names are overwritten):
  PREFIX.p0.tsv, PREFIX.p1.tsv, PREFIX.p2.tsv
  PREFIX.p0.normalized.tsv, PREFIX.p1.normalized.tsv, PREFIX.p2.normalized.tsv
''')
    parser.add_argument('--gtf', required=True, type=Path,
                        help='ORF GTF: CDS features with orf_id/transcript_id attributes')
    parser.add_argument('--bedgraphs', required=True, type=Path,
                        help='three-column bedGraph list: path, tissue, strand')
    parser.add_argument('--out-prefix', required=True,
                        help='output path prefix, e.g. results/nicos_psites')
    parser.add_argument('--target-total', type=float, default=1_000_000,
                        help='normalized total across all ORFs and frames per tissue (default: 1000000)')
    args = parser.parse_args()
    if not math.isfinite(args.target_total) or args.target_total <= 0:
        parser.error('--target-total must be finite and greater than zero')
    try:
        samples = read_bedgraph_list(args.bedgraphs)
        rows, index = read_gtf(args.gtf)
        print(f'Loaded {len(rows)} ORF/transcript rows and {len(samples)} tissues.', file=sys.stderr)
        counts = {}
        for tissue, entries in samples.items():
            counts[tissue] = array('d', [0.0]) * (3 * len(rows))
            for path, strand in entries:
                print(f'Counting {tissue} ({strand}): {path}', file=sys.stderr)
                count_bedgraph(path, strand, index, counts[tissue])
            report_cds_frames(tissue, rows, counts[tissue])
        write_tables(args.out_prefix, rows, samples, counts, args.target_total)
    except (OSError, ValueError, OverflowError) as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()
