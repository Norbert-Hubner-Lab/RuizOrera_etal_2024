#!/usr/bin/env python3
"""Count strand-specific ribosome P-sites or RNA-seq reads in spliced ORFs.

Python 3.9+; PDF plots require NumPy and Matplotlib (pip install matplotlib).
Run with --help for arguments and an example.
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
SPECIFICITY = ('max_tissue_translation', 'tau_tissue_translation')
RNA_SPECIFICITY = ('max_tissue_expression', 'tau_tissue_expression')
BIN_SIZE = 16384
ATTRIBUTE = re.compile(r'(\w+)\s+"([^"\r\n]*)"')


def open_text(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    return opener(path, 'rt', encoding='utf-8-sig')


def read_gtf(path, reads='ribo'):
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
    rows, lengths = [], []
    index = defaultdict(lambda: defaultdict(list))
    for row_number, (key, (metadata, chrom, strand, segments)) in enumerate(orfs.items()):
        ordered = sorted(segments)
        if any(left[1] > right[0] for left, right in zip(ordered, ordered[1:])):
            raise ValueError(f'{path}: overlapping CDS segments within {key}')
        if strand == '-':
            ordered.reverse()
        offset = 0
        for start, end in ordered:
            block = (start, end, offset, row_number)
            for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
                index[(chrom, strand)][bin_id].append(block)
            offset += end - start
        rows.append(metadata)
        lengths.append(offset)
    repeated = sum(n > 1 for n in Counter(row[0] for row in rows).values())
    if repeated:
        print(f'Note: {repeated} ORF IDs occur in multiple transcripts; '
              'keeping separate (orf_id, transcript_id) rows.', file=sys.stderr)
    incomplete = sum(length % 3 != 0 for length in lengths)
    if incomplete and reads == 'ribo':
        print(f'Warning: {incomplete} ORFs have CDS length not divisible by three; '
              'coverage uses the actual number of p0 positions (ceil(length/3)).',
              file=sys.stderr)
    return rows, index, lengths


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
            if tissue in METADATA + SPECIFICITY + RNA_SPECIFICITY or ';' in tissue or tissue == 'NA':
                raise ValueError(f'{where}: tissue name conflicts with output columns '
                                 'or reserved specificity notation (NA or ;)')
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


def read_bam_list(path):
    """Read whitespace-separated BAM path and tissue; pool files per tissue."""
    samples, seen = {}, set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            fields = line.split()
            where = f'{path}:{line_number}'
            if len(fields) != 2:
                raise ValueError(f'{where}: expected bam_path tissue')
            filename, tissue = fields
            if tissue in METADATA + SPECIFICITY + RNA_SPECIFICITY or ';' in tissue or tissue == 'NA':
                raise ValueError(f'{where}: reserved tissue name')
            bam = Path(filename).expanduser()
            if not bam.is_absolute():
                bam = Path(path).resolve().parent / bam
            if not bam.is_file():
                raise ValueError(f'{where}: BAM not found: {bam}')
            identity = (bam.resolve(), tissue)
            if identity in seen:
                raise ValueError(f'{where}: duplicate BAM/tissue entry')
            seen.add(identity)
            samples.setdefault(tissue, []).append(bam)
    if not samples:
        raise ValueError(f'{path}: empty BAM list')
    return samples


def load_bam_reader():
    try:
        import pysam
    except ImportError as error:
        raise ValueError('RNA mode requires pysam; install with python -m pip install pysam') from error
    return pysam


def count_bam(path, index, counts, pysam, strandedness='unstranded'):
    """Count each primary alignment once per overlapping ORF, across all CDS bases."""
    assigned = 0
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            strand = '-' if read.is_reverse else '+'
            if strandedness != 'unstranded':
                # Forward: read 1/single read agrees with transcript; read 2 is opposite.
                if (read.is_paired and read.is_read2) != (strandedness == 'reverse'):
                    strand = '+' if strand == '-' else '-'
                strands = (strand,)
            else:
                strands = ('+', '-')
            hits = set()
            for start, end in read.get_blocks():
                for target_strand in strands:
                    bins = index.get((read.reference_name, target_strand), {})
                    for bin_id in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
                        for exon_start, exon_end, _, row in bins.get(bin_id, ()):
                            if start < exon_end and end > exon_start:
                                hits.add(row)
            for row in hits:
                counts[row] += 1
            assigned += len(hits)
    if not assigned:
        print(f'Warning: {path}: no RNA reads overlapped CDS; check chromosome names '
              'and strandedness.', file=sys.stderr)


def count_bedgraph(path, strand, index, counts, covered=None):
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
                    transcript_start = offset + (lo - exon_start if strand == '+' else exon_end - hi)
                    first = transcript_start % 3
                    length = hi - lo
                    for frame in range(3):
                        bases = (length + 2 - (frame - first) % 3) // 3
                        counts[3 * row + frame] += bases * value
                    if covered is not None:
                        # Pool positive p0 positions across overlapping records/files.
                        begin = (transcript_start + 2) // 3
                        stop = (transcript_start + length + 2) // 3
                        covered[row][begin:stop] = b'\x01' * (stop - begin)
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


def tissue_specificity(values, tissues):
    """Tau on linear normalized counts; retain every exactly tied maximum."""
    maximum = max(values)
    if maximum == 0:
        return (), math.nan
    winners = tuple(tissue for tissue, value in zip(tissues, values) if value == maximum)
    tau = (math.fsum(1 - value / maximum for value in values) / (len(values) - 1)
           if len(values) > 1 else math.nan)
    return winners, min(1.0, max(0.0, tau)) if math.isfinite(tau) else tau


def format_metric(value):
    return format(value, '.15g') if math.isfinite(value) else 'NA'


def write_tables(prefix, rows, samples, counts, target_total, reads='ribo'):
    """Normalize assigned counts per tissue; ribosome frames share one factor."""
    stride = 3 if reads == 'ribo' else 1
    specificity_columns = SPECIFICITY if reads == 'ribo' else RNA_SPECIFICITY
    factors = {}
    for tissue in samples:
        total = math.fsum(counts[tissue])
        factors[tissue] = target_total / total if total else 0.0
        if not total:
            print(f'Warning: {tissue}: zero total; normalized column stays zero.', file=sys.stderr)
        print(f'{tissue}: assigned total={total:.12g}; normalization factor='
              f'{factors[tissue]:.12g}', file=sys.stderr)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    p0_specificity = []
    for frame in range(stride):
        for normalized in (False, True):
            suffix = '.normalized' if normalized else ''
            component = f'p{frame}' if reads == 'ribo' else 'counts'
            path = f'{prefix}.{component}{suffix}.tsv'
            with open(path, 'w', newline='', encoding='utf-8') as handle:
                writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                writer.writerow((*METADATA, *samples, *(specificity_columns if normalized else ())))
                for row_number, metadata in enumerate(rows):
                    values = [counts[tissue][stride * row_number + frame] *
                              (factors[tissue] if normalized else 1.0) for tissue in samples]
                    extra = ()
                    if normalized:
                        winners, tau = tissue_specificity(values, samples)
                        extra = (';'.join(winners) if winners else 'NA', format_metric(tau))
                        if frame == 0:
                            p0_specificity.append((winners, tau))
                    writer.writerow((*metadata, *(format_metric(v) for v in values), *extra))
            print(f'Wrote {path}', file=sys.stderr)
    return p0_specificity


def write_quality_table(prefix, rows, samples, counts, coverage):
    """Write p0/all-frame periodicity and distinct-p0-position coverage."""
    periodicity = {}
    for tissue in samples:
        values = array('d')
        for row in range(len(rows)):
            total = math.fsum(counts[tissue][3 * row:3 * row + 3])
            values.append(counts[tissue][3 * row] / total if total else math.nan)
        periodicity[tissue] = values
    path = f'{prefix}.periodicity_coverage.tsv'
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow((*METADATA, *(f'{tissue}_{metric}' for tissue in samples
                                      for metric in ('periodicity', 'coverage'))))
        for row, metadata in enumerate(rows):
            writer.writerow((*metadata, *(format_metric(metric[tissue][row])
                                          for tissue in samples
                                          for metric in (periodicity, coverage))))
    print(f'Wrote {path}', file=sys.stderr)
    return periodicity


def load_plotting():
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as error:
        raise ValueError('PDF plots require NumPy and Matplotlib; install with '
                         'python -m pip install matplotlib, or use --skip-plots') from error
    return np, plt


def plot_density_panel(ax, values, np):
    """Binned Gaussian KDE, reflected at 0 and 1; no density for point masses."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.text(0.97, 0.95, f'n={len(values):,}', transform=ax.transAxes,
            ha='right', va='top', fontsize=8)
    if not len(values):
        ax.text(0.5, 0.5, 'No defined values', transform=ax.transAxes,
                ha='center', fontsize=8)
        return
    if np.all(values == values[0]):
        ax.axvline(values[0], color='#2878a5', linewidth=2, clip_on=False)
        ax.text(0.5, 0.5, f'All = {values[0]:.3g}', transform=ax.transAxes,
                ha='center', fontsize=8)
        return
    bins = 256
    hist, edges = np.histogram(values, bins=bins, range=(0, 1), density=True)
    centers = (edges[:-1] + edges[1:]) / 2
    # Scott bandwidth with a bin-width floor. Binning keeps large panels fast.
    bandwidth = max(float(np.std(values, ddof=1)) * len(values) ** (-0.2), 1 / bins)
    radius = min(bins, int(math.ceil(4 * bandwidth * bins)))
    offsets = np.arange(-radius, radius + 1) / bins
    kernel = np.exp(-0.5 * (offsets / bandwidth) ** 2)
    kernel /= kernel.sum()
    padded = np.concatenate((hist[::-1], hist, hist[::-1]))
    density = np.convolve(padded, kernel, mode='same')[bins:2 * bins]
    x = np.concatenate(([0], centers, [1]))
    y = np.concatenate(([density[0]], density, [density[-1]]))
    ax.plot(x, y, color='#2878a5', linewidth=1.2)
    ax.fill_between(x, y, color='#2878a5', alpha=0.25)
    ax.set_ylim(bottom=0)


def write_density_plots(prefix, rows, samples, periodicity, coverage, specificity, reads='ribo'):
    np, plt = load_plotting()
    tissues = list(samples)
    biotypes = list(dict.fromkeys(row[1] for row in rows))
    by_biotype = {biotype: [] for biotype in biotypes}
    tau_groups = defaultdict(list)
    for i, row in enumerate(rows):
        by_biotype[row[1]].append(i)
        _, tau = specificity[i]
        if math.isfinite(tau):
            tau_groups[row[1]].append(tau)
    metrics = (
            ('periodicity', 'Periodicity (p0 / all frames)', periodicity),
            ('coverage', 'Coverage (covered p0 positions / all p0 positions)', coverage),
            ('tau_tissue_translation', 'Tissue specificity (tau, normalized p0)', None))
    if reads == 'rna':
        metrics = (('tau_tissue_expression', 'Tissue specificity (tau, normalized RNA counts)', None),)
    for metric, label, data in metrics:
        columns = ['All tissues'] if data is None else tissues
        fig, axes = plt.subplots(len(biotypes), len(columns), squeeze=False,
                                 figsize=(max(4, 2.6 * len(columns)), 2.1 * len(biotypes) + 1),
                                 sharex=True)
        for r, biotype in enumerate(biotypes):
            for c, tissue in enumerate(columns):
                ax = axes[r, c]
                values = ([data[tissue][i] for i in by_biotype[biotype]] if data is not None
                          else tau_groups[biotype])
                plot_density_panel(ax, values, np)
                if r == 0:
                    ax.set_title(tissue, fontsize=10)
                if c == 0:
                    ax.set_ylabel(f'{biotype}\nDensity', fontsize=9)
                if r == len(biotypes) - 1:
                    ax.set_xlabel('Tau' if data is None else metric.capitalize())
        fig.suptitle(label, fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        path = f'{prefix}.{metric}.density.pdf'
        fig.savefig(path, bbox_inches='tight')
        plt.close(fig)
        print(f'Wrote {path}', file=sys.stderr)


def count_tissue_enrichment(rows, specificity):
    """Count nested strict tau thresholds in each maximum-translation tissue."""
    enriched, specific = Counter(), Counter()
    for metadata, (winners, tau) in zip(rows, specificity):
        if not math.isfinite(tau):
            continue
        for tissue in winners:
            key = (metadata[1], tissue)
            if tau > 0.8:
                enriched[key] += 1
            if tau > 0.95:
                specific[key] += 1
    return enriched, specific


def write_enrichment_barplot(prefix, rows, samples, specificity, reads='ribo'):
    np, plt = load_plotting()
    from matplotlib.ticker import MaxNLocator

    tissues = list(samples)
    biotypes = list(dict.fromkeys(row[1] for row in rows))
    enriched, specific = count_tissue_enrichment(rows, specificity)
    fig, axes = plt.subplots(len(biotypes), 1, squeeze=False, sharex=True,
                             figsize=(max(6, 1.2 * len(tissues)), 2.5 * len(biotypes) + 1.5))
    x = np.arange(len(tissues))
    for r, biotype in enumerate(biotypes):
        ax = axes[r, 0]
        for offset, counts, label, color in (
                (-0.2, enriched, 'Tissue-enriched (tau > 0.8)', '#2878a5'),
                (0.2, specific, 'Tissue-specific (tau > 0.95)', '#e58b38')):
            bars = ax.bar(x + offset, [counts[(biotype, tissue)] for tissue in tissues],
                          width=0.4, label=label, color=color)
            ax.bar_label(bars, padding=3, fontsize=8)
        ax.set_ylabel(f'{biotype}\nNumber of ORFs')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        maximum = max((enriched[(biotype, tissue)] for tissue in tissues), default=0)
        ax.set_ylim(0, max(1, maximum * 1.2))
        ax.set_xticks(x)
        ax.set_xticklabels(tissues, rotation=45, ha='right')
    axes[-1, 0].set_xlabel('Tissue of maximum normalized ' +
                           ('p0 translation' if reads == 'ribo' else 'RNA count'))
    fig.suptitle('Tissue-enriched and tissue-specific ORFs by biotype', fontsize=12)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.96), ncol=2)
    fig.text(0.5, 0.01, 'Specific ORFs are included in enriched counts; tied maxima count in each tied tissue.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 0.91))
    path = f'{prefix}.tissue_enrichment.barplot.pdf'
    fig.savefig(path, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote {path}', file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=r'''Example (from this script's directory):
  python count_orf_psites.py \
    --gtf ../human_nicos_20231005_pooled/collapsed/nicos_20231005.orfs.gtf \
    --bedgraphs checked_tissues.txt --out-prefix nicos_psites

RNA-seq mode:
  python count_orf_psites.py --reads rna --gtf orfs.gtf --input rna_bams.txt \
    --out-prefix rna_counts
  RNA list: one whitespace-separated 'bam_path tissue' per line. Relative paths
  resolve against the list directory; multiple BAMs per tissue are pooled.
  Requires pysam (python -m pip install pysam). BAM sorting/indexing is not required.
  Each primary, mapped, QC-passing alignment counts once per overlapping ORF,
  regardless of reading frame. Paired mates count separately. Duplicate-marked
  reads are retained; secondary and supplementary alignments are excluded.
  Only aligned bases count: skipped introns and deletions do not create overlaps.
  Shared ORFs receive independent assignments, including for multimapping primary
  alignments. Default is unstranded. --rna-strandedness forward means read 1 or
  single-end reads agree with transcript strand; reverse means they oppose it.
  Read 2 uses the opposite orientation to read 1.
  Normalization scales total ORF-assigned reads per tissue to --target-total,
  without length normalization. Tau and enrichment use normalized RNA counts.
  RNA outputs: PREFIX.counts.tsv, PREFIX.counts.normalized.tsv (adds
  max_tissue_expression and tau_tissue_expression),
  PREFIX.tau_tissue_expression.density.pdf, PREFIX.tissue_enrichment.barplot.pdf.
  No frames, periodicity, or coverage are computed in RNA mode.

Ribo mode (default; remaining details describe ribo outputs):
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

Tissue specificity and translation quality:
  Each normalized frame table adds max_tissue_translation (tissue name) and
  tau_tissue_translation, calculated independently from that frame's normalized
  linear counts: tau = sum(1 - x_tissue / max(x)) / (number_of_tissues - 1).
  Tau ranges from 0 (uniform) to 1 (one tissue). All-zero rows have NA for both
  columns; tau is also NA with only one tissue. Tied maxima are joined with ';'.
  Raw count tables retain their original columns.
  Periodicity = p0 / (p0+p1+p2), or NA when no P-sites are assigned.
  Coverage = distinct p0 positions with positive counts / total p0 positions.
  Positive fractional bedGraph values count as observed; repeated positions
  across intervals/files count once. Coverage is zero when no p0 sites occur.
  The denominator is CDS nucleotide length / 3 for complete ORFs, or ceil(length/3)
  for incomplete codons, counting all actual p0 positions. Introns are excluded.
  Quality tables keep the same (orf_id, transcript_id) rows as the count tables.

Density PDFs:
  Periodicity/coverage have tissue columns and orf_biotype rows. Each panel
  reports its defined-value count.
  Binned Gaussian densities are reflected at the [0,1] boundaries; constant
  groups are drawn as vertical lines. NA values are omitted, zeros retained.
  Tau uses normalized p0 counts across tissues, with one panel per biotype.
  Each ORF contributes its single tau score once, without tissue subdivision.
  The barplot shows ORF counts by tissue (x axis) and biotype (panel rows):
  tissue-enriched = tau > 0.8; tissue-specific = tau > 0.95 (strict thresholds).
  These are nested groups: specific ORFs are also included in enriched counts.
  ORFs are assigned to their maximum normalized p0 tissue; tied maxima count
  in each tied tissue. Undefined tau values are excluded. Counts use the same
  (orf_id, transcript_id) rows as the tables.
  Install plotting dependencies with: python -m pip install matplotlib
  Use --skip-plots to generate only the TSV tables without plotting dependencies.

Progress report (stderr):
  After all files for each tissue are counted, report p0/p1/p2 percentages from
  summed raw counts over rows with orf_biotype exactly CDS, combining strands.
  Each percentage is 100 * CDS frame total / (CDS p0 + CDS p1 + CDS p2).
  These are pooled count percentages, not averages of per-ORF percentages.
  Raw frame totals are included; percentages are NA when the CDS total is zero.

Outputs (existing files with these names are overwritten):
  PREFIX.p0.tsv, PREFIX.p1.tsv, PREFIX.p2.tsv
  PREFIX.p0.normalized.tsv, PREFIX.p1.normalized.tsv, PREFIX.p2.normalized.tsv
  PREFIX.periodicity_coverage.tsv (metadata, TISSUE_periodicity, TISSUE_coverage, ...)
  PREFIX.periodicity.density.pdf, PREFIX.coverage.density.pdf
  PREFIX.tau_tissue_translation.density.pdf, PREFIX.tissue_enrichment.barplot.pdf
''')
    parser.add_argument('--gtf', required=True, type=Path,
                        help='ORF GTF: CDS features with orf_id/transcript_id attributes')
    parser.add_argument('--reads', choices=('ribo', 'rna'), default='ribo',
                        help='input read type (default: ribo)')
    parser.add_argument('--bedgraphs', '--input', dest='bedgraphs', required=True, type=Path,
                        help='input list: path tissue strand for ribo; BAM_path tissue for rna')
    parser.add_argument('--rna-strandedness', choices=('unstranded', 'forward', 'reverse'),
                        default='unstranded', help='RNA library orientation (default: unstranded)')
    parser.add_argument('--out-prefix', required=True,
                        help='output path prefix, e.g. results/nicos_psites')
    parser.add_argument('--target-total', type=float, default=1_000_000,
                        help='normalized total across all ORFs and frames per tissue (default: 1000000)')
    parser.add_argument('--skip-plots', action='store_true',
                        help='write all TSVs without requiring plotting dependencies')
    args = parser.parse_args()
    if not math.isfinite(args.target_total) or args.target_total <= 0:
        parser.error('--target-total must be finite and greater than zero')
    try:
        if not args.skip_plots:
            load_plotting()  # Fail before expensive counting if dependencies are missing.
        pysam = load_bam_reader() if args.reads == 'rna' else None
        samples = (read_bam_list(args.bedgraphs) if args.reads == 'rna'
                   else read_bedgraph_list(args.bedgraphs))
        rows, index, lengths = read_gtf(args.gtf, args.reads)
        print(f'Loaded {len(rows)} ORF/transcript rows and {len(samples)} tissues.', file=sys.stderr)
        counts, coverage = {}, {}
        p0_positions = [(length + 2) // 3 for length in lengths]
        for tissue, entries in samples.items():
            if args.reads == 'rna':
                counts[tissue] = array('d', [0.0]) * len(rows)
                for path in entries:
                    print(f'Counting {tissue}: {path}', file=sys.stderr)
                    count_bam(path, index, counts[tissue], pysam, args.rna_strandedness)
                continue
            counts[tissue] = array('d', [0.0]) * (3 * len(rows))
            covered = [bytearray(size) for size in p0_positions]
            for path, strand in entries:
                print(f'Counting {tissue} ({strand}): {path}', file=sys.stderr)
                count_bedgraph(path, strand, index, counts[tissue], covered)
            coverage[tissue] = array('d', (bitmap.count(1) / size
                                           for bitmap, size in zip(covered, p0_positions)))
            del covered
            report_cds_frames(tissue, rows, counts[tissue])
        specificity = write_tables(args.out_prefix, rows, samples, counts, args.target_total, args.reads)
        periodicity = (write_quality_table(args.out_prefix, rows, samples, counts, coverage)
                       if args.reads == 'ribo' else None)
        if not args.skip_plots:
            write_density_plots(args.out_prefix, rows, samples, periodicity, coverage,
                                specificity, args.reads)
            write_enrichment_barplot(args.out_prefix, rows, samples, specificity, args.reads)
    except (OSError, ValueError, OverflowError) as error:
        parser.exit(1, f'Error: {error}\n')


if __name__ == '__main__':
    main()
