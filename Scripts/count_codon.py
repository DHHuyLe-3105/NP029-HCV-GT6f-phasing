#!/usr/bin/env python3
"""
Read-level codon (haplotype) counting from BAM alignments.

For a given codon, tally the three-nucleotide haplotype carried by every read
that spans all three positions, requiring a minimum base quality at each
position. This preserves intra-codon phase and avoids the false amino-acid
calls produced when per-position frequencies are multiplied under an
independence assumption.

Reference: DQ835760.1 (HCV genotype 6f).
Loci used in this study (1-based):
    NS5B-282  nt 8463-8465
    NS5A-93   nt 6531-6533
    NS5A-28   nt 6336-6338

Usage:
    python count_codon.py --bam sample.bam --codon 8463 8464 8465 \
        --name NS5B-282 --timepoint Baseline
"""

import argparse
import csv
import sys
from collections import Counter

import pysam

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'AGT': 'S', 'AGC': 'S', 'CCT': 'P', 'CCC': 'P',
    'CCA': 'P', 'CCG': 'P', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CGT': 'R',
    'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R',
    'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': '*', 'TAG': '*', 'TGA': '*',
}


def translate(codon):
    return CODON_TABLE.get(codon, 'X')


def count_codon(bam_path, contig, positions, min_bq=20):
    """Return a Counter of codon -> read count for reads spanning `positions`.

    positions is a tuple of three 1-based reference coordinates (the codon).
    A read contributes only if it has a called, non-indel base at all three
    positions and each of those bases has quality >= min_bq.
    """
    pos0 = [p - 1 for p in positions]          # pysam is 0-based
    lo, hi = min(pos0), max(pos0)
    counts = Counter()

    bam = pysam.AlignmentFile(bam_path, 'rb')
    for read in bam.fetch(contig, lo, hi + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        seq = read.query_sequence
        qual = read.query_qualities
        if seq is None or qual is None:
            continue

        # map reference position -> position in the read (skips indels)
        ref_to_query = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
        if not all(p in ref_to_query for p in pos0):
            continue

        bases = []
        ok = True
        for p in pos0:
            q = ref_to_query[p]
            if qual[q] < min_bq:
                ok = False
                break
            bases.append(seq[q].upper())
        if ok:
            counts[''.join(bases)] += 1

    bam.close()
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bam', required=True, help='indexed BAM file')
    ap.add_argument('--contig', default='DQ835760.1', help='reference name in the BAM')
    ap.add_argument('--codon', required=True, type=int, nargs=3, metavar=('P1', 'P2', 'P3'),
                    help='three 1-based reference positions of the codon')
    ap.add_argument('--min-bq', type=int, default=20, help='minimum base quality (default 20)')
    ap.add_argument('--name', default='', help='locus label for the output')
    ap.add_argument('--timepoint', default='', help='sample/time-point label for the output')
    ap.add_argument('--out', help='write CSV here (default: print to stdout)')
    args = ap.parse_args()

    counts = count_codon(args.bam, args.contig, tuple(args.codon), args.min_bq)
    total = sum(counts.values())
    if total == 0:
        sys.exit('no reads span {} at {}'.format(args.codon, args.bam))

    rows = []
    for codon, n in counts.most_common():
        rows.append([args.timepoint, codon, translate(codon), n, total,
                     round(100 * n / total, 3)])

    header = ['timepoint', 'codon', 'amino_acid', 'count', 'total_spanning_reads', 'percent']
    out = open(args.out, 'w', newline='') if args.out else sys.stdout
    w = csv.writer(out)
    if args.name:
        w.writerow(['# locus', args.name, 'positions', '-'.join(map(str, args.codon)),
                    'min_bq', args.min_bq])
    w.writerow(header)
    w.writerows(rows)
    if args.out:
        out.close()
        sys.stderr.write('{}: {} codons, {} spanning reads\n'.format(
            args.name or args.bam, len(rows), total))


if __name__ == '__main__':
    main()
