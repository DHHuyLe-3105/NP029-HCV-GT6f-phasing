# HCV GT6f NP029 — read-level haplotype phasing

Data and code supporting the case report on haplotype phasing at NS5B codon 282
in a hepatitis C virus genotype 6f treatment failure (patient NP029). The
alignments, Sanger traces, counting script and per-codon count tables here let
the read-level analyses in the paper (Table 1, Tables S4–S6, Figure 1) be
reproduced from the primary data.

All reads are aligned to the HCV genotype 6f reference **DQ835760.1**. Codon
positions used throughout (1-based, on DQ835760.1):

| Locus | Codon | Nucleotides | Substitution |
|-------|-------|-------------|--------------|
| NS5B  | 282   | 8463–8465   | S282C / S282T (RAS) |
| NS5A  | 93    | 6531–6533   | T93S |
| NS5A  | 28    | 6336–6338   | A28T |

NS5A and NS5B were amplified and sequenced separately, so the two genes cannot
be phased against each other.

## Layout

```
.
├── data/
│   ├── bam/            # 6 alignments (NS5A, NS5B × 3 time points) + .bai indexes
│   │   ├── NS5A/
│   │   └── NS5B/
│   ├── sanger/         # confirmatory Sanger traces (.ab1)
│   └── counts/         # per-codon, per-time-point haplotype counts (CSV)
└── scripts/
    └── count_codon.py  # read-level codon counting
```

## Alignments

Each library is one gene at one time point. Reads were aligned to DQ835760.1 and
the BAMs indexed. Spanning-read depth at the reported codon (base quality ≥ 20)
is given for reference.

| File | Gene | Time point | Spanning reads |
|------|------|-----------|----------------|
| `NS5B/L1_NS5B_Baseline.bam` | NS5B | Baseline | 3385 (codon 282) |
| `NS5B/L2_NS5B_Week24.bam`   | NS5B | Week 24  | 2297 (codon 282) |
| `NS5B/L3_NS5B_Week72.bam`   | NS5B | Week 72  | 2963 (codon 282) |
| `NS5A/L4_NS5A_Baseline.bam` | NS5A | Baseline | 3017 / 1974 (codon 93 / 28) |
| `NS5A/L5_NS5A_Week24.bam`   | NS5A | Week 24  | 2825 / 2888 (codon 93 / 28) |
| `NS5A/L6_NS5A_Week72.bam`   | NS5A | Week 72  | 1735 / 1697 (codon 93 / 28) |

## Requirements

Python 3 with [pysam](https://pysam.readthedocs.io):

```
pip install pysam
```

## Reproducing the counts

`scripts/count_codon.py` tallies the three-nucleotide haplotype carried by every
read spanning a codon, keeping only positions with base quality ≥ 20. Counting
whole codons (rather than multiplying per-position frequencies) preserves phase
and is what distinguishes the real serine sweep from the fabricated S282C/S282T
codons.

```
python scripts/count_codon.py --bam data/bam/NS5B/L1_NS5B_Baseline.bam \
    --codon 8463 8464 8465 --name NS5B-282 --timepoint Baseline
```

Run once per BAM and codon (positions above), concatenating the output to
rebuild the tables in `data/counts/`. `--min-bq` changes the quality cut-off
(default 20); `--out file.csv` writes to disk instead of stdout.

## Sanger confirmation

`data/sanger/` holds ABI trace files (`.ab1`) used to confirm the NGS calls at
codon 282 — at baseline the trace shows a clean AGT (serine) with no secondary
peak, i.e. no S282C/S282T.

## Data availability

- Consensus sequences (6): GenBank **PZ191478–PZ191483**
- This repository, archived with a permanent DOI: **[Zenodo DOI]**
- Raw sequencing reads (FASTQ): SRA **BioProject [PRJNA######]**

## Citation

If you use these data, please cite the paper:

> Le,D.H.H., Puenpa,J., Nilyanimit,P., Honsawek,S. and Poovorawan,Y. Center of Excellence in Clinical Virology,
            Chulalongkorn University, 1873, Rama 4 Road, Pathumwan, Bangkok
            10330, Thailand

## License

Code is released under the MIT License; data files under CC-BY-4.0. See `LICENSE`.
