# DNA Sequence Tools

Python utilities for DNA sequence manipulation:

- **Reverse complement** of a DNA sequence
- **Calculate DNA weight** (ng) from moles (fmol) and size (kb)
- **Find non-cutting enzymes** in SnapGene `.dna` files
- **Convert `.dna` files** to FASTA and report sizes

## Installation

```bash
pip install snapgene-reader
```

## Usage
```bash
python dna_tools.py [options]
```

## Options
| Flag    | Description                                            | Example                  |
|---------|--------------------------------------------------------|--------------------------|
| `-r`    | Reverse complement sequence                            | `-r ATCGGA`              |
| `-m2w`  | Calculate weight: moles (fmol) and size (kb)          | `-m2w 50 2`              |
| `-enc`  | Find enzymes that don’t cut sequence in `.dna` file   | `-enc plasmid.dna`       |
| `-fc`   | Convert `.dna` files to FASTA and print sizes         | `-fc plasmid.dna`       |
