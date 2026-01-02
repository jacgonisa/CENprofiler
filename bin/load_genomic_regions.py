#!/usr/bin/env python3
"""
Load genomic region annotations (centromeres, pericentromeres, rDNA).
Handles chromosome name mapping between different formats (Chr1 vs CP accessions).

Outputs a unified regions file for downstream analysis.
"""

import sys
import os
import argparse
from pathlib import Path

# Chromosome name mapping: CP accession → Chr format
CHROM_NAME_MAP = {
    'CP116280.1': 'Chr1',
    'CP116281.2': 'Chr2',
    'CP116282.1': 'Chr3',
    'CP116283.2': 'Chr4',
    'CP116284.1': 'Chr5',
}


def load_bed_file(bed_path, region_type):
    """
    Load regions from BED file.
    Returns list of (chrom, start, end, type) tuples.
    Handles chromosome name mapping.
    """
    regions = []

    if not os.path.exists(bed_path):
        print(f"Warning: BED file not found: {bed_path}", file=sys.stderr)
        return regions

    with open(bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            # Map chromosome names if needed
            if chrom in CHROM_NAME_MAP:
                chrom = CHROM_NAME_MAP[chrom]

            regions.append((chrom, start, end, region_type))

    return regions


def load_all_regions(annotation_dir):
    """
    Load all genomic regions from annotation directory.

    Expected files:
    - centromere.bed
    - pericentromere_clean.bed (or pericentromere_raw.bed)
    - 5s_rdna_regions.bed
    - 45s_rdna_regions.bed

    Returns: list of (chrom, start, end, region_type) tuples
    """
    anno_path = Path(annotation_dir)
    all_regions = []

    # Load centromeres
    cen_file = anno_path / 'centromere.bed'
    if cen_file.exists():
        regions = load_bed_file(cen_file, 'centromere')
        all_regions.extend(regions)
        print(f"Loaded {len(regions)} centromere regions", file=sys.stderr)

    # Load pericentromeres (prefer clean version)
    pericen_file = anno_path / 'pericentromere_clean.bed'
    if not pericen_file.exists():
        pericen_file = anno_path / 'pericentromere_raw.bed'

    if pericen_file.exists():
        regions = load_bed_file(pericen_file, 'pericentromere')
        all_regions.extend(regions)
        print(f"Loaded {len(regions)} pericentromere regions", file=sys.stderr)

    # Load 5S rDNA
    rdna_5s_file = anno_path / '5s_rdna_regions.bed'
    if rdna_5s_file.exists():
        regions = load_bed_file(rdna_5s_file, '5s_rdna')
        all_regions.extend(regions)
        print(f"Loaded {len(regions)} 5S rDNA regions", file=sys.stderr)

    # Load 45S rDNA
    rdna_45s_file = anno_path / '45s_rdna_regions.bed'
    if rdna_45s_file.exists():
        regions = load_bed_file(rdna_45s_file, '45s_rdna')
        all_regions.extend(regions)
        print(f"Loaded {len(regions)} 45S rDNA regions", file=sys.stderr)

    return all_regions


def write_unified_regions(regions, output_file):
    """
    Write all regions to a unified TSV file.
    Format: chrom  start  end  region_type
    """
    with open(output_file, 'w') as f:
        f.write("chrom\tstart\tend\tregion_type\n")
        for chrom, start, end, region_type in sorted(regions):
            f.write(f"{chrom}\t{start}\t{end}\t{region_type}\n")

    print(f"Wrote {len(regions)} regions to {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Load genomic region annotations from BED files'
    )
    parser.add_argument('annotation_dir',
                       help='Directory containing BED annotation files')
    parser.add_argument('-o', '--output', default='genomic_regions.tsv',
                       help='Output unified regions file (default: genomic_regions.tsv)')

    args = parser.parse_args()

    # Check annotation directory exists
    if not os.path.isdir(args.annotation_dir):
        print(f"Error: Annotation directory not found: {args.annotation_dir}",
              file=sys.stderr)
        sys.exit(1)

    # Load all regions
    print(f"Loading annotations from: {args.annotation_dir}", file=sys.stderr)
    regions = load_all_regions(args.annotation_dir)

    if not regions:
        print("Warning: No regions loaded!", file=sys.stderr)
        sys.exit(1)

    # Write unified output
    write_unified_regions(regions, args.output)

    print("✅ Genomic regions loaded successfully", file=sys.stderr)


if __name__ == '__main__':
    main()
