#!/usr/bin/env python3
"""
Extract reads with large indels from BAM file.

Parses CIGAR strings to identify insertions and deletions ≥ min_size.
Classifies indels by genomic region (centromere, pericentromere, arms, rDNA).
Outputs:  - reads FASTA file
- indel catalog TSV
- extraction statistics

Usage:
    python extract_reads_from_bam.py <bam_file> <regions_file> <min_size> <output_prefix>
"""

import sys
import os
import pysam
import csv
from collections import defaultdict


def load_genomic_regions(regions_file):
    """
    Load genomic regions from TSV file.
    Returns dict: {chrom: [(start, end, region_type), ...]}
    """
    regions = defaultdict(list)

    with open(regions_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            region_type = fields[3]

            regions[chrom].append((start, end, region_type))

    # Sort regions by start position for efficient lookup
    for chrom in regions:
        regions[chrom].sort()

    return regions


def classify_position(chrom, pos, regions):
    """
    Classify genomic position by region.
    Priority: 5s_rdna > 45s_rdna > centromere > pericentromere > arms
    """
    if chrom not in regions:
        return 'other'

    # Check regions in priority order
    priority_order = ['5s_rdna', '45s_rdna', 'centromere', 'pericentromere']

    for region_type in priority_order:
        for start, end, rtype in regions[chrom]:
            if rtype == region_type and start <= pos < end:
                return region_type

    # Default to arms if not in any special region
    return 'arms'


def extract_reads_with_large_indels(bam_file, regions, min_indel_size, reads_output, catalog_output):
    """
    Extract reads containing large indels from BAM file.

    Returns: number of reads extracted, total indels found
    """
    try:
        bamfile = pysam.AlignmentFile(bam_file, "rb")
    except Exception as e:
        print(f"Error opening BAM file: {e}", file=sys.stderr)
        return 0, 0

    # Track extracted reads (avoid duplicates)
    extracted_reads = set()
    indel_catalog = []

    # Statistics
    total_reads_scanned = 0
    reads_with_large_indels = 0
    total_large_indels = 0

    # Open output files
    reads_fasta = open(reads_output, 'w')

    print("Scanning BAM for reads with large indels...", file=sys.stderr)

    # Process each chromosome
    for chrom in bamfile.references:
        print(f"  Processing {chrom}...", file=sys.stderr)
        chrom_reads = 0
        chrom_indels = 0

        for read in bamfile.fetch(chrom):
            # Skip unmapped, secondary, supplementary
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            total_reads_scanned += 1

            # Parse CIGAR for indels
            read_pos = 0
            ref_pos = read.reference_start
            has_large_indel = False

            if read.cigartuples:
                for op, length in read.cigartuples:
                    if op == 0:  # M (match/mismatch)
                        read_pos += length
                        ref_pos += length

                    elif op == 1:  # I (insertion)
                        if length >= min_indel_size:
                            has_large_indel = True
                            total_large_indels += 1
                            chrom_indels += 1

                            region = classify_position(chrom, ref_pos, regions)

                            # Calculate CEN178 multiple info
                            multiple_of_178 = (length % 178 == 0)
                            closest_178 = round(length / 178, 2)
                            closest_multiple = round(length / 178) * 178
                            distance_to_multiple = abs(length - closest_multiple)

                            indel_catalog.append({
                                'read_id': read.query_name,
                                'chromosome': chrom,
                                'ref_pos': ref_pos,
                                'read_pos': read_pos,
                                'type': 'Insertion',
                                'size': length,
                                'region': region,
                                'mapping_quality': read.mapping_quality,
                                'multiple_of_178': multiple_of_178,
                                'closest_178_multiple': closest_178,
                                'distance_to_178_multiple': distance_to_multiple
                            })

                        read_pos += length

                    elif op == 2:  # D (deletion)
                        if length >= min_indel_size:
                            has_large_indel = True
                            total_large_indels += 1
                            chrom_indels += 1

                            region = classify_position(chrom, ref_pos, regions)

                            # Calculate CEN178 multiple info
                            multiple_of_178 = (length % 178 == 0)
                            closest_178 = round(length / 178, 2)
                            closest_multiple = round(length / 178) * 178
                            distance_to_multiple = abs(length - closest_multiple)

                            indel_catalog.append({
                                'read_id': read.query_name,
                                'chromosome': chrom,
                                'ref_pos': ref_pos,
                                'read_pos': read_pos,
                                'type': 'Deletion',
                                'size': length,
                                'region': region,
                                'mapping_quality': read.mapping_quality,
                                'multiple_of_178': multiple_of_178,
                                'closest_178_multiple': closest_178,
                                'distance_to_178_multiple': distance_to_multiple
                            })

                        ref_pos += length

                    elif op == 4:  # S (soft clip)
                        read_pos += length

            # Extract read if it has large indels and not already extracted
            if has_large_indel and read.query_name not in extracted_reads:
                extracted_reads.add(read.query_name)
                reads_with_large_indels += 1
                chrom_reads += 1

                # Write to FASTA
                sequence = read.query_sequence
                if sequence:
                    reads_fasta.write(f">{read.query_name}\n{sequence}\n")

        print(f"    Found {chrom_reads} reads with {chrom_indels} large indels", file=sys.stderr)

    # Close files
    reads_fasta.close()
    bamfile.close()

    # Write indel catalog
    if indel_catalog:
        fieldnames = ['read_id', 'chromosome', 'ref_pos', 'read_pos', 'type',
                     'size', 'region', 'mapping_quality', 'multiple_of_178',
                     'closest_178_multiple', 'distance_to_178_multiple']

        with open(catalog_output, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(indel_catalog)

    # Print summary
    print(f"\n=== Extraction Summary ===", file=sys.stderr)
    print(f"Total reads scanned: {total_reads_scanned:,}", file=sys.stderr)
    print(f"Reads with large indels (≥{min_indel_size}bp): {reads_with_large_indels:,}", file=sys.stderr)
    print(f"Total large indels: {total_large_indels:,}", file=sys.stderr)
    print(f"Unique reads extracted: {len(extracted_reads):,}", file=sys.stderr)

    return len(extracted_reads), total_large_indels


def main():
    if len(sys.argv) != 5:
        print("Usage: python extract_reads_from_bam.py <bam_file> <regions_file> <min_size> <output_prefix>",
              file=sys.stderr)
        sys.exit(1)

    bam_file = sys.argv[1]
    regions_file = sys.argv[2]
    min_indel_size = int(sys.argv[3])
    output_prefix = sys.argv[4]

    reads_output = f"{output_prefix}_reads.fa"
    catalog_output = f"{output_prefix}_indel_catalog.tsv"

    # Load genomic regions
    print(f"Loading genomic regions from: {regions_file}", file=sys.stderr)
    regions = load_genomic_regions(regions_file)

    region_counts = sum(len(v) for v in regions.values())
    print(f"Loaded {region_counts} regions across {len(regions)} chromosomes", file=sys.stderr)

    # Extract reads
    print(f"\nExtracting reads with indels ≥{min_indel_size}bp from: {bam_file}", file=sys.stderr)
    n_reads, n_indels = extract_reads_with_large_indels(
        bam_file, regions, min_indel_size, reads_output, catalog_output
    )

    # Write stats file
    stats_output = f"{output_prefix}_stats.txt"
    with open(stats_output, 'w') as f:
        f.write(f"Reads extracted: {n_reads}\n")
        f.write(f"Large indels found: {n_indels}\n")
        f.write(f"Minimum indel size: {min_indel_size}bp\n")

    print(f"\n✅ Extraction complete!", file=sys.stderr)
    print(f"   Reads FASTA: {reads_output}", file=sys.stderr)
    print(f"   Indel catalog: {catalog_output}", file=sys.stderr)
    print(f"   Statistics: {stats_output}", file=sys.stderr)


if __name__ == '__main__':
    main()
