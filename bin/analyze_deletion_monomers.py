#!/usr/bin/env python3
"""
Analyze monomers in deletion regions by extracting reference sequence.

For deletions, the sequence exists in the reference but not in the read.
This script extracts reference sequence from deletions and classifies monomers.
"""

import pandas as pd
import pysam
import argparse
import subprocess
import tempfile
from pathlib import Path

def extract_deletion_sequences(bam_file, ref_fasta, read_id, output_fasta):
    """Extract reference sequences for all large deletions in a read."""

    bam = pysam.AlignmentFile(bam_file, 'rb')

    # Find the read
    alignment = None
    for read in bam.fetch():
        if read.query_name == read_id:
            alignment = read
            break

    if alignment is None:
        print(f"Read {read_id} not found")
        bam.close()
        return []

    # Get reference name and parse CIGAR
    ref_name = alignment.reference_name
    ref_start = alignment.reference_start

    bam.close()

    # Open reference FASTA
    ref_fa = pysam.FastaFile(ref_fasta)

    deletions = []
    ref_pos = ref_start

    for op, length in alignment.cigartuples:
        if op == 0:  # M
            ref_pos += length
        elif op == 2:  # D - deletion
            if length >= 100:  # Only large deletions
                # Extract reference sequence
                del_seq = ref_fa.fetch(ref_name, ref_pos, ref_pos + length)

                deletions.append({
                    'read_id': read_id,
                    'ref_name': ref_name,
                    'ref_start': ref_pos,
                    'ref_end': ref_pos + length,
                    'length': length,
                    'sequence': del_seq
                })
            ref_pos += length

    ref_fa.close()

    # Write to FASTA
    if len(deletions) > 0:
        with open(output_fasta, 'w') as f:
            for i, d in enumerate(deletions):
                f.write(f">{read_id}_del{i}_{d['ref_name']}:{d['ref_start']}-{d['ref_end']}\n")
                f.write(f"{d['sequence']}\n")

    return deletions

def run_fastan_on_deletions(fasta_file):
    """Run FASTAN on deletion sequences."""

    if not Path(fasta_file).exists():
        return None

    # Run FASTAN
    output_1aln = str(fasta_file).replace('.fa', '.1aln')
    subprocess.run([
        '/home/jg2070/bin/FasTAN',
        str(fasta_file),
        output_1aln
    ], check=True, capture_output=True)

    # Convert to BED with tanbed
    output_bed = str(fasta_file).replace('.fa', '.bed')
    with open(output_bed, 'w') as out:
        subprocess.run([
            '/home/jg2070/alntools/tanbed',
            output_1aln
        ], stdout=out, check=True)

    return output_bed

def classify_deletion_monomers(bed_file, fasta_file, ref_monomers, cluster_file):
    """Classify monomers from deletion regions."""

    if not Path(bed_file).exists():
        return pd.DataFrame()

    # Read FASTAN BED output
    bed_df = pd.read_csv(bed_file, sep='\t', header=None,
                         names=['seq_id', 'start', 'end', 'period', 'quality'])

    # Filter for CEN178-like periods (160-200bp)
    cen178_arrays = bed_df[(bed_df['period'] >= 160) & (bed_df['period'] <= 200)].copy()

    if len(cen178_arrays) == 0:
        return pd.DataFrame()

    # Extract monomers
    fasta = pysam.FastaFile(fasta_file)

    all_monomers = []

    for _, array in cen178_arrays.iterrows():
        seq_id = array['seq_id']
        start = array['start']
        end = array['end']
        period = int(array['period'])

        # Get sequence
        array_seq = fasta.fetch(seq_id, start, end)

        # Extract individual monomers
        pos = 0
        mon_idx = 0
        while pos + period <= len(array_seq):
            monomer_seq = array_seq[pos:pos + period]

            monomer_id = f"{seq_id}_array{0}_mon{mon_idx}"
            all_monomers.append({
                'monomer_id': monomer_id,
                'seq_id': seq_id,
                'array_start': start,
                'monomer_start': start + pos,
                'monomer_end': start + pos + period,
                'period': period,
                'sequence': monomer_seq
            })

            pos += period
            mon_idx += 1

    fasta.close()

    if len(all_monomers) == 0:
        return pd.DataFrame()

    # Write monomers to FASTA for classification
    monomer_fasta = str(fasta_file).replace('.fa', '_monomers.fa')
    with open(monomer_fasta, 'w') as f:
        for mon in all_monomers:
            f.write(f">{mon['monomer_id']}\n{mon['sequence']}\n")

    # Map to representative monomers with minimap2 (same settings as classify_fastan_monomers_v2.py)
    paf_file = monomer_fasta.replace('.fa', '.paf')
    subprocess.run([
        'minimap2',
        '-x', 'asm20',
        '-t', '4',
        '--eqx',
        ref_monomers,
        monomer_fasta
    ], stdout=open(paf_file, 'w'), stderr=subprocess.DEVNULL, check=True)

    # Parse PAF and assign families
    cluster_map = {}
    with open(cluster_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                monomer = parts[0]
                family = parts[1]
                cluster_map[monomer] = family

    # Read PAF results
    paf_results = {}
    with open(paf_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            query = fields[0]
            target = fields[5]
            mapq = int(fields[11])

            # Calculate identity
            matches = int([x.split(':')[2] for x in fields[12:] if x.startswith('de:f:')][0] if any(x.startswith('de:f:') for x in fields[12:]) else 0)
            aln_len = int(fields[10])
            identity = (1 - matches / aln_len) * 100 if aln_len > 0 else 0

            if query not in paf_results or identity > paf_results[query][1]:
                paf_results[query] = (target, identity, mapq)

    # Assign families to monomers
    for mon in all_monomers:
        mon_id = mon['monomer_id']
        if mon_id in paf_results:
            target, identity, mapq = paf_results[mon_id]
            if identity >= 60:  # Minimum identity threshold (lowered from 70% for deletion monomers)
                mon['best_match'] = target
                mon['alignment_identity'] = identity
                mon['mapq'] = mapq
                mon['monomer_family'] = cluster_map.get(target, None)
            else:
                mon['best_match'] = None
                mon['alignment_identity'] = None
                mon['mapq'] = None
                mon['monomer_family'] = None
        else:
            mon['best_match'] = None
            mon['alignment_identity'] = None
            mon['mapq'] = None
            mon['monomer_family'] = None

    return pd.DataFrame(all_monomers)

def main():
    parser = argparse.ArgumentParser(description='Analyze monomers in deletion regions')
    parser.add_argument('--bam', required=True, help='BAM file')
    parser.add_argument('--ref-fasta', required=True, help='Reference genome FASTA')
    parser.add_argument('--read-id', required=True, help='Read ID')
    parser.add_argument('--ref-monomers', required=True, help='Representative monomers FASTA')
    parser.add_argument('--cluster-file', required=True, help='Monomer to family cluster file')
    parser.add_argument('--output', required=True, help='Output TSV file')

    args = parser.parse_args()

    # Create temp directory
    temp_dir = Path(args.output).parent / 'temp_deletions'
    temp_dir.mkdir(exist_ok=True)

    # Extract deletion sequences
    del_fasta = temp_dir / f'{args.read_id}_deletions.fa'
    deletions = extract_deletion_sequences(args.bam, args.ref_fasta, args.read_id, del_fasta)

    if len(deletions) == 0:
        print(f"No large deletions found for {args.read_id}")
        return

    print(f"Found {len(deletions)} large deletions")

    # Run FASTAN
    bed_file = run_fastan_on_deletions(del_fasta)

    # Classify monomers
    monomers_df = classify_deletion_monomers(
        bed_file, del_fasta,
        args.ref_monomers, args.cluster_file
    )

    if len(monomers_df) > 0:
        monomers_df.to_csv(args.output, sep='\t', index=False)
        print(f"Classified {len(monomers_df)} monomers from deletions")
        print(f"Saved to {args.output}")
    else:
        print("No CEN178 monomers found in deletions")

if __name__ == '__main__':
    main()
