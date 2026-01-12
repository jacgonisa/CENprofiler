#!/usr/bin/env python3
"""
Extract and organize monomer sequences by family.

This script:
- Reads monomer classifications and sequences
- Organizes sequences by family
- Calculates consensus sequences per family
- Performs within-family diversity analysis
- Generates sequence logos (requires logomaker)

Usage:
    python extract_monomer_sequences.py <classifications.tsv> <monomers.fa> <output_dir>
"""

import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import subprocess
import numpy as np

def load_classifications(tsv_file):
    """Load monomer classifications."""
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"Loaded {len(df)} monomer classifications")
    return df

def load_sequences(fasta_file):
    """Load monomer sequences into dictionary."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        sequences[record.id] = str(record.seq)
    print(f"Loaded {len(sequences)} monomer sequences")
    return sequences

def organize_by_family(df, sequences):
    """Organize sequences by family."""
    family_seqs = defaultdict(list)

    classified = df[df['monomer_family'].notna()]

    for _, row in classified.iterrows():
        monomer_id = row['monomer_id']
        family = int(row['monomer_family'])

        if monomer_id in sequences:
            family_seqs[family].append({
                'id': monomer_id,
                'seq': sequences[monomer_id],
                'seq_id': row['seq_id'],
                'array_idx': row['array_idx'],
                'identity': row['alignment_identity'],
                'length': len(sequences[monomer_id])
            })

    print(f"\nOrganized sequences into {len(family_seqs)} families:")
    for family in sorted(family_seqs.keys()):
        print(f"  Family {family}: {len(family_seqs[family])} monomers")

    return family_seqs

def save_family_fastas(family_seqs, output_dir):
    """Save FASTA files per family."""
    fasta_dir = output_dir / 'family_fastas'
    fasta_dir.mkdir(parents=True, exist_ok=True)

    for family, seqs in sorted(family_seqs.items()):
        fasta_file = fasta_dir / f'family_{family}.fa'

        records = []
        for i, seq_info in enumerate(seqs):
            record = SeqRecord(
                Seq(seq_info['seq']),
                id=seq_info['id'],
                description=f"family={family} array={seq_info['array_idx']} identity={seq_info['identity']:.1f}"
            )
            records.append(record)

        SeqIO.write(records, fasta_file, 'fasta')

    print(f"\nSaved {len(family_seqs)} family FASTA files to {fasta_dir}")

def calculate_sequence_diversity(family_seqs, output_dir):
    """Calculate within-family sequence diversity."""
    diversity_stats = []

    for family, seqs in sorted(family_seqs.items()):
        if len(seqs) < 2:
            continue

        # Calculate pairwise identities (simple approach)
        seq_strs = [s['seq'] for s in seqs[:100]]  # Limit to 100 for speed

        # Simple identity calculation
        identities = []
        for i in range(len(seq_strs)):
            for j in range(i+1, len(seq_strs)):
                s1, s2 = seq_strs[i], seq_strs[j]
                # Align at same length
                min_len = min(len(s1), len(s2))
                matches = sum(c1 == c2 for c1, c2 in zip(s1[:min_len], s2[:min_len]))
                identity = matches / min_len * 100
                identities.append(identity)

        stats = {
            'family': family,
            'n_sequences': len(seqs),
            'mean_pairwise_identity': np.mean(identities) if identities else 0,
            'std_pairwise_identity': np.std(identities) if identities else 0,
            'min_pairwise_identity': np.min(identities) if identities else 0,
            'max_pairwise_identity': np.max(identities) if identities else 0,
            'mean_length': np.mean([s['length'] for s in seqs]),
            'std_length': np.std([s['length'] for s in seqs])
        }
        diversity_stats.append(stats)

    df_diversity = pd.DataFrame(diversity_stats)

    # Save
    diversity_file = output_dir / 'family_diversity.tsv'
    df_diversity.to_csv(diversity_file, sep='\t', index=False)
    print(f"\nSaved diversity statistics to {diversity_file}")

    # Save report
    report_file = output_dir / 'family_diversity_report.txt'
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("WITHIN-FAMILY SEQUENCE DIVERSITY\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"{'Family':<8} {'N':>6} {'Mean ID':>10} {'Std ID':>10} {'Min ID':>10} {'Max ID':>10}\n")
        f.write("-" * 80 + "\n")

        for _, row in df_diversity.iterrows():
            f.write(f"F{int(row['family']):<7} {int(row['n_sequences']):>6} "
                   f"{row['mean_pairwise_identity']:>10.1f} {row['std_pairwise_identity']:>10.1f} "
                   f"{row['min_pairwise_identity']:>10.1f} {row['max_pairwise_identity']:>10.1f}\n")

        f.write("\n")
        f.write("Note: Diversity calculated from pairwise sequence identity (up to 100 seqs per family)\n")
        f.write("=" * 80 + "\n")

    print(f"Saved diversity report to {report_file}")

    return df_diversity

def calculate_consensus(family_seqs, output_dir, min_seqs=10):
    """Calculate consensus sequences for families with sufficient data."""
    consensus_dir = output_dir / 'consensus'
    consensus_dir.mkdir(parents=True, exist_ok=True)

    consensus_records = []

    for family, seqs in sorted(family_seqs.items()):
        if len(seqs) < min_seqs:
            print(f"  Family {family}: Too few sequences ({len(seqs)} < {min_seqs}), skipping consensus")
            continue

        print(f"  Family {family}: Calculating consensus from {len(seqs)} sequences...")

        # Take up to 50 sequences for alignment (for speed)
        selected_seqs = seqs[:50]

        # Create temporary FASTA
        temp_fa = consensus_dir / f'temp_family_{family}.fa'
        records = [SeqRecord(Seq(s['seq']), id=f"seq{i}", description="")
                  for i, s in enumerate(selected_seqs)]
        SeqIO.write(records, temp_fa, 'fasta')

        # Try to run MUSCLE alignment
        try:
            aligned_fa = consensus_dir / f'temp_family_{family}_aligned.fa'

            # Run muscle
            result = subprocess.run(['muscle', '-align', str(temp_fa),
                                   '-output', str(aligned_fa)],
                                  capture_output=True, text=True, timeout=60)

            if result.returncode != 0:
                print(f"    MUSCLE failed for family {family}, skipping")
                temp_fa.unlink()
                continue

            # Read alignment and calculate consensus
            alignment = AlignIO.read(aligned_fa, 'fasta')

            # Simple consensus: most common base at each position
            consensus_seq = []
            for i in range(alignment.get_alignment_length()):
                column = alignment[:, i]
                # Count non-gap bases
                bases = [b for b in column if b != '-']
                if bases:
                    most_common = max(set(bases), key=bases.count)
                    consensus_seq.append(most_common)

            consensus_str = ''.join(consensus_seq)

            # Save consensus
            consensus_record = SeqRecord(
                Seq(consensus_str),
                id=f"Family_{family}_consensus",
                description=f"n={len(selected_seqs)} length={len(consensus_str)}bp"
            )
            consensus_records.append(consensus_record)

            # Cleanup
            temp_fa.unlink()
            aligned_fa.unlink()

        except subprocess.TimeoutExpired:
            print(f"    MUSCLE timeout for family {family}, skipping")
            if temp_fa.exists():
                temp_fa.unlink()
        except FileNotFoundError:
            print(f"    MUSCLE not found, skipping consensus calculation")
            break
        except Exception as e:
            print(f"    Error processing family {family}: {e}")
            if temp_fa.exists():
                temp_fa.unlink()

    # Save all consensus sequences
    if consensus_records:
        consensus_file = consensus_dir / 'all_consensus.fa'
        SeqIO.write(consensus_records, consensus_file, 'fasta')
        print(f"\nSaved {len(consensus_records)} consensus sequences to {consensus_file}")
    else:
        print("\nNo consensus sequences generated (MUSCLE may not be installed)")

def generate_sequence_summary(family_seqs, df_diversity, output_dir):
    """Generate comprehensive sequence summary."""
    summary_file = output_dir / 'sequence_summary.txt'

    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MONOMER SEQUENCE ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Total families with sequences: {len(family_seqs)}\n")
        f.write(f"Total monomer sequences: {sum(len(seqs) for seqs in family_seqs.values())}\n\n")

        f.write("PER-FAMILY SUMMARY\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Family':<8} {'Count':>8} {'Mean Len':>10} {'Std Len':>10} {'Diversity':>12}\n")
        f.write("-" * 80 + "\n")

        for family in sorted(family_seqs.keys()):
            seqs = family_seqs[family]
            mean_len = np.mean([s['length'] for s in seqs])
            std_len = np.std([s['length'] for s in seqs])

            # Get diversity if available
            div_row = df_diversity[df_diversity['family'] == family]
            if len(div_row) > 0:
                diversity = f"{div_row.iloc[0]['mean_pairwise_identity']:.1f}%"
            else:
                diversity = "N/A"

            f.write(f"F{family:<7} {len(seqs):>8} {mean_len:>10.1f} {std_len:>10.1f} {diversity:>12}\n")

        f.write("\n")
        f.write("=" * 80 + "\n")

    print(f"\nSaved sequence summary to {summary_file}")

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_monomer_sequences.py <classifications.tsv> <monomers.fa> <output_dir>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_dir = Path(sys.argv[3])
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*80)
    print("MONOMER SEQUENCE EXTRACTION AND ANALYSIS")
    print("="*80 + "\n")

    # Load data
    print("Loading classifications...")
    df = load_classifications(tsv_file)

    print("Loading sequences...")
    sequences = load_sequences(fasta_file)

    # Organize by family
    print("\nOrganizing by family...")
    family_seqs = organize_by_family(df, sequences)

    # Save family FASTAs
    print("\nSaving family FASTA files...")
    save_family_fastas(family_seqs, output_dir)

    # Calculate diversity
    print("\nCalculating within-family diversity...")
    df_diversity = calculate_sequence_diversity(family_seqs, output_dir)

    # Calculate consensus (if MUSCLE available)
    print("\nCalculating consensus sequences...")
    calculate_consensus(family_seqs, output_dir, min_seqs=10)

    # Generate summary
    print("\nGenerating summary...")
    generate_sequence_summary(family_seqs, df_diversity, output_dir)

    print("\n" + "="*80)
    print("SEQUENCE EXTRACTION COMPLETE")
    print(f"Results saved to: {output_dir}")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
