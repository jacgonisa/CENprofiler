#!/usr/bin/env python3
"""
Large-scale analysis of indels at monomer-level resolution.

Processes multiple reads with large SVs and generates summary statistics.
"""

import pandas as pd
import subprocess
from pathlib import Path
import argparse

def main():
    parser = argparse.ArgumentParser(description='Large-scale monomer-level indel analysis')
    parser.add_argument('--bam', required=True, help='BAM file')
    parser.add_argument('--ref-fasta', required=True, help='Reference genome FASTA')
    parser.add_argument('--ref-monomers', required=True, help='Representative monomers FASTA')
    parser.add_argument('--cluster-file', required=True, help='Monomer family cluster file')
    parser.add_argument('--sv-reads-fa', required=True, help='FASTA file with SV reads')
    parser.add_argument('--sv-info', required=True, help='TSV with SV positions')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--max-reads', type=int, default=25, help='Max reads to process')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read SV info to get read IDs
    sv_info = pd.read_csv(args.sv_info, sep='\t')

    # Get top reads by number of large SVs
    sv_counts = sv_info.groupby('read_id').size().reset_index(name='sv_count')
    sv_counts = sv_counts.sort_values('sv_count', ascending=False)
    top_reads = sv_counts.head(args.max_reads)['read_id'].tolist()

    print(f"\n{'='*80}")
    print(f"Large-Scale Monomer-Level Indel Analysis")
    print(f"{'='*80}")
    print(f"Processing {len(top_reads)} reads with most large SVs (≥100bp)")
    print(f"Output directory: {output_dir}\n")

    # Step 1: Run FASTAN on all reads
    print("[1/5] Running FASTAN on reads...")
    fastan_bed = output_dir / 'all_reads.bed'
    if not fastan_bed.exists():
        subprocess.run([
            '/home/jg2070/bin/FasTAN',
            args.sv_reads_fa,
            str(output_dir / 'all_reads.1aln')
        ], check=True, capture_output=True)

        with open(fastan_bed, 'w') as out:
            subprocess.run([
                '/home/jg2070/alntools/tanbed',
                str(output_dir / 'all_reads.1aln')
            ], stdout=out, check=True)
        print(f"  ✓ FASTAN complete: {fastan_bed}")
    else:
        print(f"  ✓ Using existing: {fastan_bed}")

    # Step 2: Classify read monomers
    print("\n[2/5] Classifying read monomers...")
    read_monomers_tsv = output_dir / 'all_read_monomers.tsv'
    if not read_monomers_tsv.exists():
        subprocess.run([
            'python', 'CEN178profiler/scripts/classify_fastan_monomers_v2.py',
            '--fasta', args.sv_reads_fa,
            '--fastan-bed', str(fastan_bed),
            '--ref-monomers', args.ref_monomers,
            '--cluster-file', args.cluster_file,
            '--output', str(read_monomers_tsv)
        ], check=True)
        print(f"  ✓ Classification complete: {read_monomers_tsv}")
    else:
        print(f"  ✓ Using existing: {read_monomers_tsv}")

    # Step 3: Extract and classify deletion monomers for each read
    print("\n[3/5] Analyzing deletion monomers...")
    deletion_results = []
    for i, read_id in enumerate(top_reads[:10], 1):  # Do first 10 for deletions (slow)
        print(f"  [{i}/10] {read_id[:16]}...")
        del_output = output_dir / f'deletion_monomers_{read_id[:8]}.tsv'

        if not del_output.exists():
            try:
                subprocess.run([
                    'python', 'CEN178profiler/scripts/analyze_deletion_monomers.py',
                    '--bam', args.bam,
                    '--ref-fasta', args.ref_fasta,
                    '--read-id', read_id,
                    '--ref-monomers', args.ref_monomers,
                    '--cluster-file', args.cluster_file,
                    '--output', str(del_output)
                ], check=True, capture_output=True)
                deletion_results.append(read_id)
            except:
                print(f"    ⚠ Failed to analyze deletions for {read_id[:16]}")
        else:
            deletion_results.append(read_id)

    print(f"  ✓ Analyzed deletions for {len(deletion_results)} reads")

    # Step 4: Compute summary statistics
    print("\n[4/5] Computing summary statistics...")
    stats_script = output_dir / 'compute_stats.py'
    stats_script.write_text('''
import pandas as pd
from pathlib import Path
import sys

output_dir = Path(sys.argv[1])
sv_info_file = sys.argv[2]

# Load all data
read_monomers = pd.read_csv(output_dir / 'all_read_monomers.tsv', sep='\\t')
sv_info = pd.read_csv(sv_info_file, sep='\\t')

# Load deletion monomers
deletion_files = list(output_dir.glob('deletion_monomers_*.tsv'))
if deletion_files:
    deletion_dfs = [pd.read_csv(f, sep='\\t') for f in deletion_files]
    deletion_monomers = pd.concat(deletion_dfs, ignore_index=True)
else:
    deletion_monomers = pd.DataFrame()

print("\\n" + "="*80)
print("SUMMARY STATISTICS - Monomer-Level Indel Analysis")
print("="*80)

print(f"\\nDataset Overview:")
print(f"  Total reads analyzed: {read_monomers['read_id'].nunique()}")
print(f"  Total SVs (≥100bp): {len(sv_info)}")
print(f"  - Insertions: {len(sv_info[sv_info['sv_type'] == 'INS'])}")
print(f"  - Deletions: {len(sv_info[sv_info['sv_type'] == 'DEL'])}")

print(f"\\nRead Monomers (from insertions):")
print(f"  Total monomers detected: {len(read_monomers)}")
print(f"  Classified: {len(read_monomers[~read_monomers['monomer_family'].isna()])} ({100*len(read_monomers[~read_monomers['monomer_family'].isna()])/len(read_monomers):.1f}%)")
print(f"  Unclassified: {len(read_monomers[read_monomers['monomer_family'].isna()])}")

if len(deletion_monomers) > 0:
    print(f"\\nDeletion Monomers (from reference):")
    print(f"  Total monomers detected: {len(deletion_monomers)}")
    print(f"  Classified: {len(deletion_monomers[~deletion_monomers['monomer_family'].isna()])} ({100*len(deletion_monomers[~deletion_monomers['monomer_family'].isna()])/len(deletion_monomers):.1f}%)")
    print(f"  Unclassified: {len(deletion_monomers[deletion_monomers['monomer_family'].isna()])}")

print(f"\\nFamily Distribution (Read Monomers):")
classified_reads = read_monomers[~read_monomers['monomer_family'].isna()]
if len(classified_reads) > 0:
    family_counts = classified_reads['monomer_family'].value_counts().sort_index()
    for fam, count in family_counts.items():
        pct = 100 * count / len(classified_reads)
        print(f"  Family {int(fam):2d}: {count:5d} monomers ({pct:5.1f}%)")

if len(deletion_monomers) > 0:
    print(f"\\nFamily Distribution (Deletion Monomers):")
    classified_dels = deletion_monomers[~deletion_monomers['monomer_family'].isna()]
    if len(classified_dels) > 0:
        family_counts_del = classified_dels['monomer_family'].value_counts().sort_index()
        for fam, count in family_counts_del.items():
            pct = 100 * count / len(classified_dels)
            print(f"  Family {int(fam):2d}: {count:5d} monomers ({pct:5.1f}%)")

# Save summary
with open(output_dir / 'SUMMARY_STATS.txt', 'w') as f:
    f.write("MONOMER-LEVEL INDEL ANALYSIS - SUMMARY STATISTICS\\n")
    f.write("="*80 + "\\n\\n")
    f.write(f"Total reads: {read_monomers['read_id'].nunique()}\\n")
    f.write(f"Total large SVs: {len(sv_info)}\\n")
    f.write(f"Read monomers: {len(read_monomers)} ({len(read_monomers[~read_monomers['monomer_family'].isna()])} classified)\\n")
    if len(deletion_monomers) > 0:
        f.write(f"Deletion monomers: {len(deletion_monomers)} ({len(deletion_monomers[~deletion_monomers['monomer_family'].isna()])} classified)\\n")

print(f"\\n✓ Summary saved to {output_dir / 'SUMMARY_STATS.txt'}")
''')

    subprocess.run([
        'python', str(stats_script),
        str(output_dir),
        args.sv_info
    ], check=True)

    # Step 5: Generate visualizations for top 5 examples
    print("\n[5/5] Generating visualizations for top examples...")
    for i, read_id in enumerate(top_reads[:5], 1):
        print(f"  [{i}/5] {read_id[:16]}...")
        subprocess.run([
            'python', 'CEN178profiler/scripts/visualize_indel_families_v2.py',
            '--bam', args.bam,
            '--sv-info', args.sv_info,
            '--monomers', str(read_monomers_tsv),
            '--read-ids', read_id,
            '--output-dir', str(output_dir)
        ], check=True, capture_output=True)

    print(f"\n{'='*80}")
    print(f"✅ Analysis complete! Results in: {output_dir}")
    print(f"{'='*80}\n")

if __name__ == '__main__':
    main()
