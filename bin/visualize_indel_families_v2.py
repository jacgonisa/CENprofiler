#!/usr/bin/env python3
"""
Simplified visualization showing all monomers colored by family on both tracks.

Highlights which CEN178 families are present in insertions and deletions,
while also showing all monomers (even those not in indels) colored by family.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle, Polygon
import numpy as np
import argparse
import pysam
from pathlib import Path

# Family color scheme
FAMILY_COLORS = {
    1: '#e41a1c',   # Red
    2: '#377eb8',   # Blue
    3: '#4daf4a',   # Green
    4: '#984ea3',   # Purple
    7: '#ff7f00',   # Orange
    8: '#ffff33',   # Yellow
    11: '#a65628',  # Brown
    13: '#f781bf',  # Pink
    18: '#00ffff',  # Cyan
}

def parse_cigar(cigar_tuples, ref_start):
    """Parse CIGAR string into alignment blocks."""
    blocks = []
    ref_pos = ref_start
    read_pos = 0

    for op, length in cigar_tuples:
        if op == 0:  # M - match/mismatch
            blocks.append({
                'type': 'M',
                'ref_start': ref_pos,
                'ref_end': ref_pos + length,
                'read_start': read_pos,
                'read_end': read_pos + length,
            })
            ref_pos += length
            read_pos += length
        elif op == 1:  # I - insertion
            blocks.append({
                'type': 'I',
                'ref_start': ref_pos,
                'ref_end': ref_pos,
                'read_start': read_pos,
                'read_end': read_pos + length,
            })
            read_pos += length
        elif op == 2:  # D - deletion
            blocks.append({
                'type': 'D',
                'ref_start': ref_pos,
                'ref_end': ref_pos + length,
                'read_start': read_pos,
                'read_end': read_pos,
            })
            ref_pos += length
        elif op == 4:  # S - soft clip
            read_pos += length

    return blocks

def plot_indel_families_with_all_monomers(bam_file, sv_info_file, monomers_file, read_id, output_file, deletion_monomers_file=None):
    """Visualization showing all monomers colored by family, with indels highlighted."""

    # Load data
    sv_info = pd.read_csv(sv_info_file, sep='\t')
    monomers = pd.read_csv(monomers_file, sep='\t')

    # Load deletion monomers if available
    deletion_monomers = None
    if deletion_monomers_file and Path(deletion_monomers_file).exists():
        deletion_monomers = pd.read_csv(deletion_monomers_file, sep='\t')

    # Filter for this read
    read_svs = sv_info[sv_info['read_id'] == read_id]
    read_monomers = monomers[monomers['read_id'] == read_id]

    if len(read_svs) == 0:
        print(f"No SVs found for read {read_id}")
        return

    # Get alignment from BAM
    bam = pysam.AlignmentFile(bam_file, 'rb')
    alignment = None
    for read in bam.fetch():
        if read.query_name == read_id:
            alignment = read
            break
    bam.close()

    if alignment is None or alignment.is_unmapped:
        print(f"Read {read_id} not found or unmapped")
        return

    # Parse CIGAR
    blocks = parse_cigar(alignment.cigartuples, alignment.reference_start)

    # Get alignment region
    ref_start = alignment.reference_start
    ref_end = alignment.reference_end
    read_length = alignment.query_length

    # Create figure
    fig, ax = plt.subplots(figsize=(20, 5))

    # Track positions
    y_ref = 1.0
    y_read = 0.0
    track_height = 0.25

    # Draw reference and read backbone (light gray)
    ax.plot([0, ref_end - ref_start], [y_ref, y_ref], 'lightgray', linewidth=4, alpha=0.3, zorder=1)
    ax.plot([0, read_length], [y_read, y_read], 'lightgray', linewidth=4, alpha=0.3, zorder=1)

    # First, draw ALL monomers on READ track (colored by family)
    # Each monomer is drawn separately, even if consecutive monomers have same family
    for _, mon in read_monomers.iterrows():
        start = mon['monomer_start']
        end = mon['monomer_end']
        family = mon['monomer_family']

        if not pd.isna(family):
            color = FAMILY_COLORS.get(int(family), 'gray')
            # Draw monomer on read track - SOLID, not transparent
            ax.add_patch(Rectangle((start, y_read - track_height/2),
                                   end - start, track_height,
                                   facecolor=color, edgecolor='black',
                                   linewidth=0.8, alpha=1.0, zorder=3))
        else:
            # Unclassified - draw in visible gray - SOLID
            ax.add_patch(Rectangle((start, y_read - track_height/2),
                                   end - start, track_height,
                                   facecolor='#CCCCCC', edgecolor='#666666',
                                   linewidth=0.8, alpha=1.0, zorder=3))

    # Create read position to reference position mapping
    read_to_ref = {}
    for block in blocks:
        if block['type'] == 'M':
            # Match regions map 1:1
            ref_offset = block['ref_start'] - ref_start
            for i in range(block['read_end'] - block['read_start']):
                read_to_ref[block['read_start'] + i] = ref_offset + i

    # Draw monomers on REFERENCE track (projected from read monomers)
    # Each monomer is drawn separately with visible borders
    for _, mon in read_monomers.iterrows():
        start = mon['monomer_start']
        end = mon['monomer_end']
        family = mon['monomer_family']

        # Map monomer to reference coordinates (if possible)
        ref_positions = []
        for pos in range(start, end):
            if pos in read_to_ref:
                ref_positions.append(read_to_ref[pos])

        if len(ref_positions) > 0:
            ref_s = min(ref_positions)
            ref_e = max(ref_positions) + 1

            if not pd.isna(family):
                color = FAMILY_COLORS.get(int(family), 'gray')
                # Draw monomer on reference track - SOLID
                ax.add_patch(Rectangle((ref_s, y_ref - track_height/2),
                                       ref_e - ref_s, track_height,
                                       facecolor=color, edgecolor='black',
                                       linewidth=0.8, alpha=1.0, zorder=3))
            else:
                # Unclassified - draw in visible gray - SOLID
                ax.add_patch(Rectangle((ref_s, y_ref - track_height/2),
                                       ref_e - ref_s, track_height,
                                       facecolor='#CCCCCC', edgecolor='#666666',
                                       linewidth=0.8, alpha=1.0, zorder=3))

    # Draw deletion monomers on REFERENCE track (from deletion analysis)
    if deletion_monomers is not None:
        # Filter deletion monomers for this read
        # Deletion monomer IDs format: readid_del#_Chr#:start-end_array0_mon0
        read_deletion_monomers = deletion_monomers[
            deletion_monomers['monomer_id'].str.startswith(read_id + '_del')
        ]

        if len(read_deletion_monomers) > 0:
            print(f"  Loaded {len(read_deletion_monomers)} deletion monomers")

        for _, del_mon in read_deletion_monomers.iterrows():
            # Parse the seq_id to get deletion coordinates
            # Format: readid_del#_Chr#:start-end
            seq_id = del_mon['seq_id']
            if '_del' in seq_id and ':' in seq_id:
                coord_part = seq_id.split('_del')[1].split('_')[1]  # Get "Chr3:14727812-14729948"
                chrom, pos_range = coord_part.split(':')
                del_ref_start, del_ref_end = map(int, pos_range.split('-'))

                # Get monomer position within deletion
                mon_offset = del_mon['monomer_start']
                mon_length = del_mon['monomer_end'] - del_mon['monomer_start']

                # Calculate absolute reference position
                abs_ref_start = del_ref_start + mon_offset
                abs_ref_pos = abs_ref_start - ref_start  # Normalize to plot coordinates

                # Get family if available
                family = del_mon.get('monomer_family')

                if not pd.isna(family):
                    color = FAMILY_COLORS.get(int(family), 'gray')
                else:
                    # Unclassified deletion monomer
                    color = '#FFCCCC'  # Light pink to distinguish from read monomers

                # Draw monomer on reference track (solid, not transparent)
                ax.add_patch(Rectangle((abs_ref_pos, y_ref - track_height/2),
                                       mon_length, track_height,
                                       facecolor=color, edgecolor='darkred',
                                       linewidth=0.8, alpha=1.0, zorder=4))

    # Now overlay indels on top (only those >= 100bp)
    MIN_SV_SIZE = 100

    for block in blocks:
        block_type = block['type']

        if block_type == 'I':
            # INSERTION - highlight with darker border
            read_s = block['read_start']
            read_e = block['read_end']
            size = read_e - read_s

            # Only show if >= 100bp
            if size < MIN_SV_SIZE:
                continue

            # Draw subtle highlight for insertion (thin line on top)
            ax.plot([read_s, read_e], [y_read + track_height/2, y_read + track_height/2],
                   color='steelblue', linewidth=1.5, alpha=0.7, zorder=6)

            # Label as INS (smaller, lighter)
            ax.text((read_s + read_e) / 2, y_read - track_height/2 - 0.08,
                   f'{size}bp', fontsize=6, ha='center', va='top',
                   color='steelblue', alpha=0.7, zorder=7)

        elif block_type == 'D':
            # DELETION - highlight on reference
            ref_s = block['ref_start'] - ref_start
            ref_e = block['ref_end'] - ref_start
            size = ref_e - ref_s

            # Only show if >= 100bp
            if size < MIN_SV_SIZE:
                continue

            # Draw subtle highlight for deletion (thin line on top)
            ax.plot([ref_s, ref_e], [y_ref + track_height/2, y_ref + track_height/2],
                   color='indianred', linewidth=1.5, alpha=0.7, zorder=6)

            # Label as DEL (smaller, lighter)
            ax.text((ref_s + ref_e) / 2, y_ref + track_height/2 + 0.08,
                   f'{size}bp', fontsize=6, ha='center', va='bottom',
                   color='indianred', alpha=0.7, zorder=7)

    # Draw gray ribbons ONLY for aligned regions (match blocks)
    # Gaps show where indels are
    for block in blocks:
        if block['type'] == 'M':
            ref_s = block['ref_start'] - ref_start
            ref_e = block['ref_end'] - ref_start
            read_s = block['read_start']
            read_e = block['read_end']

            # Draw simple gray ribbon connecting aligned regions
            from matplotlib.patches import Polygon
            poly = Polygon([
                [ref_s, y_ref - track_height/2],
                [ref_e, y_ref - track_height/2],
                [read_e, y_read + track_height/2],
                [read_s, y_read + track_height/2]
            ], facecolor='darkgray', edgecolor='none', alpha=0.4, zorder=2)
            ax.add_patch(poly)

    # Labels
    ax.text(-read_length * 0.015, y_ref, 'Reference', ha='right', va='center',
           fontsize=12, weight='bold')
    ax.text(-read_length * 0.015, y_read, 'Read', ha='right', va='center',
           fontsize=12, weight='bold')

    # Title
    chrom = read_svs.iloc[0]['chromosome']
    n_insertions = sum(1 for b in blocks if b['type'] == 'I' and (b['read_end'] - b['read_start']) >= MIN_SV_SIZE)
    n_deletions = sum(1 for b in blocks if b['type'] == 'D' and (b['ref_end'] - b['ref_start']) >= MIN_SV_SIZE)
    ax.set_title(f'{read_id}\n{chrom}:{ref_start:,}-{ref_end:,} bp | Read: {read_length:,} bp | Large indels (≥100bp): {n_insertions} INS, {n_deletions} DEL',
                fontsize=11, pad=10)

    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='none', edgecolor='blue', linewidth=2, linestyle='--',
                      label='Insertion (dashed blue box)'),
        mpatches.Patch(facecolor='none', edgecolor='red', linewidth=2, linestyle='--',
                      label='Deletion (dashed red box)'),
    ]

    # Add family colors present in this read
    families_present = read_monomers['monomer_family'].dropna().unique()
    if len(families_present) > 0:
        legend_elements.append(mpatches.Patch(facecolor='white', label='CEN178 Families:'))
        for fam in sorted([int(f) for f in families_present]):
            color = FAMILY_COLORS.get(fam, 'gray')
            legend_elements.append(
                mpatches.Patch(facecolor=color, edgecolor='black',
                             label=f'Family {fam}')
            )

    ax.legend(handles=legend_elements, loc='upper right', fontsize=9, ncol=2,
             framealpha=0.95)

    # Formatting
    ax.set_xlim(-read_length * 0.03, max(ref_end - ref_start, read_length) * 1.02)
    ax.set_ylim(-0.5, 1.5)
    ax.set_xlabel('Position (bp)', fontsize=11)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Visualize all monomers with indels highlighted')
    parser.add_argument('--bam', required=True, help='BAM file')
    parser.add_argument('--sv-info', required=True, help='SV info TSV')
    parser.add_argument('--monomers', required=True, help='Monomer classifications TSV')
    parser.add_argument('--read-ids', nargs='+', required=True,
                       help='Specific read IDs to plot (space-separated)')
    parser.add_argument('--output-dir', default='CEN178profiler/results_v2_test',
                       help='Output directory')
    parser.add_argument('--deletion-monomers', default=None,
                       help='Combined deletion monomers TSV file (optional)')

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Plot each read
    for i, read_id in enumerate(args.read_ids, 1):
        output_file = output_dir / f'indel_families_{i}_{read_id[:8]}.png'
        print(f"\n[{i}/{len(args.read_ids)}] Processing {read_id}...")
        plot_indel_families_with_all_monomers(args.bam, args.sv_info, args.monomers, read_id, output_file, args.deletion_monomers)

    print(f"\n✅ Created {len(args.read_ids)} indel family plots in {output_dir}")

if __name__ == '__main__':
    main()
