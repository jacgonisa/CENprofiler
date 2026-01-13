#!/usr/bin/env python3
"""
Detailed visualization of large HOR duplications showing monomer→HOR calling process

This script creates 3-track plots showing:
1. Individual monomers by family (colored)
2. HOR calls overlaid
3. Position scale

This demonstrates visually how individual monomer families are organized into HORs.

Usage:
    python plot_large_duplication_detail.py <hors_file> <classifications_file> <output_dir> [--top-n N]

Arguments:
    hors_file: TSV file with HOR detections
    classifications_file: TSV file with monomer classifications
    output_dir: Directory to save output plots
    --top-n: Number of top large duplications to plot (default: 3)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys
import argparse
from pathlib import Path

# Family colors (consistent across all plots)
FAMILY_COLORS = {
    1: '#e74c3c', 2: '#3498db', 3: '#2ecc71', 4: '#f39c12',
    5: '#9b59b6', 6: '#1abc9c', 7: '#e67e22', 8: '#34495e',
    9: '#16a085', 10: '#27ae60', 11: '#2980b9', 12: '#8e44ad',
    13: '#c0392b', 14: '#d35400', 15: '#7f8c8d', 16: '#2c3e50',
    17: '#f1c40f', 18: '#95a5a6', 19: '#ecf0f1', 20: '#bdc3c7'
}


def plot_large_duplication_detail(large_hor, hors, classified, output_path, rank):
    """
    Create detailed 3-track visualization for a single large HOR duplication

    Track 1: Individual monomers colored by family
    Track 2: HOR calls overlaid (highlighting the large duplication)
    Track 3: Position axis
    """

    print(f"  Processing rank {rank}: {large_hor['hor_unit']} at {large_hor['seq_id']}")

    # Get region info
    seq_id = large_hor['seq_id']
    hor_start = large_hor['hor_start']
    hor_end = large_hor['hor_end']

    # Add padding (50kb on each side to show context)
    padding = 50000
    region_start = max(0, hor_start - padding)
    region_end = hor_end + padding

    # Get monomers in this region
    region_monomers = classified[
        (classified['seq_id'] == seq_id) &
        (classified['monomer_start'] >= region_start) &
        (classified['monomer_end'] <= region_end)
    ].copy().sort_values('monomer_start')

    # Get ALL HORs in this region
    region_hors = hors[
        (hors['seq_id'] == seq_id) &
        (hors['hor_start'] < region_end) &
        (hors['hor_end'] > region_start)
    ].copy()

    # Create figure with 3 tracks
    fig, axes = plt.subplots(3, 1, figsize=(24, 9),
                            gridspec_kw={'height_ratios': [1, 1, 0.5]})

    # Track 1: Individual monomers
    ax1 = axes[0]
    for _, mono in region_monomers.iterrows():
        if pd.notna(mono['monomer_family']):
            fam = int(mono['monomer_family'])
            color = FAMILY_COLORS.get(fam, '#95a5a6')
            ax1.add_patch(mpatches.Rectangle(
                (mono['monomer_start'] / 1000, 0.2),
                (mono['monomer_end'] - mono['monomer_start']) / 1000,
                0.6,
                facecolor=color,
                edgecolor='black',
                linewidth=0.2
            ))

    ax1.set_xlim(region_start / 1000, region_end / 1000)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Monomers', fontweight='bold')
    ax1.set_yticks([])
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_xticks([])

    # Track 2: HORs
    ax2 = axes[1]

    # Highlight the large HOR
    ax2.add_patch(mpatches.Rectangle(
        (large_hor['hor_start'] / 1000, 0.25),
        (large_hor['hor_end'] - large_hor['hor_start']) / 1000,
        0.5,
        facecolor='#2ecc71',
        edgecolor='darkgreen',
        linewidth=3,
        alpha=0.8,
        label='Large duplication'
    ))

    # Add label to large HOR
    mid = (large_hor['hor_start'] + large_hor['hor_end']) / 2 / 1000
    ax2.text(mid, 0.5,
            f"{large_hor['hor_unit']} × {int(large_hor['hor_copies'])}\n" +
            f"{int(large_hor['total_monomers'])} monomers\n" +
            f"{large_hor['hor_length_kb']:.1f} kb",
            ha='center', va='center', fontsize=9, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    # Show other HORs in region
    for _, hor in region_hors.iterrows():
        if hor['hor_start'] != large_hor['hor_start']:  # Don't redraw the large one
            color = '#e74c3c' if hor['hor_type'] == 'homHOR' else '#3498db'
            ax2.add_patch(mpatches.Rectangle(
                (hor['hor_start'] / 1000, 0.25),
                (hor['hor_end'] - hor['hor_start']) / 1000,
                0.5,
                facecolor=color,
                edgecolor='black',
                linewidth=1,
                alpha=0.5
            ))

    ax2.set_xlim(region_start / 1000, region_end / 1000)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel('HORs', fontweight='bold')
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xticks([])

    # Track 3: Position axis
    ax3 = axes[2]
    ax3.set_xlim(region_start / 1000, region_end / 1000)
    ax3.set_ylim(0, 1)
    ax3.set_xlabel('Position (kb)', fontweight='bold', fontsize=12)
    ax3.set_yticks([])
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_visible(False)

    # Extract chromosome from seq_id
    chrom = seq_id.split('_')[0] if '_' in seq_id else seq_id

    # Title
    fig.suptitle(f'Large HOR Duplication - {chrom}\n' +
                f'{large_hor["hor_unit"]} repeated {int(large_hor["hor_copies"])} times = ' +
                f'{int(large_hor["total_monomers"])} monomers ({large_hor["hor_length_kb"]:.1f} kb)',
                fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"    Saved: {output_path.name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hors_file', help='TSV file with HOR detections')
    parser.add_argument('classifications_file', help='TSV file with monomer classifications')
    parser.add_argument('output_dir', help='Directory to save output plots')
    parser.add_argument('--top-n', type=int, default=3,
                       help='Number of top large duplications to plot (default: 3)')
    parser.add_argument('--min-size-kb', type=float, default=40,
                       help='Minimum HOR size in kb to be considered "large" (default: 40)')

    args = parser.parse_args()

    # Validate inputs
    hors_path = Path(args.hors_file)
    classifications_path = Path(args.classifications_file)
    output_dir = Path(args.output_dir)

    if not hors_path.exists():
        print(f"ERROR: HORs file not found: {hors_path}")
        sys.exit(1)

    if not classifications_path.exists():
        print(f"ERROR: Classifications file not found: {classifications_path}")
        sys.exit(1)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("LARGE DUPLICATION DETAIL PLOTS")
    print("=" * 70)

    # Load data
    print(f"\nLoading HORs from: {hors_path}")
    hors = pd.read_csv(hors_path, sep='\t')

    print(f"Loading classifications from: {classifications_path}")
    classified = pd.read_csv(classifications_path, sep='\t')

    # Calculate HOR lengths
    hors['hor_length_kb'] = (hors['hor_end'] - hors['hor_start']) / 1000

    # Get large duplications
    large_hors = hors[hors['hor_length_kb'] >= args.min_size_kb].copy()
    large_hors = large_hors.sort_values('hor_length_kb', ascending=False)

    print(f"\nFound {len(large_hors)} large duplications (≥{args.min_size_kb} kb)")

    if len(large_hors) == 0:
        print("No large duplications found. Exiting.")
        return

    print(f"\nTop {min(args.top_n, len(large_hors))} largest:")
    print(large_hors.head(args.top_n)[['seq_id', 'hor_unit', 'hor_copies',
                                        'total_monomers', 'hor_length_kb']].to_string(index=False))

    # Plot top N
    top_n = large_hors.head(args.top_n)

    print(f"\nGenerating detail plots for top {len(top_n)} large duplications...")

    for rank, (_, large_hor) in enumerate(top_n.iterrows(), 1):
        # Extract chromosome for filename
        chrom = large_hor['seq_id'].split('_')[0] if '_' in large_hor['seq_id'] else 'region'
        output_path = output_dir / f'large_duplication_{chrom}_rank{rank}.png'

        plot_large_duplication_detail(
            large_hor=large_hor,
            hors=hors,
            classified=classified,
            output_path=output_path,
            rank=rank
        )

    print("\n" + "=" * 70)
    print(f"All {len(top_n)} large duplication detail plots generated!")
    print(f"Output directory: {output_dir}")
    print("=" * 70)


if __name__ == '__main__':
    main()
