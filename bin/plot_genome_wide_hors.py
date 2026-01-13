#!/usr/bin/env python3
"""
Genome-wide HOR distribution visualization

Creates a multi-track figure showing HOR distribution across all chromosomes.
Each chromosome has:
1. Individual HOR boxes colored by type (homHOR/hetHOR)
2. Coverage density plot showing HOR-rich regions

Usage:
    python plot_genome_wide_hors.py <hors_file> <output_file>

Arguments:
    hors_file: TSV file with HOR detections (must have: seq_id, hor_start, hor_end, hor_type, etc.)
    output_file: Output PNG file path
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys
import argparse
from pathlib import Path


def plot_genome_wide_hors(hors_file, output_file):
    """
    Create genome-wide HOR distribution plot
    """

    print("=" * 70)
    print("GENOME-WIDE HOR DISTRIBUTION")
    print("=" * 70)

    # Load HORs
    print(f"\nLoading HORs from: {hors_file}")
    hors = pd.read_csv(hors_file, sep='\t')

    # Extract chromosome from seq_id
    hors['chromosome'] = hors['seq_id'].str.extract(r'(Chr\d+)')[0]

    # Filter for valid chromosomes
    valid_chroms = [f'Chr{i}' for i in range(1, 6)]
    hors = hors[hors['chromosome'].isin(valid_chroms)].copy()

    print(f"Total HORs: {len(hors)}")
    print(f"Chromosomes: {sorted(hors['chromosome'].unique())}")

    # Chromosome colors
    chr_colors = {
        'Chr1': '#e74c3c',
        'Chr2': '#3498db',
        'Chr3': '#2ecc71',
        'Chr4': '#f39c12',
        'Chr5': '#9b59b6'
    }

    # HOR type colors
    colors_hor_type = {
        'homHOR': '#e74c3c',
        'hetHOR': '#3498db'
    }

    # Create figure with 6 panels (5 chromosomes + legend)
    fig = plt.figure(figsize=(26, 14))
    gs = fig.add_gridspec(6, 1, height_ratios=[1, 1, 1, 1, 1, 0.3], hspace=0.4)

    # Plot each chromosome
    chromosomes = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']

    for chr_idx, chrom in enumerate(chromosomes):
        chr_hors = hors[hors['chromosome'] == chrom].copy()

        ax = fig.add_subplot(gs[chr_idx, 0])

        # Get chromosome length from data
        if len(chr_hors) > 0:
            chr_length = chr_hors['hor_end'].max()
        else:
            chr_length = 1000000  # Default 1Mb

        # Track 1: Individual HOR boxes colored by type
        y_pos_type = 0.7
        for _, hor in chr_hors.iterrows():
            color = colors_hor_type.get(hor['hor_type'], '#95a5a6')
            ax.add_patch(mpatches.Rectangle(
                (hor['hor_start'] / 1000, y_pos_type),
                (hor['hor_end'] - hor['hor_start']) / 1000,
                0.25,
                facecolor=color,
                edgecolor='black',
                linewidth=0.5,
                alpha=0.9
            ))

        # Track 2: Coverage density plot
        y_pos_cov = 0.3
        if len(chr_hors) > 0:
            # Create coverage array
            coverage = np.zeros(int(chr_length) + 1)
            for _, hor in chr_hors.iterrows():
                start_idx = int(hor['hor_start'])
                end_idx = int(hor['hor_end'])
                coverage[start_idx:end_idx] = 1

            # Smooth coverage for visualization (10kb windows)
            window_size = 10000
            smoothed = np.convolve(coverage, np.ones(window_size)/window_size, mode='same')

            # Plot as filled area
            x_coords = np.arange(len(smoothed)) / 1000  # Convert to kb
            y_coords = y_pos_cov + smoothed * 0.25
            ax.fill_between(x_coords, y_pos_cov, y_coords,
                            color=chr_colors[chrom], alpha=0.6)

        # Formatting
        ax.set_xlim(0, chr_length / 1000)
        ax.set_ylim(0, 1)
        ax.set_ylabel(chrom, fontsize=14, fontweight='bold', rotation=0,
                     labelpad=40, va='center')
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        if chr_idx == len(chromosomes) - 1:
            ax.set_xlabel('Position (kb)', fontsize=12, fontweight='bold')
        else:
            ax.set_xticks([])
            ax.spines['bottom'].set_visible(False)

        # Add HOR count summary
        homhor_count = len(chr_hors[chr_hors['hor_type'] == 'homHOR'])
        hethor_count = len(chr_hors[chr_hors['hor_type'] == 'hetHOR'])
        ax.text(0.02, 0.95, f'n={len(chr_hors)} ({homhor_count} homHOR, {hethor_count} hetHOR)',
               transform=ax.transAxes, fontsize=11,
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
               verticalalignment='top')

    # Legend panel
    ax_legend = fig.add_subplot(gs[5, 0])
    ax_legend.axis('off')

    homhor_total = len(hors[hors['hor_type'] == 'homHOR'])
    hethor_total = len(hors[hors['hor_type'] == 'hetHOR'])

    legend_elements = [
        mpatches.Patch(facecolor=colors_hor_type['homHOR'],
                      edgecolor='black', label=f'homHOR (n={homhor_total})', alpha=0.9),
        mpatches.Patch(facecolor=colors_hor_type['hetHOR'],
                      edgecolor='black', label=f'hetHOR (n={hethor_total})', alpha=0.9)
    ]

    ax_legend.legend(handles=legend_elements, loc='center', ncol=2,
                    frameon=True, fontsize=12, title='HOR Type')

    # Title
    n_unique = hors['hor_unit'].nunique()
    fig.suptitle('Genome-Wide HOR Distribution - MONOMER-LEVEL Detection\n' +
                f'(min_copies ≥ 3 AND monomers_per_unit ≥ 3)\n' +
                f'Total: {len(hors)} HORs | {n_unique} unique patterns | {homhor_total} homHORs + {hethor_total} hetHORs',
                fontsize=16, fontweight='bold', y=0.98)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {output_file}")
    print("=" * 70)

    # Print statistics
    print("\n=== STATISTICS ===")
    print(f"Total HORs: {len(hors)}")
    print(f"Unique patterns: {n_unique}")
    print(f"\nPer chromosome:")
    for chrom in chromosomes:
        chr_hors = hors[hors['chromosome'] == chrom]
        homhor = len(chr_hors[chr_hors['hor_type'] == 'homHOR'])
        hethor = len(chr_hors[chr_hors['hor_type'] == 'hetHOR'])
        print(f"  {chrom}: {len(chr_hors):3d} HORs ({homhor:3d} homHOR, {hethor:3d} hetHOR)")

    print(f"\nHOR type distribution:")
    print(f"  homHOR: {homhor_total} ({homhor_total/len(hors)*100:.1f}%)")
    print(f"  hetHOR: {hethor_total} ({hethor_total/len(hors)*100:.1f}%)")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hors_file', help='TSV file with HOR detections')
    parser.add_argument('output_file', help='Output PNG file path')

    args = parser.parse_args()

    # Validate inputs
    hors_path = Path(args.hors_file)
    if not hors_path.exists():
        print(f"ERROR: HORs file not found: {hors_path}")
        sys.exit(1)

    plot_genome_wide_hors(
        hors_file=args.hors_file,
        output_file=args.output_file
    )


if __name__ == '__main__':
    main()
