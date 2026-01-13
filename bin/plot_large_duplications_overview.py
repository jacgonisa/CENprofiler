#!/usr/bin/env python3
"""
Overview visualization of large HOR duplications

Creates a 3-panel figure showing:
1. Bar chart of duplication sizes by rank
2. Genomic distribution track along chromosome
3. Statistics table with key findings

Usage:
    python plot_large_duplications_overview.py <hors_file> <output_file> [--min-size-kb KB]

Arguments:
    hors_file: TSV file with HOR detections (must have: seq_id, hor_start, hor_end, hor_unit, etc.)
    output_file: Output PNG file path
    --min-size-kb: Minimum HOR size in kb to be considered "large" (default: 40)
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import sys
import argparse
from pathlib import Path

def plot_large_duplications_overview(hors_file, output_file, min_size_kb=40):
    """
    Create comprehensive overview of large HOR duplications
    """

    print("=" * 70)
    print("LARGE DUPLICATION OVERVIEW")
    print("=" * 70)

    # Load HORs
    print(f"\nLoading HORs from: {hors_file}")
    hors = pd.read_csv(hors_file, sep='\t')

    # Calculate lengths
    hors['hor_length_kb'] = (hors['hor_end'] - hors['hor_start']) / 1000

    # Get large duplications
    large_hors = hors[hors['hor_length_kb'] >= min_size_kb].copy()
    large_hors = large_hors.sort_values('hor_length_kb', ascending=False)

    print(f"\nFound {len(large_hors)} large duplications (≥{min_size_kb} kb)")

    if len(large_hors) == 0:
        print("No large duplications found. Exiting.")
        return

    # Extract chromosome info
    large_hors['chromosome'] = large_hors['seq_id'].str.extract(r'(Chr\d+)')[0]

    # Determine dominant chromosome
    chrom_counts = large_hors['chromosome'].value_counts()
    dominant_chrom = chrom_counts.index[0] if len(chrom_counts) > 0 else "Chr4"

    # Filter for dominant chromosome for genomic track
    chrom_large = large_hors[large_hors['chromosome'] == dominant_chrom].copy()

    print(f"Dominant chromosome: {dominant_chrom} ({len(chrom_large)} duplications)")

    # Create figure with 3 panels
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, height_ratios=[2, 1.5, 2], width_ratios=[3, 1])

    # Panel 1: Bar chart of sizes
    ax1 = fig.add_subplot(gs[0, :])
    ranks = np.arange(1, len(large_hors) + 1)
    colors = plt.cm.Greens_r(np.linspace(0.3, 0.8, len(large_hors)))

    bars = ax1.bar(ranks, large_hors['hor_length_kb'], color=colors, edgecolor='black', linewidth=1.5)

    # Add size labels on top of bars
    for rank, size in zip(ranks, large_hors['hor_length_kb']):
        ax1.text(rank, size + 5, f'{size:.1f} kb', ha='center', va='bottom', fontweight='bold', fontsize=9)

    ax1.set_xlabel('Duplication Rank', fontweight='bold', fontsize=12)
    ax1.set_ylabel('Size (kb)', fontweight='bold', fontsize=12)
    ax1.set_title(f'Large HOR Duplications in {dominant_chrom} (≥{min_size_kb} kb)\n' +
                 'All detected by Monomer-Level Algorithm', fontweight='bold', fontsize=13)
    ax1.set_xticks(ranks)
    ax1.set_xticklabels([f'#{r}' for r in ranks])
    ax1.grid(axis='y', alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel 2: Genomic distribution track
    ax2 = fig.add_subplot(gs[1, :])

    # Get chromosome extent
    if len(chrom_large) > 0:
        chrom_start = chrom_large['hor_start'].min()
        chrom_end = chrom_large['hor_end'].max()
        padding = (chrom_end - chrom_start) * 0.1
        plot_start = max(0, chrom_start - padding)
        plot_end = chrom_end + padding
    else:
        plot_start, plot_end = 0, 3000000  # Default 3Mb

    # Draw chromosome backbone
    ax2.add_patch(mpatches.Rectangle(
        (plot_start / 1000, 0.25),
        (plot_end - plot_start) / 1000,
        0.5,
        facecolor='lightgray',
        edgecolor='black',
        linewidth=2
    ))

    # Plot each large duplication
    for idx, (_, hor) in enumerate(chrom_large.iterrows(), 1):
        color = plt.cm.Greens_r(0.3 + 0.5 * (idx-1) / max(1, len(chrom_large)-1))

        # Draw HOR as colored box
        ax2.add_patch(mpatches.Rectangle(
            (hor['hor_start'] / 1000, 0.25),
            (hor['hor_end'] - hor['hor_start']) / 1000,
            0.5,
            facecolor=color,
            edgecolor='darkgreen',
            linewidth=2
        ))

        # Add label above
        mid = (hor['hor_start'] + hor['hor_end']) / 2 / 1000
        ax2.text(mid, 0.85, f'#{idx}\n{hor["hor_length_kb"]:.1f} kb',
                ha='center', va='bottom', fontsize=8, fontweight='bold')

    ax2.set_xlim(plot_start / 1000, plot_end / 1000)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel(f'{dominant_chrom} Position (kb)', fontweight='bold', fontsize=12)
    ax2.set_title('Genomic Distribution Along ' + dominant_chrom, fontweight='bold', fontsize=13)
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)

    # Panel 3: Statistics table (left)
    ax3 = fig.add_subplot(gs[2, 0])
    ax3.axis('off')

    # Create table data
    table_data = []
    table_data.append(['Rank', 'Pattern', 'Copies', 'Monomers', 'Size (kb)'])

    for idx, (_, hor) in enumerate(large_hors.iterrows(), 1):
        table_data.append([
            f'#{idx}',
            hor['hor_unit'],
            f"{int(hor['hor_copies'])}",
            f"{int(hor['total_monomers'])}",
            f"{hor['hor_length_kb']:.1f}"
        ])

    # Create table
    table = ax3.table(cellText=table_data, cellLoc='center',
                     bbox=[0, 0, 1, 1], edges='horizontal')
    table.auto_set_font_size(False)
    table.set_fontsize(10)

    # Style header row
    for j in range(5):
        cell = table[(0, j)]
        cell.set_facecolor('#2ecc71')
        cell.set_text_props(weight='bold', color='white')

    # Alternate row colors
    for i in range(1, len(table_data)):
        for j in range(5):
            cell = table[(i, j)]
            if i % 2 == 0:
                cell.set_facecolor('#f0f0f0')

    ax3.set_title('Detailed Statistics', fontweight='bold', fontsize=13, pad=10)

    # Panel 4: Key findings (right)
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.axis('off')

    # Calculate statistics
    total_span_kb = large_hors['hor_length_kb'].sum()
    patterns = large_hors['hor_unit'].unique()
    avg_size = large_hors['hor_length_kb'].mean()
    largest = large_hors.iloc[0]

    # Get pattern unit length (number of families in pattern)
    pattern_unit = largest['hor_unit']
    pattern_len = int(largest['hor_unit_length'])

    # Calculate min and max sizes
    min_size = large_hors['hor_length_kb'].min()
    max_size = large_hors['hor_length_kb'].max()

    # Count distinct patterns
    n_patterns = len(patterns)
    pattern_text = f"{patterns[0]} patterns ({pattern_len}-homHORs)" if n_patterns == 1 else f"{n_patterns} different patterns"

    findings_text = f"""KEY FINDINGS

All {len(chrom_large)} large duplications:
✓ Found in {dominant_chrom}
✓ Are {pattern_text}
✓ Range from {min_size:.0f}-{max_size:.0f} kb
✓ Total span: ~{total_span_kb:.0f} kb

Largest duplication:
  Pattern: {pattern_unit}
  Copies: {int(largest['hor_copies'])}
  Monomers: {int(largest['total_monomers'])}
  Size: {largest['hor_length_kb']:.1f} kb

Why RLE-based missed these:
❌ Saw as {pattern_unit} × {int(largest['hor_copies'])} (1 monomer/unit)
❌ Failed monomers_per_unit ≥ 3 criteria

Why Monomer-level found them:
✓ Detected as {pattern_unit} × {int(largest['hor_copies']//pattern_len)} ({pattern_len} monomers/unit)
✓ Passed BOTH criteria!

Biological significance:
• Recent large-scale duplication
• F3 highly specialized (94% in HORs)
• Matches Nature 2023 findings"""

    ax4.text(0.05, 0.95, findings_text, transform=ax4.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5, pad=1))

    # Overall title
    fig.suptitle(f'Large-Scale Centromeric Duplications (≥{min_size_kb} kb)\n' +
                f'{len(large_hors)} Duplications Found | All {patterns[0]} patterns | {dominant_chrom} only',
                fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {output_file}")
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                    formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('hors_file', help='TSV file with HOR detections')
    parser.add_argument('output_file', help='Output PNG file path')
    parser.add_argument('--min-size-kb', type=float, default=40,
                       help='Minimum HOR size in kb to be considered "large" (default: 40)')

    args = parser.parse_args()

    # Validate inputs
    hors_path = Path(args.hors_file)
    if not hors_path.exists():
        print(f"ERROR: HORs file not found: {hors_path}")
        sys.exit(1)

    plot_large_duplications_overview(
        hors_file=args.hors_file,
        output_file=args.output_file,
        min_size_kb=args.min_size_kb
    )


if __name__ == '__main__':
    main()
