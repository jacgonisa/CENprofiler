#!/usr/bin/env python3
"""
Genome-wide visualization for monomer-level HOR detection
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Load monomer-level HORs
hors = pd.read_csv('reference_genome_hors_MONOMER_LEVEL.tsv', sep='\t')
# Handle both 'read_id' and 'seq_id' column names
id_col = 'read_id' if 'read_id' in hors.columns else 'seq_id'
hors['chromosome'] = hors[id_col].str.extract(r'(Chr\d+)')[0]

print(f"Loaded {len(hors)} monomer-level HORs")

# Chromosome order and colors
chromosomes = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']
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

# Create figure
fig = plt.figure(figsize=(26, 14))
gs = fig.add_gridspec(6, 1, height_ratios=[1, 1, 1, 1, 1, 0.3], hspace=0.4)

# Track for each chromosome
for chr_idx, chrom in enumerate(chromosomes):
    chr_hors = hors[hors['chromosome'] == chrom].copy()

    ax = fig.add_subplot(gs[chr_idx, 0])

    # Get chromosome length
    chr_length = chr_hors['array_end'].max() if len(chr_hors) > 0 else 1000000

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

    # Track 2: Coverage density
    y_pos_cov = 0.3
    if len(chr_hors) > 0:
        # Create coverage array
        coverage = np.zeros(int(chr_length) + 1)
        for _, hor in chr_hors.iterrows():
            coverage[int(hor['hor_start']):int(hor['hor_end'])] = 1

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

    # Add HOR count
    homhor_count = len(chr_hors[chr_hors['hor_type'] == 'homHOR'])
    hethor_count = len(chr_hors[chr_hors['hor_type'] == 'hetHOR'])
    ax.text(0.02, 0.95, f'n={len(chr_hors)} ({homhor_count} homHOR, {hethor_count} hetHOR)',
           transform=ax.transAxes, fontsize=11,
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
           verticalalignment='top')

# Legend
ax_legend = fig.add_subplot(gs[5, 0])
ax_legend.axis('off')

legend_elements = [
    mpatches.Patch(facecolor=colors_hor_type['homHOR'],
                  edgecolor='black', label=f'homHOR (n={len(hors[hors["hor_type"]=="homHOR"])})', alpha=0.9),
    mpatches.Patch(facecolor=colors_hor_type['hetHOR'],
                  edgecolor='black', label=f'hetHOR (n={len(hors[hors["hor_type"]=="hetHOR"])})', alpha=0.9)
]

ax_legend.legend(handles=legend_elements, loc='center', ncol=2,
                frameon=True, fontsize=12, title='HOR Type')

# Title
fig.suptitle('Genome-Wide HOR Distribution - MONOMER-LEVEL Detection\n(min_copies ≥ 3 AND monomers_per_unit ≥ 3)\n' +
             f'Total: {len(hors)} HORs | {hors["hor_unit"].nunique()} unique patterns | 187 homHORs + 247 hetHORs',
             fontsize=16, fontweight='bold', y=0.98)

plt.savefig('genome_wide_HORs_monomer_level.png', dpi=300, bbox_inches='tight')
plt.close()

print("✅ Saved: genome_wide_HORs_monomer_level.png")

# Statistics
print(f"\n=== MONOMER-LEVEL HOR STATISTICS ===")
print(f"Total HORs: {len(hors)}")
print(f"Unique patterns: {hors['hor_unit'].nunique()}")

print(f"\nPer chromosome:")
for chrom in chromosomes:
    chr_hors = hors[hors['chromosome'] == chrom]
    homhor = len(chr_hors[chr_hors['hor_type'] == 'homHOR'])
    hethor = len(chr_hors[chr_hors['hor_type'] == 'hetHOR'])
    print(f"  {chrom}: {len(chr_hors)} ({homhor} homHOR, {hethor} hetHOR)")

print(f"\nHOR type:")
print(hors['hor_type'].value_counts())
