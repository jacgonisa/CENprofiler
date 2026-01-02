#!/usr/bin/env python3
"""
Detailed visualization of large 3F3 duplications in Chr4
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# Load data
hors = pd.read_csv('reference_genome_hors_MONOMER_LEVEL.tsv', sep='\t')
hors['chromosome'] = hors['read_id'].str.extract(r'(Chr\d+)')[0]
hors['hor_length_kb'] = (hors['hor_end'] - hors['hor_start']) / 1000

classified = pd.read_csv('../monomer_classifications.tsv', sep='\t')

# Get large duplications
large_hors = hors[hors['hor_length_kb'] >= 40].copy()
large_hors = large_hors.sort_values('hor_length_kb', ascending=False)

print(f"Found {len(large_hors)} large duplications (≥40 kb)")
print(large_hors[['chromosome', 'hor_unit', 'hor_copies', 'total_monomers', 'hor_length_kb']])

# Family colors
family_colors = {
    1: '#e74c3c', 2: '#3498db', 3: '#2ecc71', 4: '#f39c12',
    5: '#9b59b6', 6: '#1abc9c', 7: '#e67e22', 8: '#34495e',
    9: '#16a085', 10: '#27ae60', 11: '#2980b9', 12: '#8e44ad',
    13: '#c0392b', 14: '#d35400', 15: '#7f8c8d', 16: '#2c3e50',
    17: '#f1c40f', 18: '#95a5a6', 19: '#ecf0f1', 20: '#bdc3c7'
}

# Select top 3 for detailed visualization
top_3 = large_hors.head(3)

for idx, large_hor in top_3.iterrows():
    print(f"\nProcessing: {large_hor['hor_unit']} at {large_hor['read_id']}")

    # Get region info
    read_id = large_hor['read_id']
    hor_start = large_hor['hor_start']
    hor_end = large_hor['hor_end']

    # Add padding (50kb on each side to show context)
    padding = 50000
    region_start = max(0, hor_start - padding)
    region_end = hor_end + padding

    # Get monomers in this region
    region_monomers = classified[
        (classified['read_id'] == read_id) &
        (classified['monomer_start'] >= region_start) &
        (classified['monomer_end'] <= region_end)
    ].copy().sort_values('monomer_start')

    # Get ALL HORs in this region
    region_hors = hors[
        (hors['read_id'] == read_id) &
        (hors['hor_start'] < region_end) &
        (hors['hor_end'] > region_start)
    ].copy()

    # Create figure
    fig, axes = plt.subplots(3, 1, figsize=(24, 9),
                            gridspec_kw={'height_ratios': [1, 1, 0.5]})

    # Track 1: Individual monomers
    ax1 = axes[0]
    for _, mono in region_monomers.iterrows():
        if pd.notna(mono['monomer_family']):
            fam = int(mono['monomer_family'])
            color = family_colors.get(fam, '#95a5a6')
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

    # Add label
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

    # Title
    fig.suptitle(f'Large 3F3 Duplication - Chr4\n' +
                f'{large_hor["hor_unit"]} repeated {int(large_hor["hor_copies"])} times = ' +
                f'{int(large_hor["total_monomers"])} monomers ({large_hor["hor_length_kb"]:.1f} kb)',
                fontsize=14, fontweight='bold')

    plt.tight_layout()

    # Save
    rank = list(top_3.index).index(idx) + 1
    filename = f'large_duplication_Chr4_rank{rank}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  ✅ Saved: {filename}")

print("\n✅ All large duplication detail plots generated!")
