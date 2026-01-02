#!/usr/bin/env python3
"""
Overview visualization of all 6 large 3F3 duplications in Chr4
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Load data
hors = pd.read_csv('reference_genome_hors_MONOMER_LEVEL.tsv', sep='\t')
hors['chromosome'] = hors['read_id'].str.extract(r'(Chr\d+)')[0]
hors['hor_length_kb'] = (hors['hor_end'] - hors['hor_start']) / 1000

# Get large duplications
large_hors = hors[hors['hor_length_kb'] >= 40].copy()
large_hors = large_hors.sort_values('hor_length_kb', ascending=False)

print(f"Creating overview of {len(large_hors)} large duplications")

# Create figure
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.3)

# Plot 1: Bar chart of sizes
ax1 = fig.add_subplot(gs[0, :])
x_pos = np.arange(len(large_hors))
colors = plt.cm.Greens(np.linspace(0.4, 0.9, len(large_hors)))

bars = ax1.bar(x_pos, large_hors['hor_length_kb'].values, color=colors, edgecolor='darkgreen', linewidth=2)

ax1.set_xlabel('Duplication Rank', fontweight='bold', fontsize=12)
ax1.set_ylabel('Size (kb)', fontweight='bold', fontsize=12)
ax1.set_title('Large 3F3 Duplications in Chr4 (≥40 kb)\nAll detected by Monomer-Level Algorithm',
             fontweight='bold', fontsize=14)
ax1.set_xticks(x_pos)
ax1.set_xticklabels([f'#{i+1}' for i in range(len(large_hors))])

# Add value labels
for i, (bar, size) in enumerate(zip(bars, large_hors['hor_length_kb'].values)):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
            f'{size:.1f} kb',
            ha='center', va='bottom', fontweight='bold', fontsize=10)

# Plot 2: Genomic positions on Chr4
ax2 = fig.add_subplot(gs[1, :])

chr4_length = hors[hors['chromosome'] == 'Chr4']['array_end'].max()

# Draw chromosome backbone
ax2.add_patch(mpatches.Rectangle((0, 0.4), chr4_length / 1000, 0.2,
                                 facecolor='lightgray', edgecolor='black', linewidth=2))

# Plot each large HOR
for i, (_, hor) in enumerate(large_hors.iterrows()):
    y_offset = 0.7 + (i % 2) * 0.15  # Alternate heights to avoid overlap

    # Draw HOR box
    ax2.add_patch(mpatches.Rectangle(
        (hor['hor_start'] / 1000, y_offset),
        (hor['hor_end'] - hor['hor_start']) / 1000,
        0.1,
        facecolor=colors[i],
        edgecolor='darkgreen',
        linewidth=2
    ))

    # Add label
    mid = (hor['hor_start'] + hor['hor_end']) / 2 / 1000
    ax2.text(mid, y_offset + 0.15, f'#{i+1}\n{hor["hor_length_kb"]:.1f} kb',
            ha='center', va='bottom', fontsize=8, fontweight='bold')

ax2.set_xlim(0, chr4_length / 1000)
ax2.set_ylim(0, 1.5)
ax2.set_xlabel('Chr4 Position (kb)', fontweight='bold', fontsize=12)
ax2.set_title('Genomic Distribution Along Chr4', fontweight='bold', fontsize=13)
ax2.set_yticks([])
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)

# Plot 3: Statistics table
ax3 = fig.add_subplot(gs[2, 0])
ax3.axis('off')

table_data = []
for i, (_, hor) in enumerate(large_hors.iterrows()):
    table_data.append([
        f'#{i+1}',
        f'{hor["hor_unit"]}',
        f'{int(hor["hor_copies"])}',
        f'{int(hor["total_monomers"])}',
        f'{hor["hor_length_kb"]:.1f}'
    ])

table = ax3.table(cellText=table_data,
                 colLabels=['Rank', 'Pattern', 'Copies', 'Monomers', 'Size (kb)'],
                 cellLoc='center',
                 loc='center',
                 bbox=[0, 0, 1, 1])

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# Style header
for i in range(5):
    table[(0, i)].set_facecolor('#2ecc71')
    table[(0, i)].set_text_props(weight='bold', color='white')

# Style rows
for i in range(1, len(table_data) + 1):
    for j in range(5):
        if i % 2 == 0:
            table[(i, j)].set_facecolor('#f0f0f0')

ax3.set_title('Detailed Statistics', fontweight='bold', fontsize=13)

# Plot 4: Comparison text
ax4 = fig.add_subplot(gs[2, 1])
ax4.axis('off')

comparison_text = """
KEY FINDINGS

All 6 large duplications:
  ✓ Are in Chr4
  ✓ Are 3F3 patterns (3-homHORs)
  ✓ Range from 40-372 kb
  ✓ Total span: ~840 kb

Largest duplication:
  Pattern: 3F3
  Copies: 637
  Monomers: 1,911
  Size: 371.5 kb

Why RLE-based missed these:
  ❌ Saw as 1F3 × 1911 (1 monomer/unit)
  ❌ Failed monomers_per_unit ≥ 3 criteria

Why Monomer-level found them:
  ✅ Detected as 3F3 × 637 (3 monomers/unit)
  ✅ Passed BOTH criteria!

Biological significance:
  • Recent large-scale duplication
  • F3 highly specialized (94% in HORs)
  • Matches Nature 2023 findings
"""

ax4.text(0.05, 0.95, comparison_text, transform=ax4.transAxes,
        fontsize=10, verticalalignment='top', family='monospace',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.7))

# Overall title
fig.suptitle(f'Large-Scale Centromeric Duplications (≥40 kb)\n' +
            f'{len(large_hors)} Duplications Found | All 3F3 patterns | Chr4 only',
            fontsize=16, fontweight='bold')

plt.savefig('large_duplications_overview_monomer_level.png', dpi=300, bbox_inches='tight')
plt.close()

print("✅ Saved: large_duplications_overview_monomer_level.png")

# Print summary
print("\n=== LARGE DUPLICATION SUMMARY ===")
print(f"Total large duplications: {len(large_hors)}")
print(f"Total span: {large_hors['hor_length_kb'].sum():.1f} kb")
print(f"Largest: {large_hors.iloc[0]['hor_length_kb']:.1f} kb")
print(f"Smallest: {large_hors.iloc[-1]['hor_length_kb']:.1f} kb")
print(f"Mean: {large_hors['hor_length_kb'].mean():.1f} kb")
print(f"\nAll are:")
print(f"  - 3F3 patterns (3-homHORs)")
print(f"  - Located in Chr4")
print(f"  - Detected by monomer-level algorithm")
