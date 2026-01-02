#!/usr/bin/env python3
"""
Generate HOR schematic diagrams for monomer-level detection
Shows TOP homHORs AND TOP hetHORs
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# Load data
hors = pd.read_csv('reference_genome_hors_MONOMER_LEVEL.tsv', sep='\t')

# Family colors
family_colors = {
    1: '#e74c3c', 2: '#3498db', 3: '#2ecc71', 4: '#f39c12',
    5: '#9b59b6', 6: '#1abc9c', 7: '#e67e22', 8: '#34495e',
    9: '#16a085', 10: '#27ae60', 11: '#2980b9', 12: '#8e44ad',
    13: '#c0392b', 14: '#d35400', 15: '#7f8c8d', 16: '#2c3e50',
    17: '#f1c40f', 18: '#95a5a6', 19: '#ecf0f1', 20: '#bdc3c7'
}

def parse_hor_unit(hor_unit_str):
    """Parse HOR unit string into list of (count, family) tuples"""
    elements = hor_unit_str.split('-')
    parsed = []
    for elem in elements:
        count, family = elem.split('F')
        parsed.append((int(count), int(family)))
    return parsed

def draw_hor_schematic(ax, hor_unit_str, occurrences, y_pos, is_homhor=True):
    """Draw schematic representation of a HOR pattern"""
    parsed = parse_hor_unit(hor_unit_str)

    # Box parameters
    box_width = 0.4
    spacing = 0.05
    x_start = 1.5

    # Draw individual monomer blocks
    # IMPORTANT: If count > 1, draw multiple separate boxes!
    x_pos = x_start
    for i, (count, family) in enumerate(parsed):
        color = family_colors.get(family, '#95a5a6')

        # Draw 'count' number of boxes for this family
        for c in range(count):
            # Draw box
            box = FancyBboxPatch(
                (x_pos, y_pos), box_width, 0.6,
                boxstyle="round,pad=0.05",
                facecolor=color,
                edgecolor='black',
                linewidth=2
            )
            ax.add_patch(box)

            # Add family label
            ax.text(x_pos + box_width/2, y_pos + 0.3, f'F{family}',
                   ha='center', va='center', fontweight='bold', fontsize=10)

            x_pos += box_width + spacing

    # Draw repeat bracket
    bracket_y = y_pos - 0.15
    bracket_start = x_start - 0.1
    bracket_end = x_pos - spacing + 0.1

    # Horizontal line
    ax.plot([bracket_start, bracket_end], [bracket_y, bracket_y],
           'k-', linewidth=2)

    # Left vertical tick
    ax.plot([bracket_start, bracket_start], [bracket_y, bracket_y + 0.1],
           'k-', linewidth=2)

    # Right vertical tick
    ax.plot([bracket_end, bracket_end], [bracket_y, bracket_y + 0.1],
           'k-', linewidth=2)

    # Add repeat annotation
    bracket_mid = (bracket_start + bracket_end) / 2
    ax.text(bracket_mid, bracket_y - 0.15, f'× many copies',
           ha='center', va='top', fontweight='bold', fontsize=9, style='italic')

    # Add occurrence count on the right
    ax.text(x_pos + 0.5, y_pos + 0.3, f'n={occurrences}',
           ha='left', va='center', fontsize=10,
           bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    return x_pos + 1.0

# Separate homHORs and hetHORs
homhors = hors[hors['hor_type'] == 'homHOR']
hethors = hors[hors['hor_type'] == 'hetHOR']

# Get top patterns
top_homhors = homhors['hor_unit'].value_counts().head(10)
top_hethors = hethors['hor_unit'].value_counts().head(10)

print(f"Generating schematics for:")
print(f"  - Top {len(top_homhors)} homHOR patterns")
print(f"  - Top {len(top_hethors)} hetHOR patterns")

# Create two separate figures: one for homHORs, one for hetHORs

# FIGURE 1: homHORs
fig1, ax1 = plt.subplots(figsize=(18, 10))
ax1.set_xlim(0, 15)
ax1.set_ylim(0, len(top_homhors) + 1)
ax1.axis('off')

# Title
ax1.text(7.5, len(top_homhors) + 0.5,
       'homHOR Architecture - Monomer-Level Detection\n(min_copies ≥ 3 AND monomers_per_unit ≥ 3)',
       ha='center', va='bottom', fontsize=16, fontweight='bold')

# Draw each homHOR pattern
for idx, (pattern, occurrences) in enumerate(top_homhors.items()):
    y_pos = len(top_homhors) - idx - 0.5

    # Pattern number
    ax1.text(0.3, y_pos + 0.3, f'{idx+1}.',
           ha='right', va='center', fontsize=12, fontweight='bold')

    # Draw schematic
    draw_hor_schematic(ax1, pattern, occurrences, y_pos, is_homhor=True)

plt.savefig('hor_schematics_homHOR_monomer_level.png', dpi=300, bbox_inches='tight')
plt.close()
print("✅ Saved: hor_schematics_homHOR_monomer_level.png")

# FIGURE 2: hetHORs
fig2, ax2 = plt.subplots(figsize=(18, 10))
ax2.set_xlim(0, 15)
ax2.set_ylim(0, len(top_hethors) + 1)
ax2.axis('off')

# Title
ax2.text(7.5, len(top_hethors) + 0.5,
       'hetHOR Architecture - Monomer-Level Detection\n(min_copies ≥ 3 AND monomers_per_unit ≥ 3)',
       ha='center', va='bottom', fontsize=16, fontweight='bold')

# Draw each hetHOR pattern
for idx, (pattern, occurrences) in enumerate(top_hethors.items()):
    y_pos = len(top_hethors) - idx - 0.5

    # Pattern number
    ax2.text(0.3, y_pos + 0.3, f'{idx+1}.',
           ha='right', va='center', fontsize=12, fontweight='bold')

    # Draw schematic
    draw_hor_schematic(ax2, pattern, occurrences, y_pos, is_homhor=False)

plt.savefig('hor_schematics_hetHOR_monomer_level.png', dpi=300, bbox_inches='tight')
plt.close()
print("✅ Saved: hor_schematics_hetHOR_monomer_level.png")

# Print statistics
print("\n=== HOR PATTERN STATISTICS ===")
print(f"Total monomer-level HORs: {len(hors)}")
print(f"homHORs: {len(homhors)} ({len(homhors)/len(hors)*100:.1f}%)")
print(f"hetHORs: {len(hethors)} ({len(hethors)/len(hors)*100:.1f}%)")

print(f"\nTop 10 homHOR patterns:")
print(top_homhors)

print(f"\nTop 10 hetHOR patterns:")
print(top_hethors)
