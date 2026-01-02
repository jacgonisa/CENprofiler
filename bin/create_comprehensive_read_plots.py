#!/usr/bin/env python3
"""
Comprehensive read visualization module - replicates results_v2 plots

Generates:
1. Family distribution bar chart
2. Family transition heatmap
3. Top arrays combined view
4. Individual array detail plots
5. Summary statistics
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import argparse
from pathlib import Path

# Family color scheme (consistent with visualize_indel_families_v2.py)
FAMILY_COLORS = {
    1: '#e41a1c',   # Red
    2: '#377eb8',   # Blue
    3: '#4daf4a',   # Green
    4: '#984ea3',   # Purple
    5: '#999999',   # Gray
    6: '#cccccc',   # Light gray
    7: '#ff7f00',   # Orange
    8: '#ffff33',   # Yellow
    9: '#a65628',   # Brown
    10: '#f781bf',  # Pink
    11: '#a65628',  # Brown
    12: '#999999',  # Gray
    13: '#f781bf',  # Pink
    14: '#999999',  # Gray
    15: '#999999',  # Gray
    16: '#999999',  # Gray
    17: '#999999',  # Gray
    18: '#00ffff',  # Cyan
    19: '#999999',  # Gray
    20: '#999999',  # Gray
}

def plot_family_summary(monomers_df, output_file):
    """Generate two-panel summary: distribution + transition heatmap"""

    # Filter classified monomers
    classified = monomers_df[monomers_df['monomer_family'].notna()].copy()

    if len(classified) == 0:
        print("No classified monomers for family summary")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

    # Panel A: Family distribution
    family_counts = classified['monomer_family'].value_counts().sort_index()
    families = [int(f) for f in family_counts.index]
    colors = [FAMILY_COLORS.get(f, '#999999') for f in families]

    ax1.bar(range(len(families)), family_counts.values,
            color=colors, edgecolor='black', linewidth=1.5)
    ax1.set_xticks(range(len(families)))
    ax1.set_xticklabels([f'Fam {f}' for f in families], rotation=45, ha='right')
    ax1.set_ylabel('Monomer count', fontweight='bold', fontsize=12)
    ax1.set_title('CEN178 Family Distribution', fontweight='bold', fontsize=14)
    ax1.grid(axis='y', alpha=0.3)

    # Add counts on bars
    for i, (f, count) in enumerate(zip(families, family_counts.values)):
        pct = count / family_counts.sum() * 100
        ax1.text(i, count, f'{count}\n({pct:.1f}%)',
                ha='center', va='bottom', fontsize=9, fontweight='bold')

    # Panel B: Family transition heatmap
    # Group by read and get sequential transitions
    transitions = {}

    for read_id, read_df in classified.groupby('seq_id'):
        # Sort by array_idx and monomer_idx
        read_df = read_df.sort_values(['array_idx', 'monomer_idx'])

        # Get family sequence
        fam_seq = read_df['monomer_family'].astype(int).tolist()

        # Count transitions
        for i in range(len(fam_seq) - 1):
            from_fam = fam_seq[i]
            to_fam = fam_seq[i + 1]
            key = (from_fam, to_fam)
            transitions[key] = transitions.get(key, 0) + 1

    if len(transitions) > 0:
        # Create transition matrix
        all_fams = sorted(set([f for pair in transitions.keys() for f in pair]))
        n_fams = len(all_fams)
        trans_matrix = np.zeros((n_fams, n_fams))

        fam_to_idx = {f: i for i, f in enumerate(all_fams)}
        for (from_fam, to_fam), count in transitions.items():
            i = fam_to_idx[from_fam]
            j = fam_to_idx[to_fam]
            trans_matrix[i, j] = count

        # Plot heatmap
        sns.heatmap(trans_matrix, annot=True, fmt='.0f', cmap='YlOrRd',
                   xticklabels=[f'F{f}' for f in all_fams],
                   yticklabels=[f'F{f}' for f in all_fams],
                   cbar_kws={'label': 'Transition count'},
                   ax=ax2, linewidths=0.5, linecolor='gray')
        ax2.set_xlabel('To family', fontweight='bold', fontsize=12)
        ax2.set_ylabel('From family', fontweight='bold', fontsize=12)
        ax2.set_title('Sequential Family Transitions', fontweight='bold', fontsize=14)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Generated {output_file}")

def plot_top_arrays_combined(monomers_df, output_file, n_arrays=3):
    """Plot top N largest arrays side-by-side showing family composition"""

    classified = monomers_df[monomers_df['monomer_family'].notna()].copy()

    # Find largest arrays
    array_sizes = classified.groupby(['seq_id', 'array_idx']).size().reset_index(name='size')
    top_arrays = array_sizes.nlargest(n_arrays, 'size')

    if len(top_arrays) == 0:
        print("No arrays found for combined plot")
        return

    fig, axes = plt.subplots(n_arrays, 1, figsize=(16, 3 * n_arrays))
    if n_arrays == 1:
        axes = [axes]

    for idx, (_, arr_info) in enumerate(top_arrays.iterrows()):
        ax = axes[idx]

        seq_id = arr_info['seq_id']
        array_idx = arr_info['array_idx']

        # Get monomers for this array
        arr_monomers = classified[
            (classified['seq_id'] == seq_id) &
            (classified['array_idx'] == array_idx)
        ].sort_values('monomer_idx')

        # Plot each monomer as a vertical bar
        for _, mon in arr_monomers.iterrows():
            family = int(mon['monomer_family'])
            color = FAMILY_COLORS.get(family, '#999999')

            mon_idx = mon['monomer_idx']
            ax.add_patch(mpatches.Rectangle((mon_idx, 0), 1, 1,
                                           facecolor=color, edgecolor='black',
                                           linewidth=0.8))

        # Format axis
        ax.set_xlim(0, len(arr_monomers))
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel('Monomer index', fontsize=11)

        # Title with info
        array_length_kb = (arr_monomers['monomer_end'].max() - arr_monomers['monomer_start'].min()) / 1000
        title = f'Array {idx+1}: {len(arr_monomers)} monomers, {array_length_kb:.1f} kb'
        title += f' | Read: {seq_id[:16]}...'
        ax.set_title(title, fontweight='bold', fontsize=12)

    # Add legend
    legend_elements = []
    families_present = sorted(classified['monomer_family'].unique())
    for fam in families_present:
        fam_int = int(fam)
        color = FAMILY_COLORS.get(fam_int, '#999999')
        legend_elements.append(
            mpatches.Patch(facecolor=color, edgecolor='black',
                         label=f'Family {fam_int}')
        )

    fig.legend(handles=legend_elements, loc='upper right', ncol=2,
              bbox_to_anchor=(0.98, 0.98), framealpha=0.9)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Generated {output_file}")

def plot_individual_array(arr_monomers, output_file, array_name):
    """Plot a single array showing detailed monomer composition"""

    fig, ax = plt.subplots(figsize=(16, 4))

    # Plot each monomer
    for _, mon in arr_monomers.iterrows():
        if pd.notna(mon['monomer_family']):
            family = int(mon['monomer_family'])
            color = FAMILY_COLORS.get(family, '#999999')
        else:
            color = '#CCCCCC'  # Unclassified

        mon_idx = mon['monomer_idx']
        ax.add_patch(mpatches.Rectangle((mon_idx, 0), 1, 1,
                                       facecolor=color, edgecolor='black',
                                       linewidth=1))

    # Format
    ax.set_xlim(0, len(arr_monomers))
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel('Monomer index', fontweight='bold', fontsize=12)

    # Title
    array_length_kb = (arr_monomers['monomer_end'].max() - arr_monomers['monomer_start'].min()) / 1000
    period = arr_monomers['array_period'].iloc[0] if 'array_period' in arr_monomers.columns else 178
    quality = arr_monomers['array_quality'].iloc[0] if 'array_quality' in arr_monomers.columns else 'N/A'

    title = f'{array_name}: {len(arr_monomers)} monomers × {period:.0f}bp = {array_length_kb:.1f} kb'
    title += f' | Quality: {quality}'
    ax.set_title(title, fontweight='bold', fontsize=13)

    # Legend
    legend_elements = []
    families_present = sorted(arr_monomers['monomer_family'].dropna().unique())
    for fam in families_present:
        fam_int = int(fam)
        color = FAMILY_COLORS.get(fam_int, '#999999')
        legend_elements.append(
            mpatches.Patch(facecolor=color, edgecolor='black',
                         label=f'Family {fam_int}')
        )

    if pd.isna(arr_monomers['monomer_family']).any():
        legend_elements.append(
            mpatches.Patch(facecolor='#CCCCCC', edgecolor='black',
                         label='Unclassified')
        )

    ax.legend(handles=legend_elements, loc='upper right', ncol=3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Generated {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive read visualizations')
    parser.add_argument('--classifications', required=True, help='Monomer classifications TSV')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--n-arrays', type=int, default=5, help='Number of individual arrays to plot')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading classifications from {args.classifications}...")
    monomers_df = pd.read_csv(args.classifications, sep='\t')
    print(f"  Loaded {len(monomers_df)} monomers")

    classified = monomers_df[monomers_df['monomer_family'].notna()]
    print(f"  {len(classified)} classified ({len(classified)/len(monomers_df)*100:.1f}%)")

    # 1. Family summary (distribution + transitions)
    print("\n[1/4] Generating family summary with transition heatmap...")
    plot_family_summary(monomers_df, output_dir / 'family_summary.png')

    # 2. Top arrays combined
    print("\n[2/4] Generating top arrays combined plot...")
    plot_top_arrays_combined(monomers_df, output_dir / 'top_arrays_combined.png', n_arrays=3)

    # 3. Individual array plots
    print(f"\n[3/4] Generating individual array plots (top {args.n_arrays})...")

    # Find largest arrays
    array_info = []
    for (seq_id, array_idx), arr_df in monomers_df.groupby(['seq_id', 'array_idx']):
        array_info.append({
            'seq_id': seq_id,
            'array_idx': array_idx,
            'n_monomers': len(arr_df),
            'monomers': arr_df.sort_values('monomer_idx')
        })

    array_info.sort(key=lambda x: x['n_monomers'], reverse=True)

    for i, arr_info in enumerate(array_info[:args.n_arrays]):
        array_name = f"Array {i+1} (idx{arr_info['array_idx']})"
        output_file = output_dir / f'array_{i+1}_idx{arr_info["array_idx"]}.png'
        plot_individual_array(arr_info['monomers'], output_file, array_name)

    # 4. Summary statistics
    print("\n[4/4] Writing summary statistics...")
    with open(output_dir / 'ARRAY_SUMMARY.txt', 'w') as f:
        f.write("="*80 + "\n")
        f.write("CENprofiler Read Analysis - Array Summary\n")
        f.write("="*80 + "\n\n")

        f.write(f"Total monomers: {len(monomers_df):,}\n")
        f.write(f"Classified: {len(classified):,} ({len(classified)/len(monomers_df)*100:.1f}%)\n\n")

        f.write("Family distribution:\n")
        family_counts = classified['monomer_family'].value_counts().sort_index()
        for fam, count in family_counts.items():
            pct = count / len(classified) * 100
            f.write(f"  Family {int(fam):2d}: {count:4d} ({pct:5.1f}%)\n")

        f.write(f"\nTop {min(args.n_arrays, len(array_info))} arrays:\n")
        for i, arr_info in enumerate(array_info[:args.n_arrays]):
            arr_df = arr_info['monomers']
            length_kb = (arr_df['monomer_end'].max() - arr_df['monomer_start'].min()) / 1000
            n_classified = arr_df['monomer_family'].notna().sum()
            f.write(f"  Array {i+1}: {arr_info['n_monomers']:3d} monomers, {length_kb:.1f} kb, ")
            f.write(f"{n_classified} classified ({n_classified/arr_info['n_monomers']*100:.0f}%)\n")

    print("\n✅ All comprehensive visualizations generated!")
    print(f"Output directory: {output_dir}")

if __name__ == '__main__':
    main()
