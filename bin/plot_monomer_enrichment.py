#!/usr/bin/env python3
"""
Monomer family HOR enrichment analysis and visualization

This script analyzes which monomer families are preferentially found in HORs
and generates a comprehensive 6-panel visualization showing enrichment patterns.

Usage:
    python plot_monomer_enrichment.py <hors_file> <classifications_file> <output_prefix>

Arguments:
    hors_file: TSV file with HOR detections (must have: seq_id, hor_start, hor_end, hor_type)
    classifications_file: TSV file with monomer classifications (must have: seq_id, monomer_start, monomer_end, monomer_family)
    output_prefix: Prefix for output files (will create <prefix>.png and <prefix>.tsv)
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
from pathlib import Path

def calculate_enrichment(hors_file, classifications_file):
    """Calculate monomer family HOR enrichment statistics"""

    print("Loading data...")
    hors = pd.read_csv(hors_file, sep='\t')
    classified = pd.read_csv(classifications_file, sep='\t')

    print(f"Loaded {len(hors):,} HORs and {len(classified):,} monomers")

    # Mark which monomers are in HORs
    classified['in_hor'] = False

    print("Marking monomers in HORs...")
    for _, hor in hors.iterrows():
        mask = (classified['seq_id'] == hor['seq_id']) & \
               (classified['monomer_start'] >= hor['hor_start']) & \
               (classified['monomer_end'] <= hor['hor_end'])
        classified.loc[mask, 'in_hor'] = True

    # Family-level statistics
    family_counts = classified['monomer_family'].value_counts().sort_index()
    family_in_hor = classified[classified['in_hor']]['monomer_family'].value_counts().sort_index()

    # Build enrichment dataframe
    enrichment_data = []
    for fam in family_counts.index:
        total_count = family_counts[fam]
        in_hor_count = family_in_hor.get(fam, 0)
        enrichment_pct = (in_hor_count / total_count) * 100

        # Chromosome distribution
        fam_monomers = classified[classified['monomer_family'] == fam]
        chr_dist = fam_monomers['seq_id'].str.extract(r'(Chr\d+)')[0].value_counts()

        # Get unique HORs this family participates in
        fam_in_hors = classified[(classified['monomer_family'] == fam) & (classified['in_hor'])]

        unique_hors_set = set()
        for _, mono in fam_in_hors.iterrows():
            matching_hors = hors[
                (hors['seq_id'] == mono['seq_id']) &
                (hors['hor_start'] <= mono['monomer_start']) &
                (hors['hor_end'] >= mono['monomer_end'])
            ]
            for _, hor_row in matching_hors.iterrows():
                unique_hors_set.add((hor_row['seq_id'], hor_row['hor_start'], hor_row['hor_end']))

        num_unique_hors = len(unique_hors_set)

        enrichment_data.append({
            'family': int(fam),
            'total_monomers': total_count,
            'in_hor': in_hor_count,
            'enrichment_pct': enrichment_pct,
            'unique_hors': num_unique_hors,
            'dominant_chr': chr_dist.index[0] if len(chr_dist) > 0 else 'NA'
        })

    enrichment_df = pd.DataFrame(enrichment_data).sort_values('enrichment_pct', ascending=False)

    return enrichment_df, hors, classified


def plot_enrichment(enrichment_df, hors, classified, output_prefix):
    """Create 6-panel comprehensive enrichment visualization"""

    fig, axes = plt.subplots(3, 2, figsize=(16, 14))

    # Plot 1: Enrichment by family (Top 20)
    ax1 = axes[0, 0]
    top20 = enrichment_df.head(20)
    colors = plt.cm.viridis(np.linspace(0, 1, len(top20)))
    ax1.barh(range(len(top20)), top20['enrichment_pct'], color=colors)
    ax1.set_yticks(range(len(top20)))
    ax1.set_yticklabels([f'F{int(f)}' for f in top20['family']])
    ax1.set_xlabel('% of Family Monomers in HORs', fontweight='bold')
    ax1.set_ylabel('Family', fontweight='bold')
    ax1.set_title('HOR Enrichment by Family (Top 20)', fontweight='bold')
    ax1.invert_yaxis()

    # Plot 2: Scatter - Family size vs enrichment
    ax2 = axes[0, 1]
    scatter = ax2.scatter(enrichment_df['total_monomers'],
                         enrichment_df['enrichment_pct'],
                         c=enrichment_df['family'], cmap='tab20',
                         s=100, alpha=0.7, edgecolors='black', linewidth=1)
    ax2.set_xlabel('Total Monomers', fontweight='bold')
    ax2.set_ylabel('% HORs', fontweight='bold')
    ax2.set_title('Family Size vs HOR Enrichment', fontweight='bold')
    ax2.set_xscale('log')

    # Annotate top families
    for _, row in enrichment_df.head(5).iterrows():
        ax2.annotate(f'F{int(row["family"])}',
                    (row['total_monomers'], row['enrichment_pct']),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, fontweight='bold')

    # Plot 3: Chromosome distribution heatmap
    ax3 = axes[1, 0]
    chr_family_matrix = []
    families_to_plot = enrichment_df.head(15)['family'].values
    chromosomes = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']

    for fam in families_to_plot:
        row = []
        fam_monomers = classified[classified['monomer_family'] == fam]
        total = len(fam_monomers)
        for chrom in chromosomes:
            count = len(fam_monomers[fam_monomers['seq_id'].str.contains(chrom)])
            row.append((count / total * 100) if total > 0 else 0)
        chr_family_matrix.append(row)

    sns.heatmap(chr_family_matrix, ax=ax3, cmap='YlOrRd',
               xticklabels=chromosomes,
               yticklabels=[f'F{int(f)}' for f in families_to_plot],
               cbar_kws={'label': '% of Family'},
               annot=True, fmt='.1f', annot_kws={'fontsize': 8})
    ax3.set_title('Chromosome Distribution (Top 15 Families)', fontweight='bold')
    ax3.set_ylabel('Family', fontweight='bold')

    # Plot 4: Unique HORs per family
    ax4 = axes[1, 1]
    top15_hors = enrichment_df.head(15)
    ax4.barh(range(len(top15_hors)), top15_hors['unique_hors'],
            color='steelblue', edgecolor='black')
    ax4.set_yticks(range(len(top15_hors)))
    ax4.set_yticklabels([f'F{int(f)}' for f in top15_hors['family']])
    ax4.set_xlabel('Number of Unique HORs', fontweight='bold')
    ax4.set_ylabel('Family', fontweight='bold')
    ax4.set_title('Unique HORs per Family (Top 15)', fontweight='bold')
    ax4.invert_yaxis()

    # Plot 5: Monomers in HORs vs not in HORs
    ax5 = axes[2, 0]
    top10 = enrichment_df.head(10)
    x_pos = np.arange(len(top10))
    width = 0.35

    ax5.bar(x_pos - width/2, top10['in_hor'], width,
           label='In HORs', color='#2ecc71', edgecolor='black')
    ax5.bar(x_pos + width/2, top10['total_monomers'] - top10['in_hor'], width,
           label='Not in HORs', color='#e74c3c', edgecolor='black')

    ax5.set_xlabel('Family', fontweight='bold')
    ax5.set_ylabel('Monomer Count', fontweight='bold')
    ax5.set_title('Monomers in/out of HORs (Top 10)', fontweight='bold')
    ax5.set_xticks(x_pos)
    ax5.set_xticklabels([f'F{int(f)}' for f in top10['family']])
    ax5.legend()
    ax5.set_yscale('log')

    # Plot 6: Summary statistics
    ax6 = axes[2, 1]
    ax6.axis('off')

    n_homhor = len(hors[hors['hor_type'] == 'homHOR']) if 'hor_type' in hors.columns else 0
    n_hethor = len(hors[hors['hor_type'] == 'hetHOR']) if 'hor_type' in hors.columns else 0

    summary_text = f"""
HOR ENRICHMENT SUMMARY

Total Monomers: {len(classified):,}
In HORs: {classified['in_hor'].sum():,} ({classified['in_hor'].sum() / len(classified) * 100:.1f}%)

Total HORs: {len(hors):,}
  homHORs: {n_homhor:,} ({n_homhor/len(hors)*100:.1f}%)
  hetHORs: {n_hethor:,} ({n_hethor/len(hors)*100:.1f}%)

Top Enriched Families:
  F{int(enrichment_df.iloc[0]['family'])}: {enrichment_df.iloc[0]['enrichment_pct']:.1f}%
  F{int(enrichment_df.iloc[1]['family'])}: {enrichment_df.iloc[1]['enrichment_pct']:.1f}%
  F{int(enrichment_df.iloc[2]['family'])}: {enrichment_df.iloc[2]['enrichment_pct']:.1f}%
"""

    ax6.text(0.1, 0.95, summary_text, transform=ax6.transAxes,
            fontsize=11, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Saved: {output_prefix}.png")


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    hors_file = sys.argv[1]
    classifications_file = sys.argv[2]
    output_prefix = sys.argv[3]

    # Validate input files
    if not Path(hors_file).exists():
        print(f"ERROR: HORs file not found: {hors_file}")
        sys.exit(1)

    if not Path(classifications_file).exists():
        print(f"ERROR: Classifications file not found: {classifications_file}")
        sys.exit(1)

    print("=" * 70)
    print("MONOMER FAMILY HOR ENRICHMENT ANALYSIS")
    print("=" * 70)

    # Calculate enrichment
    enrichment_df, hors, classified = calculate_enrichment(hors_file, classifications_file)

    # Save enrichment table
    enrichment_df.to_csv(f'{output_prefix}.tsv', sep='\t', index=False)
    print(f"Saved: {output_prefix}.tsv")

    # Print top enriched families
    print("\nTop 10 families by HOR enrichment:")
    print(enrichment_df.head(10).to_string(index=False))

    # Create visualization
    print("\nGenerating visualization...")
    plot_enrichment(enrichment_df, hors, classified, output_prefix)

    # Additional statistics
    print(f"\n=== ADDITIONAL STATISTICS ===")
    print(f"Families with 0% enrichment: {len(enrichment_df[enrichment_df['enrichment_pct'] == 0])}")
    print(f"Families with >10% enrichment: {len(enrichment_df[enrichment_df['enrichment_pct'] > 10])}")
    print(f"Families with >50% enrichment: {len(enrichment_df[enrichment_df['enrichment_pct'] > 50])}")
    print(f"Families with >90% enrichment: {len(enrichment_df[enrichment_df['enrichment_pct'] > 90])}")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
