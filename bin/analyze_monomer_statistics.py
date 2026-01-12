#!/usr/bin/env python3
"""
Comprehensive monomer-level statistics analysis.

This script analyzes monomer classifications to provide:
- Family composition statistics
- Quality metrics (length distribution, identity)
- Transition patterns between families
- Spatial clustering analysis
- Array-level heterogeneity metrics

Usage:
    python analyze_monomer_statistics.py <monomer_classifications.tsv> <output_dir>
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter, defaultdict
import json

# Family colors matching the visualization pipeline
FAMILY_COLORS = {
    1: '#FF0000', 2: '#FFA500', 3: '#FFFF00', 4: '#00FF00', 5: '#00FFFF',
    6: '#0000FF', 7: '#FF00FF', 8: '#8B4513', 9: '#FFC0CB', 10: '#808080',
    11: '#800000', 12: '#FF8C00', 13: '#FFD700', 14: '#32CD32', 15: '#00CED1',
    16: '#4169E1', 17: '#9370DB', 18: '#A0522D', 19: '#FFB6C1', 20: '#696969'
}

def load_monomer_data(tsv_file):
    """Load monomer classifications."""
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"Loaded {len(df)} monomers from {tsv_file}")
    return df

def calculate_basic_statistics(df):
    """Calculate basic monomer statistics."""
    stats = {}

    # Total counts
    stats['total_monomers'] = len(df)
    stats['classified_monomers'] = df['monomer_family'].notna().sum()
    stats['unclassified_monomers'] = df['monomer_family'].isna().sum()
    stats['classification_rate'] = stats['classified_monomers'] / stats['total_monomers'] * 100

    # Length statistics
    stats['mean_length'] = df['monomer_length'].mean()
    stats['median_length'] = df['monomer_length'].median()
    stats['std_length'] = df['monomer_length'].std()
    stats['min_length'] = df['monomer_length'].min()
    stats['max_length'] = df['monomer_length'].max()

    # Identity statistics (for classified monomers)
    classified = df[df['monomer_family'].notna()]
    if len(classified) > 0:
        stats['mean_identity'] = classified['alignment_identity'].mean()
        stats['median_identity'] = classified['alignment_identity'].median()
        stats['std_identity'] = classified['alignment_identity'].std()
        stats['min_identity'] = classified['alignment_identity'].min()
        stats['max_identity'] = classified['alignment_identity'].max()

    # Array statistics
    stats['total_arrays'] = df['array_idx'].nunique()
    stats['mean_monomers_per_array'] = df.groupby(['seq_id', 'array_idx']).size().mean()
    stats['median_monomers_per_array'] = df.groupby(['seq_id', 'array_idx']).size().median()

    # Read statistics
    stats['total_reads'] = df['seq_id'].nunique()
    stats['mean_arrays_per_read'] = stats['total_arrays'] / stats['total_reads']

    return stats

def calculate_family_statistics(df):
    """Calculate per-family statistics."""
    classified = df[df['monomer_family'].notna()].copy()

    family_stats = []
    for family in sorted(classified['monomer_family'].unique()):
        family_df = classified[classified['monomer_family'] == family]

        stats = {
            'family': int(family),
            'count': len(family_df),
            'percentage': len(family_df) / len(classified) * 100,
            'mean_length': family_df['monomer_length'].mean(),
            'std_length': family_df['monomer_length'].std(),
            'mean_identity': family_df['alignment_identity'].mean(),
            'std_identity': family_df['alignment_identity'].std(),
            'arrays_present': family_df.groupby(['seq_id', 'array_idx']).ngroups,
            'reads_present': family_df['seq_id'].nunique()
        }
        family_stats.append(stats)

    return pd.DataFrame(family_stats).sort_values('count', ascending=False)

def calculate_transition_matrix(df):
    """Calculate family transition matrix."""
    classified = df[df['monomer_family'].notna()].copy()

    # Group by array
    transitions = defaultdict(int)

    for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
        group = group.sort_values('monomer_idx')
        families = group['monomer_family'].values

        for i in range(len(families) - 1):
            f1, f2 = int(families[i]), int(families[i+1])
            transitions[(f1, f2)] += 1

    # Convert to matrix
    all_families = sorted(classified['monomer_family'].unique())
    n_families = len(all_families)
    matrix = np.zeros((n_families, n_families))

    family_to_idx = {int(f): i for i, f in enumerate(all_families)}

    for (f1, f2), count in transitions.items():
        i, j = family_to_idx[f1], family_to_idx[f2]
        matrix[i, j] = count

    return matrix, [int(f) for f in all_families], transitions

def calculate_array_heterogeneity(df):
    """Calculate heterogeneity metrics for each array."""
    classified = df[df['monomer_family'].notna()].copy()

    heterogeneity = []
    for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
        families = group['monomer_family'].values

        metrics = {
            'seq_id': seq_id,
            'array_idx': array_idx,
            'n_monomers': len(families),
            'n_families': len(np.unique(families)),
            'dominant_family': int(Counter(families).most_common(1)[0][0]),
            'dominant_fraction': Counter(families).most_common(1)[0][1] / len(families),
            'shannon_entropy': -sum((count/len(families)) * np.log(count/len(families))
                                   for count in Counter(families).values()),
            'simpson_diversity': 1 - sum((count/len(families))**2
                                        for count in Counter(families).values())
        }

        # Calculate run lengths (consecutive same-family stretches)
        runs = []
        current_run = 1
        for i in range(1, len(families)):
            if families[i] == families[i-1]:
                current_run += 1
            else:
                runs.append(current_run)
                current_run = 1
        runs.append(current_run)

        metrics['mean_run_length'] = np.mean(runs)
        metrics['max_run_length'] = max(runs)
        metrics['n_transitions'] = len(runs) - 1

        heterogeneity.append(metrics)

    return pd.DataFrame(heterogeneity)

def plot_length_distribution(df, output_dir):
    """Plot monomer length distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Overall length distribution
    axes[0].hist(df['monomer_length'], bins=50, edgecolor='black', alpha=0.7)
    axes[0].axvline(178, color='red', linestyle='--', linewidth=2, label='Expected 178bp')
    axes[0].set_xlabel('Monomer Length (bp)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Monomer Length Distribution')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    # Length distribution by family (top 10)
    classified = df[df['monomer_family'].notna()]
    top_families = classified['monomer_family'].value_counts().head(10).index

    for family in top_families:
        family_df = classified[classified['monomer_family'] == family]
        axes[1].hist(family_df['monomer_length'], bins=30, alpha=0.5,
                    label=f'F{int(family)}', color=FAMILY_COLORS.get(int(family), 'gray'))

    axes[1].axvline(178, color='black', linestyle='--', linewidth=2, label='Expected')
    axes[1].set_xlabel('Monomer Length (bp)')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Length Distribution by Family (Top 10)')
    axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'length_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: length_distribution.png")

def plot_identity_distribution(df, output_dir):
    """Plot alignment identity distribution."""
    classified = df[df['monomer_family'].notna()]

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Overall identity distribution
    axes[0].hist(classified['alignment_identity'], bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Alignment Identity (%)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Monomer Alignment Identity Distribution')
    axes[0].grid(True, alpha=0.3)

    # Identity by family (top 10)
    top_families = classified['monomer_family'].value_counts().head(10).index

    family_identities = []
    family_labels = []
    for family in sorted(top_families):
        family_df = classified[classified['monomer_family'] == family]
        family_identities.append(family_df['alignment_identity'].values)
        family_labels.append(f'F{int(family)}')

    bp = axes[1].boxplot(family_identities, labels=family_labels, patch_artist=True)
    for patch, family in zip(bp['boxes'], sorted(top_families)):
        patch.set_facecolor(FAMILY_COLORS.get(int(family), 'gray'))

    axes[1].set_xlabel('Family')
    axes[1].set_ylabel('Alignment Identity (%)')
    axes[1].set_title('Identity Distribution by Family (Top 10)')
    axes[1].grid(True, alpha=0.3, axis='y')
    axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(output_dir / 'identity_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: identity_distribution.png")

def plot_family_composition(family_stats, output_dir):
    """Plot family composition bar chart."""
    fig, ax = plt.subplots(figsize=(12, 6))

    families = family_stats['family'].values
    counts = family_stats['count'].values
    colors = [FAMILY_COLORS.get(int(f), 'gray') for f in families]

    bars = ax.bar(range(len(families)), counts, color=colors, edgecolor='black', linewidth=1)
    ax.set_xticks(range(len(families)))
    ax.set_xticklabels([f'F{int(f)}' for f in families], rotation=45)
    ax.set_xlabel('Family')
    ax.set_ylabel('Number of Monomers')
    ax.set_title('Monomer Family Composition')
    ax.grid(True, alpha=0.3, axis='y')

    # Add percentage labels
    for i, (count, pct) in enumerate(zip(counts, family_stats['percentage'].values)):
        ax.text(i, count, f'{pct:.1f}%', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'family_composition.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: family_composition.png")

def plot_transition_matrix(matrix, families, output_dir):
    """Plot family transition heatmap."""
    fig, ax = plt.subplots(figsize=(12, 10))

    # Normalize by row (show probability of transition)
    row_sums = matrix.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    matrix_norm = matrix / row_sums * 100

    im = ax.imshow(matrix_norm, cmap='YlOrRd', aspect='auto')

    # Set ticks
    ax.set_xticks(range(len(families)))
    ax.set_yticks(range(len(families)))
    ax.set_xticklabels([f'F{f}' for f in families])
    ax.set_yticklabels([f'F{f}' for f in families])

    ax.set_xlabel('To Family')
    ax.set_ylabel('From Family')
    ax.set_title('Family Transition Probabilities (%)')

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Probability (%)')

    # Add text annotations for significant transitions (>5%)
    for i in range(len(families)):
        for j in range(len(families)):
            if matrix_norm[i, j] > 5:
                text = ax.text(j, i, f'{matrix_norm[i, j]:.1f}',
                             ha="center", va="center", color="black", fontsize=8)

    plt.tight_layout()
    plt.savefig(output_dir / 'transition_matrix.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: transition_matrix.png")

def plot_heterogeneity_metrics(heterogeneity, output_dir):
    """Plot array heterogeneity metrics."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Number of families per array
    axes[0, 0].hist(heterogeneity['n_families'], bins=range(1, heterogeneity['n_families'].max()+2),
                    edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Number of Families in Array')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Family Diversity per Array')
    axes[0, 0].grid(True, alpha=0.3)

    # Dominant family fraction
    axes[0, 1].hist(heterogeneity['dominant_fraction'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 1].set_xlabel('Dominant Family Fraction')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Dominance of Most Abundant Family')
    axes[0, 1].grid(True, alpha=0.3)

    # Shannon entropy
    axes[1, 0].hist(heterogeneity['shannon_entropy'], bins=30, edgecolor='black', alpha=0.7)
    axes[1, 0].set_xlabel('Shannon Entropy')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('Array Complexity (Shannon Entropy)')
    axes[1, 0].grid(True, alpha=0.3)

    # Mean run length
    axes[1, 1].hist(heterogeneity['mean_run_length'], bins=30, edgecolor='black', alpha=0.7)
    axes[1, 1].set_xlabel('Mean Run Length (monomers)')
    axes[1, 1].set_ylabel('Count')
    axes[1, 1].set_title('Average Consecutive Same-Family Stretch')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'heterogeneity_metrics.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: heterogeneity_metrics.png")

def plot_array_size_vs_diversity(heterogeneity, output_dir):
    """Plot relationship between array size and diversity."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Array size vs number of families
    axes[0].scatter(heterogeneity['n_monomers'], heterogeneity['n_families'],
                   alpha=0.5, s=20)
    axes[0].set_xlabel('Array Size (monomers)')
    axes[0].set_ylabel('Number of Families')
    axes[0].set_title('Array Size vs Family Diversity')
    axes[0].grid(True, alpha=0.3)

    # Array size vs Shannon entropy
    axes[1].scatter(heterogeneity['n_monomers'], heterogeneity['shannon_entropy'],
                   alpha=0.5, s=20)
    axes[1].set_xlabel('Array Size (monomers)')
    axes[1].set_ylabel('Shannon Entropy')
    axes[1].set_title('Array Size vs Complexity')
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'array_size_vs_diversity.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: array_size_vs_diversity.png")

def save_statistics_report(basic_stats, family_stats, heterogeneity, transitions, output_dir):
    """Save comprehensive statistics report."""
    report_file = output_dir / 'monomer_statistics.txt'

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MONOMER-LEVEL STATISTICS REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Basic statistics
        f.write("BASIC STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total monomers:           {basic_stats['total_monomers']:,}\n")
        f.write(f"Classified monomers:      {basic_stats['classified_monomers']:,} ({basic_stats['classification_rate']:.1f}%)\n")
        f.write(f"Unclassified monomers:    {basic_stats['unclassified_monomers']:,}\n")
        f.write(f"Total arrays:             {basic_stats['total_arrays']:,}\n")
        f.write(f"Total reads:              {basic_stats['total_reads']:,}\n")
        f.write(f"Mean arrays per read:     {basic_stats['mean_arrays_per_read']:.2f}\n")
        f.write(f"Mean monomers per array:  {basic_stats['mean_monomers_per_array']:.1f}\n")
        f.write(f"Median monomers per array: {basic_stats['median_monomers_per_array']:.1f}\n\n")

        # Length statistics
        f.write("LENGTH STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Mean length:    {basic_stats['mean_length']:.1f} bp\n")
        f.write(f"Median length:  {basic_stats['median_length']:.1f} bp\n")
        f.write(f"Std dev:        {basic_stats['std_length']:.1f} bp\n")
        f.write(f"Min length:     {basic_stats['min_length']:.0f} bp\n")
        f.write(f"Max length:     {basic_stats['max_length']:.0f} bp\n\n")

        # Identity statistics
        if 'mean_identity' in basic_stats:
            f.write("ALIGNMENT IDENTITY STATISTICS\n")
            f.write("-" * 80 + "\n")
            f.write(f"Mean identity:    {basic_stats['mean_identity']:.1f}%\n")
            f.write(f"Median identity:  {basic_stats['median_identity']:.1f}%\n")
            f.write(f"Std dev:          {basic_stats['std_identity']:.1f}%\n")
            f.write(f"Min identity:     {basic_stats['min_identity']:.1f}%\n")
            f.write(f"Max identity:     {basic_stats['max_identity']:.1f}%\n\n")

        # Family statistics
        f.write("FAMILY COMPOSITION\n")
        f.write("-" * 80 + "\n")
        f.write(f"Number of families detected: {len(family_stats)}\n\n")
        f.write(f"{'Family':<8} {'Count':>10} {'%':>8} {'Mean Len':>10} {'Mean ID':>10} {'Arrays':>8} {'Reads':>8}\n")
        f.write("-" * 80 + "\n")
        for _, row in family_stats.iterrows():
            f.write(f"F{int(row['family']):<7} {int(row['count']):>10} {row['percentage']:>7.1f} "
                   f"{row['mean_length']:>10.1f} {row['mean_identity']:>10.1f} "
                   f"{int(row['arrays_present']):>8} {int(row['reads_present']):>8}\n")

        f.write("\n")

        # Top transitions
        f.write("TOP 20 FAMILY TRANSITIONS\n")
        f.write("-" * 80 + "\n")
        sorted_transitions = sorted(transitions.items(), key=lambda x: x[1], reverse=True)[:20]
        f.write(f"{'Transition':<15} {'Count':>10}\n")
        f.write("-" * 80 + "\n")
        for (f1, f2), count in sorted_transitions:
            f.write(f"F{f1} → F{f2:<7} {count:>10}\n")

        f.write("\n")

        # Heterogeneity statistics
        f.write("ARRAY HETEROGENEITY STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Mean families per array:        {heterogeneity['n_families'].mean():.2f}\n")
        f.write(f"Median families per array:      {heterogeneity['n_families'].median():.0f}\n")
        f.write(f"Mean dominant family fraction:  {heterogeneity['dominant_fraction'].mean():.2%}\n")
        f.write(f"Mean Shannon entropy:           {heterogeneity['shannon_entropy'].mean():.3f}\n")
        f.write(f"Mean Simpson diversity:         {heterogeneity['simpson_diversity'].mean():.3f}\n")
        f.write(f"Mean run length:                {heterogeneity['mean_run_length'].mean():.2f} monomers\n")
        f.write(f"Mean transitions per array:     {heterogeneity['n_transitions'].mean():.1f}\n")

        f.write("\n")
        f.write("=" * 80 + "\n")

    print(f"  Saved: monomer_statistics.txt")

    # Save JSON for programmatic access
    json_data = {
        'basic_statistics': {k: float(v) if isinstance(v, (np.integer, np.floating)) else v
                           for k, v in basic_stats.items()},
        'family_statistics': family_stats.to_dict(orient='records'),
        'heterogeneity_summary': {
            'mean_families_per_array': float(heterogeneity['n_families'].mean()),
            'mean_dominant_fraction': float(heterogeneity['dominant_fraction'].mean()),
            'mean_shannon_entropy': float(heterogeneity['shannon_entropy'].mean()),
            'mean_run_length': float(heterogeneity['mean_run_length'].mean())
        }
    }

    json_file = output_dir / 'monomer_statistics.json'
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    print(f"  Saved: monomer_statistics.json")

def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_monomer_statistics.py <monomer_classifications.tsv> <output_dir>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*80)
    print("MONOMER STATISTICS ANALYSIS")
    print("="*80 + "\n")

    # Load data
    print("Loading data...")
    df = load_monomer_data(tsv_file)

    # Calculate statistics
    print("\nCalculating basic statistics...")
    basic_stats = calculate_basic_statistics(df)

    print("Calculating family statistics...")
    family_stats = calculate_family_statistics(df)

    print("Calculating transition matrix...")
    matrix, families, transitions = calculate_transition_matrix(df)

    print("Calculating heterogeneity metrics...")
    heterogeneity = calculate_array_heterogeneity(df)

    # Generate plots
    print("\nGenerating plots...")
    plot_length_distribution(df, output_dir)
    plot_identity_distribution(df, output_dir)
    plot_family_composition(family_stats, output_dir)
    plot_transition_matrix(matrix, families, output_dir)
    plot_heterogeneity_metrics(heterogeneity, output_dir)
    plot_array_size_vs_diversity(heterogeneity, output_dir)

    # Save reports
    print("\nSaving statistics reports...")
    save_statistics_report(basic_stats, family_stats, heterogeneity, transitions, output_dir)

    # Save detailed tables
    family_stats.to_csv(output_dir / 'family_statistics.tsv', sep='\t', index=False)
    print(f"  Saved: family_statistics.tsv")

    heterogeneity.to_csv(output_dir / 'array_heterogeneity.tsv', sep='\t', index=False)
    print(f"  Saved: array_heterogeneity.tsv")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print(f"Results saved to: {output_dir}")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
