#!/usr/bin/env python3
"""
Analyze spatial patterns and clustering of monomer families.

This script analyzes:
- Family clustering within arrays (are same-family monomers grouped?)
- Positional preferences (do families prefer array starts/ends?)
- Family co-occurrence patterns
- Spatial autocorrelation
- Boundary analysis (what families occur at array boundaries?)

Usage:
    python analyze_monomer_positions.py <classifications.tsv> <output_dir>
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import defaultdict, Counter
import json

# Family colors
FAMILY_COLORS = {
    1: '#FF0000', 2: '#FFA500', 3: '#FFFF00', 4: '#00FF00', 5: '#00FFFF',
    6: '#0000FF', 7: '#FF00FF', 8: '#8B4513', 9: '#FFC0CB', 10: '#808080',
    11: '#800000', 12: '#FF8C00', 13: '#FFD700', 14: '#32CD32', 15: '#00CED1',
    16: '#4169E1', 17: '#9370DB', 18: '#A0522D', 19: '#FFB6C1', 20: '#696969'
}

def load_data(tsv_file):
    """Load monomer classifications."""
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"Loaded {len(df)} monomers")
    return df

def analyze_positional_preferences(df, min_array_size=20):
    """Analyze if families have positional preferences within arrays."""
    classified = df[df['monomer_family'].notna()].copy()

    # Filter to arrays with sufficient size
    array_sizes = classified.groupby(['seq_id', 'array_idx']).size()
    large_arrays = array_sizes[array_sizes >= min_array_size].index

    classified = classified.set_index(['seq_id', 'array_idx'])
    classified = classified.loc[classified.index.isin(large_arrays)]
    classified = classified.reset_index()

    # Calculate relative position (0=start, 1=end) within each array
    def calc_rel_position(group):
        group = group.sort_values('monomer_idx')
        n = len(group)
        group['rel_position'] = np.arange(n) / (n - 1) if n > 1 else 0.5
        return group

    classified = classified.groupby(['seq_id', 'array_idx']).apply(calc_rel_position).reset_index(drop=True)

    # For each family, calculate distribution of relative positions
    family_positions = {}
    for family in sorted(classified['monomer_family'].unique()):
        family_df = classified[classified['monomer_family'] == family]
        positions = family_df['rel_position'].values

        family_positions[int(family)] = {
            'mean_position': float(np.mean(positions)),
            'std_position': float(np.std(positions)),
            'median_position': float(np.median(positions)),
            'n_occurrences': len(positions)
        }

    return family_positions, classified

def analyze_boundary_enrichment(df):
    """Analyze which families are enriched at array boundaries."""
    classified = df[df['monomer_family'].notna()].copy()

    boundary_counts = defaultdict(lambda: {'start': 0, 'end': 0, 'middle': 0})

    for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
        group = group.sort_values('monomer_idx')

        if len(group) < 3:
            continue

        # First monomer
        start_family = int(group.iloc[0]['monomer_family'])
        boundary_counts[start_family]['start'] += 1

        # Last monomer
        end_family = int(group.iloc[-1]['monomer_family'])
        boundary_counts[end_family]['end'] += 1

        # Middle monomers
        for family in group.iloc[1:-1]['monomer_family']:
            boundary_counts[int(family)]['middle'] += 1

    # Calculate enrichment
    enrichment = []
    for family in sorted(boundary_counts.keys()):
        counts = boundary_counts[family]
        total = counts['start'] + counts['end'] + counts['middle']

        if total < 10:  # Skip rare families
            continue

        # Expected boundary frequency (2 boundaries per array of length n)
        # If uniform, boundary_fraction = 2/n, but average across arrays...
        # Simpler: compare observed vs expected if uniform distribution

        start_pct = counts['start'] / total * 100
        end_pct = counts['end'] / total * 100
        middle_pct = counts['middle'] / total * 100

        # Expected: if uniform, start=end=~1%, middle=~98% for large arrays
        # But this depends on array size, so we compare relative to average

        enrichment.append({
            'family': family,
            'start_count': counts['start'],
            'end_count': counts['end'],
            'middle_count': counts['middle'],
            'total': total,
            'start_pct': start_pct,
            'end_pct': end_pct,
            'middle_pct': middle_pct
        })

    df_boundary = pd.DataFrame(enrichment).sort_values('start_pct', ascending=False)

    return df_boundary

def analyze_clustering(df, window_size=10):
    """Analyze family clustering using local density."""
    classified = df[df['monomer_family'].notna()].copy()

    clustering_scores = {}

    for family in sorted(classified['monomer_family'].unique()):
        scores = []

        for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
            group = group.sort_values('monomer_idx')
            families = group['monomer_family'].values

            if len(families) < window_size:
                continue

            # Sliding window to calculate local density of target family
            for i in range(len(families) - window_size + 1):
                window = families[i:i+window_size]
                local_density = np.sum(window == family) / window_size
                scores.append(local_density)

        if scores:
            # Compare to expected density (overall frequency)
            overall_freq = np.sum(classified['monomer_family'] == family) / len(classified)

            clustering_scores[int(family)] = {
                'mean_local_density': float(np.mean(scores)),
                'overall_frequency': float(overall_freq),
                'clustering_index': float(np.mean(scores) / overall_freq) if overall_freq > 0 else 0,
                'std_local_density': float(np.std(scores))
            }

    return clustering_scores

def analyze_cooccurrence(df, distance_threshold=5):
    """Analyze which families co-occur within a distance threshold."""
    classified = df[df['monomer_family'].notna()].copy()

    cooccur = defaultdict(int)
    family_counts = Counter()

    for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
        group = group.sort_values('monomer_idx')
        families = group['monomer_family'].values

        for i, f1 in enumerate(families):
            f1 = int(f1)
            family_counts[f1] += 1

            # Look ahead within threshold
            for j in range(i+1, min(i+1+distance_threshold, len(families))):
                f2 = int(families[j])

                # Store ordered pair
                pair = tuple(sorted([f1, f2]))
                cooccur[pair] += 1

    # Calculate observed vs expected
    cooccur_analysis = []
    total_comparisons = sum(cooccur.values())

    for (f1, f2), obs_count in cooccur.items():
        # Expected if independent
        freq1 = family_counts[f1] / len(classified)
        freq2 = family_counts[f2] / len(classified)
        expected = total_comparisons * freq1 * freq2

        if expected > 0:
            enrichment = obs_count / expected

            cooccur_analysis.append({
                'family1': f1,
                'family2': f2,
                'observed': obs_count,
                'expected': expected,
                'enrichment': enrichment,
                'log2_enrichment': np.log2(enrichment)
            })

    df_cooccur = pd.DataFrame(cooccur_analysis).sort_values('enrichment', ascending=False)

    return df_cooccur

def plot_positional_preferences(family_positions, classified, output_dir):
    """Plot positional distributions for families."""
    # Filter to families with sufficient data
    families = [f for f, data in family_positions.items() if data['n_occurrences'] > 50]

    if len(families) == 0:
        print("  Not enough data for positional preference plot")
        return

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Distribution plots
    for family in families[:10]:  # Top 10
        family_df = classified[classified['monomer_family'] == family]
        positions = family_df['rel_position'].values

        axes[0].hist(positions, bins=20, alpha=0.5, label=f'F{family}',
                    color=FAMILY_COLORS.get(family, 'gray'))

    axes[0].set_xlabel('Relative Position in Array (0=start, 1=end)')
    axes[0].set_ylabel('Count')
    axes[0].set_title('Family Positional Distributions (Top 10 Families)')
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    axes[0].grid(True, alpha=0.3)

    # Mean position summary
    sorted_families = sorted(families, key=lambda f: family_positions[f]['mean_position'])

    means = [family_positions[f]['mean_position'] for f in sorted_families]
    stds = [family_positions[f]['std_position'] for f in sorted_families]
    colors = [FAMILY_COLORS.get(f, 'gray') for f in sorted_families]

    y_pos = np.arange(len(sorted_families))

    axes[1].barh(y_pos, means, xerr=stds, color=colors, edgecolor='black', alpha=0.7)
    axes[1].set_yticks(y_pos)
    axes[1].set_yticklabels([f'F{f}' for f in sorted_families])
    axes[1].set_xlabel('Mean Relative Position')
    axes[1].set_title('Family Positional Preferences')
    axes[1].axvline(0.5, color='black', linestyle='--', linewidth=1, label='Center')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    plt.savefig(output_dir / 'positional_preferences.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: positional_preferences.png")

def plot_clustering_scores(clustering_scores, output_dir):
    """Plot clustering index for each family."""
    families = sorted(clustering_scores.keys(), key=lambda f: clustering_scores[f]['clustering_index'], reverse=True)

    if len(families) == 0:
        print("  Not enough data for clustering plot")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    indices = [clustering_scores[f]['clustering_index'] for f in families]
    colors = [FAMILY_COLORS.get(f, 'gray') for f in families]

    bars = ax.bar(range(len(families)), indices, color=colors, edgecolor='black')
    ax.set_xticks(range(len(families)))
    ax.set_xticklabels([f'F{f}' for f in families], rotation=45)
    ax.set_ylabel('Clustering Index (Local / Overall Frequency)')
    ax.set_title('Family Clustering Tendency\n(>1 = clustered, <1 = dispersed)')
    ax.axhline(1, color='black', linestyle='--', linewidth=1, label='Random expectation')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(output_dir / 'clustering_scores.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: clustering_scores.png")

def plot_cooccurrence_network(df_cooccur, output_dir, top_n=30):
    """Plot co-occurrence enrichment heatmap."""
    if len(df_cooccur) == 0:
        print("  Not enough data for co-occurrence plot")
        return

    # Take top enrichments
    top = df_cooccur.head(top_n)

    # Get all families involved
    all_families = sorted(set(top['family1']) | set(top['family2']))

    # Create matrix
    n = len(all_families)
    matrix = np.zeros((n, n))
    fam_to_idx = {f: i for i, f in enumerate(all_families)}

    for _, row in top.iterrows():
        i = fam_to_idx[row['family1']]
        j = fam_to_idx[row['family2']]
        val = row['log2_enrichment']
        matrix[i, j] = val
        matrix[j, i] = val

    # Plot
    fig, ax = plt.subplots(figsize=(10, 9))

    im = ax.imshow(matrix, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([f'F{f}' for f in all_families], rotation=45)
    ax.set_yticklabels([f'F{f}' for f in all_families])

    ax.set_title(f'Family Co-occurrence Enrichment\n(log2 observed/expected within {5}bp distance)')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Log2 Enrichment')

    plt.tight_layout()
    plt.savefig(output_dir / 'cooccurrence_heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: cooccurrence_heatmap.png")

def save_position_report(family_positions, df_boundary, clustering_scores,
                        df_cooccur, output_dir):
    """Save comprehensive position analysis report."""
    report_file = output_dir / 'position_analysis.txt'

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("MONOMER POSITIONAL ANALYSIS REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Positional preferences
        f.write("POSITIONAL PREFERENCES (in arrays ≥20 monomers)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Family':<8} {'N':>8} {'Mean Pos':>10} {'Median':>10} {'Std':>10}\n")
        f.write("-" * 80 + "\n")

        for family in sorted(family_positions.keys(), key=lambda f: family_positions[f]['mean_position']):
            data = family_positions[family]
            if data['n_occurrences'] > 50:
                f.write(f"F{family:<7} {data['n_occurrences']:>8} "
                       f"{data['mean_position']:>10.3f} {data['median_position']:>10.3f} "
                       f"{data['std_position']:>10.3f}\n")

        f.write("\nNote: Position 0.0=array start, 0.5=center, 1.0=array end\n\n")

        # Boundary enrichment
        if len(df_boundary) > 0:
            f.write("BOUNDARY ENRICHMENT\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Family':<8} {'Start %':>10} {'End %':>10} {'Middle %':>10} {'Total':>10}\n")
            f.write("-" * 80 + "\n")

            for _, row in df_boundary.head(15).iterrows():
                f.write(f"F{int(row['family']):<7} {row['start_pct']:>10.2f} "
                       f"{row['end_pct']:>10.2f} {row['middle_pct']:>10.2f} "
                       f"{int(row['total']):>10}\n")

            f.write("\n")

        # Clustering
        f.write("CLUSTERING ANALYSIS\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Family':<8} {'Overall %':>10} {'Local %':>10} {'Index':>10} {'Tendency':>12}\n")
        f.write("-" * 80 + "\n")

        for family in sorted(clustering_scores.keys(),
                           key=lambda f: clustering_scores[f]['clustering_index'],
                           reverse=True):
            data = clustering_scores[family]
            index = data['clustering_index']
            tendency = 'Clustered' if index > 1.2 else ('Dispersed' if index < 0.8 else 'Random')

            f.write(f"F{family:<7} {data['overall_frequency']*100:>10.2f} "
                   f"{data['mean_local_density']*100:>10.2f} {index:>10.3f} "
                   f"{tendency:>12}\n")

        f.write("\nNote: Clustering index >1 suggests clustering, <1 suggests dispersion\n\n")

        # Top co-occurrences
        if len(df_cooccur) > 0:
            f.write("TOP 15 CO-OCCURRING FAMILY PAIRS\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Pair':<12} {'Observed':>10} {'Expected':>10} {'Enrichment':>12}\n")
            f.write("-" * 80 + "\n")

            for _, row in df_cooccur.head(15).iterrows():
                pair = f"F{row['family1']}-F{row['family2']}"
                f.write(f"{pair:<12} {int(row['observed']):>10} {row['expected']:>10.1f} "
                       f"{row['enrichment']:>12.2f}\n")

            f.write("\n")

        f.write("=" * 80 + "\n")

    print("  Saved: position_analysis.txt")

    # Save JSON
    json_data = {
        'positional_preferences': family_positions,
        'boundary_enrichment': df_boundary.to_dict(orient='records') if len(df_boundary) > 0 else [],
        'clustering_scores': clustering_scores,
        'cooccurrence': df_cooccur.head(50).to_dict(orient='records') if len(df_cooccur) > 0 else []
    }

    json_file = output_dir / 'position_analysis.json'
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    print("  Saved: position_analysis.json")

def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_monomer_positions.py <classifications.tsv> <output_dir>")
        sys.exit(1)

    tsv_file = sys.argv[1]
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*80)
    print("MONOMER POSITIONAL ANALYSIS")
    print("="*80 + "\n")

    # Load data
    print("Loading data...")
    df = load_data(tsv_file)

    # Positional preferences
    print("\nAnalyzing positional preferences...")
    family_positions, classified = analyze_positional_preferences(df)

    # Boundary enrichment
    print("Analyzing boundary enrichment...")
    df_boundary = analyze_boundary_enrichment(df)

    # Clustering
    print("Analyzing clustering...")
    clustering_scores = analyze_clustering(df)

    # Co-occurrence
    print("Analyzing family co-occurrence...")
    df_cooccur = analyze_cooccurrence(df)

    # Generate plots
    print("\nGenerating plots...")
    plot_positional_preferences(family_positions, classified, output_dir)
    plot_clustering_scores(clustering_scores, output_dir)
    plot_cooccurrence_network(df_cooccur, output_dir)

    # Save reports
    print("\nSaving reports...")
    save_position_report(family_positions, df_boundary, clustering_scores,
                        df_cooccur, output_dir)

    # Save tables
    if len(df_boundary) > 0:
        df_boundary.to_csv(output_dir / 'boundary_enrichment.tsv', sep='\t', index=False)
        print("  Saved: boundary_enrichment.tsv")

    if len(df_cooccur) > 0:
        df_cooccur.to_csv(output_dir / 'cooccurrence.tsv', sep='\t', index=False)
        print("  Saved: cooccurrence.tsv")

    print("\n" + "="*80)
    print("POSITIONAL ANALYSIS COMPLETE")
    print(f"Results saved to: {output_dir}")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
