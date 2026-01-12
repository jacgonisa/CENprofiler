#!/usr/bin/env python3
"""
Compare monomer-level statistics between two samples.

This script performs comprehensive comparative analysis:
- Family composition differences
- Statistical significance testing (chi-square, t-tests)
- Transition pattern comparisons
- Heterogeneity metric comparisons
- Side-by-side visualizations

Usage:
    python compare_samples.py <sample1_classifications.tsv> <sample2_classifications.tsv> \\
                               <sample1_name> <sample2_name> <output_dir>
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from collections import Counter, defaultdict
from scipy import stats
import json

# Family colors matching the visualization pipeline
FAMILY_COLORS = {
    1: '#FF0000', 2: '#FFA500', 3: '#FFFF00', 4: '#00FF00', 5: '#00FFFF',
    6: '#0000FF', 7: '#FF00FF', 8: '#8B4513', 9: '#FFC0CB', 10: '#808080',
    11: '#800000', 12: '#FF8C00', 13: '#FFD700', 14: '#32CD32', 15: '#00CED1',
    16: '#4169E1', 17: '#9370DB', 18: '#A0522D', 19: '#FFB6C1', 20: '#696969'
}

def load_sample_data(tsv_file, sample_name):
    """Load monomer classifications for a sample."""
    df = pd.read_csv(tsv_file, sep='\t')
    df['sample'] = sample_name
    print(f"{sample_name}: Loaded {len(df)} monomers, {df['monomer_family'].notna().sum()} classified")
    return df

def compare_family_composition(df1, df2, sample1, sample2):
    """Compare family composition between samples."""
    # Get classified monomers
    c1 = df1[df1['monomer_family'].notna()]
    c2 = df2[df2['monomer_family'].notna()]

    # Count families
    counts1 = c1['monomer_family'].value_counts()
    counts2 = c2['monomer_family'].value_counts()

    # Get all families present in either sample
    all_families = sorted(set(counts1.index) | set(counts2.index))

    # Build comparison table
    comparison = []
    for family in all_families:
        count1 = counts1.get(family, 0)
        count2 = counts2.get(family, 0)
        total1 = len(c1)
        total2 = len(c2)

        pct1 = (count1 / total1 * 100) if total1 > 0 else 0
        pct2 = (count2 / total2 * 100) if total2 > 0 else 0
        diff = pct2 - pct1

        comparison.append({
            'family': int(family),
            f'{sample1}_count': int(count1),
            f'{sample1}_pct': pct1,
            f'{sample2}_count': int(count2),
            f'{sample2}_pct': pct2,
            'pct_diff': diff,
            'fold_change': (pct2 / pct1) if pct1 > 0 else float('inf')
        })

    df_comp = pd.DataFrame(comparison).sort_values('pct_diff', key=abs, ascending=False)

    # Chi-square test for overall composition difference
    # Build contingency table
    observed = []
    for family in all_families:
        observed.append([counts1.get(family, 0), counts2.get(family, 0)])

    observed = np.array(observed)
    chi2, p_value, dof, expected = stats.chi2_contingency(observed)

    return df_comp, chi2, p_value

def compare_transitions(df1, df2, sample1, sample2):
    """Compare transition patterns between samples."""
    def get_transitions(df):
        transitions = defaultdict(int)
        classified = df[df['monomer_family'].notna()]

        for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
            group = group.sort_values('monomer_idx')
            families = group['monomer_family'].values

            for i in range(len(families) - 1):
                f1, f2 = int(families[i]), int(families[i+1])
                transitions[(f1, f2)] += 1

        return transitions

    trans1 = get_transitions(df1)
    trans2 = get_transitions(df2)

    # Get all transition pairs
    all_pairs = sorted(set(trans1.keys()) | set(trans2.keys()))

    comparison = []
    for (f1, f2) in all_pairs:
        count1 = trans1.get((f1, f2), 0)
        count2 = trans2.get((f1, f2), 0)

        # Normalize by total transitions
        total1 = sum(trans1.values())
        total2 = sum(trans2.values())

        pct1 = (count1 / total1 * 100) if total1 > 0 else 0
        pct2 = (count2 / total2 * 100) if total2 > 0 else 0

        comparison.append({
            'from_family': f1,
            'to_family': f2,
            f'{sample1}_count': count1,
            f'{sample1}_pct': pct1,
            f'{sample2}_count': count2,
            f'{sample2}_pct': pct2,
            'pct_diff': pct2 - pct1
        })

    df_trans = pd.DataFrame(comparison).sort_values('pct_diff', key=abs, ascending=False)

    return df_trans

def compare_heterogeneity(df1, df2, sample1, sample2):
    """Compare array heterogeneity metrics."""
    def calc_heterogeneity(df):
        classified = df[df['monomer_family'].notna()]
        metrics = []

        for (seq_id, array_idx), group in classified.groupby(['seq_id', 'array_idx']):
            families = group['monomer_family'].values

            n_families = len(np.unique(families))
            shannon = -sum((count/len(families)) * np.log(count/len(families))
                          for count in Counter(families).values())
            simpson = 1 - sum((count/len(families))**2
                            for count in Counter(families).values())

            metrics.append({
                'n_monomers': len(families),
                'n_families': n_families,
                'shannon_entropy': shannon,
                'simpson_diversity': simpson
            })

        return pd.DataFrame(metrics)

    het1 = calc_heterogeneity(df1)
    het2 = calc_heterogeneity(df2)

    # Statistical comparisons
    comparisons = {}

    for metric in ['n_families', 'shannon_entropy', 'simpson_diversity']:
        vals1 = het1[metric].values
        vals2 = het2[metric].values

        # t-test
        t_stat, p_val = stats.ttest_ind(vals1, vals2)

        comparisons[metric] = {
            f'{sample1}_mean': float(np.mean(vals1)),
            f'{sample1}_std': float(np.std(vals1)),
            f'{sample2}_mean': float(np.mean(vals2)),
            f'{sample2}_std': float(np.std(vals2)),
            't_statistic': float(t_stat),
            'p_value': float(p_val),
            'significant': p_val < 0.05
        }

    return comparisons, het1, het2

def plot_family_comparison(df_comp, sample1, sample2, output_dir):
    """Plot family composition comparison."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Bar plot comparison
    families = df_comp['family'].values
    x = np.arange(len(families))
    width = 0.35

    pct1 = df_comp[f'{sample1}_pct'].values
    pct2 = df_comp[f'{sample2}_pct'].values

    colors1 = [FAMILY_COLORS.get(int(f), 'gray') for f in families]
    colors2 = [FAMILY_COLORS.get(int(f), 'gray') for f in families]

    # Make sample2 colors slightly darker
    colors2 = [f"{c}CC" if c != 'gray' else 'darkgray' for c in colors2]

    axes[0].bar(x - width/2, pct1, width, label=sample1, color=colors1, edgecolor='black')
    axes[0].bar(x + width/2, pct2, width, label=sample2, color=colors2, edgecolor='black')

    axes[0].set_xlabel('Family')
    axes[0].set_ylabel('Percentage (%)')
    axes[0].set_title('Family Composition Comparison')
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([f'F{int(f)}' for f in families], rotation=45)
    axes[0].legend()
    axes[0].grid(True, alpha=0.3, axis='y')

    # Difference plot
    diffs = df_comp['pct_diff'].values
    colors_diff = ['red' if d < 0 else 'blue' for d in diffs]

    axes[1].barh(x, diffs, color=colors_diff, edgecolor='black', alpha=0.7)
    axes[1].set_yticks(x)
    axes[1].set_yticklabels([f'F{int(f)}' for f in families])
    axes[1].set_xlabel(f'Percentage Difference\n({sample2} - {sample1})')
    axes[1].set_title('Family Composition Differences')
    axes[1].axvline(0, color='black', linestyle='-', linewidth=1)
    axes[1].grid(True, alpha=0.3, axis='x')

    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label=f'Enriched in {sample2}'),
                      Patch(facecolor='red', label=f'Depleted in {sample2}')]
    axes[1].legend(handles=legend_elements, loc='best')

    plt.tight_layout()
    plt.savefig(output_dir / 'family_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: family_comparison.png")

def plot_heterogeneity_comparison(het1, het2, sample1, sample2, output_dir):
    """Plot heterogeneity metric comparisons."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    metrics = ['n_families', 'shannon_entropy', 'simpson_diversity']
    titles = ['Number of Families per Array', 'Shannon Entropy', 'Simpson Diversity']

    for ax, metric, title in zip(axes, metrics, titles):
        data = [het1[metric].values, het2[metric].values]
        bp = ax.boxplot(data, labels=[sample1, sample2], patch_artist=True)

        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')

        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(title)
        ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(output_dir / 'heterogeneity_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: heterogeneity_comparison.png")

def plot_transition_comparison(df_trans, sample1, sample2, output_dir, top_n=20):
    """Plot top transition differences."""
    top_diffs = df_trans.head(top_n)

    fig, ax = plt.subplots(figsize=(10, 8))

    y_pos = np.arange(len(top_diffs))
    diffs = top_diffs['pct_diff'].values

    colors = ['red' if d < 0 else 'blue' for d in diffs]
    bars = ax.barh(y_pos, diffs, color=colors, edgecolor='black', alpha=0.7)

    # Labels
    labels = [f"F{row['from_family']}→F{row['to_family']}"
             for _, row in top_diffs.iterrows()]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel(f'Percentage Difference\n({sample2} - {sample1})')
    ax.set_title(f'Top {top_n} Transition Pattern Differences')
    ax.axvline(0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, axis='x')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label=f'More frequent in {sample2}'),
                      Patch(facecolor='red', label=f'More frequent in {sample1}')]
    ax.legend(handles=legend_elements, loc='best')

    plt.tight_layout()
    plt.savefig(output_dir / 'transition_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: transition_comparison.png")

def save_comparison_report(df_comp, chi2, p_value, df_trans, het_comp,
                           sample1, sample2, output_dir):
    """Save comprehensive comparison report."""
    report_file = output_dir / 'comparison_report.txt'

    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SAMPLE COMPARISON REPORT\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"Sample 1: {sample1}\n")
        f.write(f"Sample 2: {sample2}\n\n")

        # Family composition
        f.write("FAMILY COMPOSITION COMPARISON\n")
        f.write("-" * 80 + "\n")
        f.write(f"Chi-square statistic: {chi2:.2f}\n")
        f.write(f"P-value: {p_value:.2e}\n")
        f.write(f"Significant difference: {'YES' if p_value < 0.05 else 'NO'}\n\n")

        f.write(f"{'Family':<8} {sample1+' %':>12} {sample2+' %':>12} {'Diff %':>12} {'Fold':>10}\n")
        f.write("-" * 80 + "\n")

        for _, row in df_comp.iterrows():
            fold = row['fold_change']
            fold_str = f"{fold:.2f}" if fold != float('inf') else "New"

            f.write(f"F{int(row['family']):<7} "
                   f"{row[f'{sample1}_pct']:>12.2f} "
                   f"{row[f'{sample2}_pct']:>12.2f} "
                   f"{row['pct_diff']:>12.2f} "
                   f"{fold_str:>10}\n")

        f.write("\n")

        # Heterogeneity metrics
        f.write("HETEROGENEITY METRICS COMPARISON\n")
        f.write("-" * 80 + "\n")

        for metric, comp in het_comp.items():
            f.write(f"\n{metric.upper().replace('_', ' ')}:\n")
            f.write(f"  {sample1}: {comp[f'{sample1}_mean']:.3f} ± {comp[f'{sample1}_std']:.3f}\n")
            f.write(f"  {sample2}: {comp[f'{sample2}_mean']:.3f} ± {comp[f'{sample2}_std']:.3f}\n")
            f.write(f"  t-statistic: {comp['t_statistic']:.3f}\n")
            f.write(f"  p-value: {comp['p_value']:.2e}\n")
            f.write(f"  Significant: {'YES' if comp['significant'] else 'NO'}\n")

        f.write("\n")

        # Top transition differences
        f.write("TOP 10 TRANSITION DIFFERENCES\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Transition':<15} {sample1+' %':>12} {sample2+' %':>12} {'Diff %':>12}\n")
        f.write("-" * 80 + "\n")

        for _, row in df_trans.head(10).iterrows():
            trans = f"F{row['from_family']}→F{row['to_family']}"
            f.write(f"{trans:<15} "
                   f"{row[f'{sample1}_pct']:>12.2f} "
                   f"{row[f'{sample2}_pct']:>12.2f} "
                   f"{row['pct_diff']:>12.2f}\n")

        f.write("\n")
        f.write("=" * 80 + "\n")

    print(f"  Saved: comparison_report.txt")

    # Save JSON
    json_data = {
        'sample1': sample1,
        'sample2': sample2,
        'family_composition': {
            'chi2_statistic': float(chi2),
            'p_value': float(p_value),
            'significant': p_value < 0.05,
            'families': df_comp.to_dict(orient='records')
        },
        'heterogeneity': het_comp,
        'transitions': df_trans.head(20).to_dict(orient='records')
    }

    json_file = output_dir / 'comparison_report.json'
    with open(json_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    print(f"  Saved: comparison_report.json")

def main():
    if len(sys.argv) != 6:
        print("Usage: python compare_samples.py <sample1.tsv> <sample2.tsv> <name1> <name2> <output_dir>")
        sys.exit(1)

    tsv1 = sys.argv[1]
    tsv2 = sys.argv[2]
    sample1 = sys.argv[3]
    sample2 = sys.argv[4]
    output_dir = Path(sys.argv[5])
    output_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*80)
    print("SAMPLE COMPARISON ANALYSIS")
    print("="*80 + "\n")

    # Load data
    print("Loading data...")
    df1 = load_sample_data(tsv1, sample1)
    df2 = load_sample_data(tsv2, sample2)

    # Compare family composition
    print("\nComparing family composition...")
    df_comp, chi2, p_value = compare_family_composition(df1, df2, sample1, sample2)
    print(f"  Chi-square: {chi2:.2f}, p-value: {p_value:.2e}")

    # Compare transitions
    print("\nComparing transition patterns...")
    df_trans = compare_transitions(df1, df2, sample1, sample2)

    # Compare heterogeneity
    print("\nComparing heterogeneity metrics...")
    het_comp, het1, het2 = compare_heterogeneity(df1, df2, sample1, sample2)

    # Generate plots
    print("\nGenerating plots...")
    plot_family_comparison(df_comp, sample1, sample2, output_dir)
    plot_heterogeneity_comparison(het1, het2, sample1, sample2, output_dir)
    plot_transition_comparison(df_trans, sample1, sample2, output_dir)

    # Save reports
    print("\nSaving reports...")
    save_comparison_report(df_comp, chi2, p_value, df_trans, het_comp,
                          sample1, sample2, output_dir)

    # Save tables
    df_comp.to_csv(output_dir / 'family_comparison.tsv', sep='\t', index=False)
    print("  Saved: family_comparison.tsv")

    df_trans.to_csv(output_dir / 'transition_comparison.tsv', sep='\t', index=False)
    print("  Saved: transition_comparison.tsv")

    print("\n" + "="*80)
    print("COMPARISON COMPLETE")
    print(f"Results saved to: {output_dir}")
    print("="*80 + "\n")

if __name__ == '__main__':
    main()
