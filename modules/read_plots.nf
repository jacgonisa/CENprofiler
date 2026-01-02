/*
========================================================================================
    READ_PLOTS: Generate read-mode visualizations
========================================================================================
    Creates plots for read-level analysis:
    - Family distribution across reads
    - Indel size distributions
    - Family-specific indel patterns
========================================================================================
*/

process READ_PLOTS {
    tag "read plots"
    publishDir "${params.outdir}/05_plots/reads", mode: 'copy'

    input:
    path classifications
    path monomer_info
    path indel_stats

    output:
    path "*.png", emit: plots, optional: true
    path "read_plots.log", emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== Read Mode Visualization ===" > read_plots.log
    echo "Classifications: ${classifications}" >> read_plots.log
    echo "Monomer info: ${monomer_info}" >> read_plots.log
    echo "Indel stats: ${indel_stats}" >> read_plots.log

    # Count reads with monomers
    n_reads=\$(tail -n +2 ${classifications} | cut -f2 | sort -u | wc -l)
    n_monomers=\$(tail -n +2 ${classifications} | wc -l)
    echo "Reads analyzed: \$n_reads" >> read_plots.log
    echo "Total monomers: \$n_monomers" >> read_plots.log

    if [ \$n_monomers -eq 0 ]; then
        echo "WARNING: No monomers found - skipping plots" >> read_plots.log
        exit 0
    fi

    # Generate family distribution plot using Python
    python3 << 'PYEOF'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df = pd.read_csv('${classifications}', sep='\\t')

# Family distribution
families = df['monomer_family'].dropna()
if len(families) > 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    family_counts = families.value_counts().sort_index()

    ax.bar(family_counts.index.astype(str), family_counts.values,
           color='steelblue', edgecolor='black', alpha=0.8)
    ax.set_xlabel('Family', fontweight='bold', fontsize=12)
    ax.set_ylabel('Count', fontweight='bold', fontsize=12)
    ax.set_title('CEN178 Family Distribution in Reads', fontweight='bold', fontsize=14)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('family_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Generated family_distribution.png")
else:
    print("No families classified")

# Per-read statistics
reads = df.groupby('seq_id').agg({
    'monomer_family': lambda x: x.notna().sum(),
    'monomer_idx': 'count'
}).reset_index()
reads.columns = ['read_id', 'classified_monomers', 'total_monomers']

if len(reads) > 0:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Total monomers per read
    ax1.hist(reads['total_monomers'], bins=30, color='forestgreen',
             edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Monomers per read', fontweight='bold', fontsize=11)
    ax1.set_ylabel('Frequency', fontweight='bold', fontsize=11)
    ax1.set_title('Array Size Distribution', fontweight='bold', fontsize=12)
    ax1.axvline(reads['total_monomers'].median(), color='red',
                linestyle='--', linewidth=2, label=f'Median: {reads["total_monomers"].median():.0f}')
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Classification rate per read
    reads['classification_rate'] = reads['classified_monomers'] / reads['total_monomers'] * 100
    ax2.hist(reads['classification_rate'], bins=20, color='coral',
             edgecolor='black', alpha=0.7)
    ax2.set_xlabel('Classification rate (%)', fontweight='bold', fontsize=11)
    ax2.set_ylabel('Frequency', fontweight='bold', fontsize=11)
    ax2.set_title('Monomer Classification Rate', fontweight='bold', fontsize=12)
    ax2.axvline(reads['classification_rate'].median(), color='darkred',
                linestyle='--', linewidth=2, label=f'Median: {reads["classification_rate"].median():.1f}%')
    ax2.legend()
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('read_statistics.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Generated read_statistics.png")
PYEOF

    # Generate indel plots if indel stats exist
    if [ -f "${indel_stats}" ]; then
        n_indels=\$(tail -n +2 ${indel_stats} | wc -l)
        echo "Indels detected: \$n_indels" >> read_plots.log

        if [ \$n_indels -gt 0 ]; then
            echo "Generating indel distribution plots..." >> read_plots.log

            python3 << 'PYEOF2'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load indel stats
indels = pd.read_csv('${indel_stats}', sep='\\t')

if len(indels) > 0:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Indel size distribution
    ax1.hist(indels['size'], bins=50, color='indianred',
             edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Indel size (bp)', fontweight='bold', fontsize=11)
    ax1.set_ylabel('Frequency', fontweight='bold', fontsize=11)
    ax1.set_title('Indel Size Distribution', fontweight='bold', fontsize=12)
    ax1.set_yscale('log')
    ax1.grid(axis='y', alpha=0.3)

    # Indel type counts
    if 'type' in indels.columns:
        type_counts = indels['type'].value_counts()
        ax2.bar(type_counts.index, type_counts.values,
                color=['steelblue', 'coral'], edgecolor='black', alpha=0.8)
        ax2.set_xlabel('Indel type', fontweight='bold', fontsize=11)
        ax2.set_ylabel('Count', fontweight='bold', fontsize=11)
        ax2.set_title('Insertions vs Deletions', fontweight='bold', fontsize=12)
        ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('indel_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("✅ Generated indel_distribution.png")
PYEOF2
        fi
    fi

    # Count outputs
    n_plots=\$(ls -1 *.png 2>/dev/null | wc -l)
    echo "Generated \$n_plots read plots" >> read_plots.log

    echo "✅ Read visualization complete" >> read_plots.log
    """
}
