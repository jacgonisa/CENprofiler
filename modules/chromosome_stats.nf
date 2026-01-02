/*
========================================================================================
    CHROMOSOME_STATS: Generate chromosome-aware statistics
========================================================================================
    Calculates per-chromosome statistics for:
    - Family distribution
    - HOR prevalence
    - Monomer counts
========================================================================================
*/

process CHROMOSOME_STATS {
    tag "chromosome stats"
    publishDir "${params.outdir}/04_stats", mode: 'copy'

    input:
    path classifications
    path hors

    output:
    path "chromosome_stats.tsv", emit: stats
    path "family_by_chromosome.tsv", emit: family_stats
    path "hor_by_chromosome.tsv", emit: hor_stats

    when:
    params.chromosome_aware

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys

    # Load data
    df_class = pd.read_csv("${classifications}", sep='\\t')
    df_hors = pd.read_csv("${hors}", sep='\\t')

    # Extract chromosome from seq_id (assuming format like "Chr1_...")
    df_class['chromosome'] = df_class['seq_id'].str.extract(r'(Chr\\d+)')
    df_hors['chromosome'] = df_hors['seq_id'].str.extract(r'(Chr\\d+)')

    # Overall stats per chromosome
    stats = []
    for chrom in sorted(df_class['chromosome'].dropna().unique()):
        chrom_monomers = df_class[df_class['chromosome'] == chrom]
        chrom_hors = df_hors[df_hors['chromosome'] == chrom]

        stats.append({
            'chromosome': chrom,
            'total_monomers': len(chrom_monomers),
            'classified_monomers': chrom_monomers['monomer_family'].notna().sum(),
            'total_hors': len(chrom_hors),
            'homHORs': (chrom_hors['hor_type'] == 'homHOR').sum() if len(chrom_hors) > 0 else 0,
            'hetHORs': (chrom_hors['hor_type'] == 'hetHOR').sum() if len(chrom_hors) > 0 else 0
        })

    df_stats = pd.DataFrame(stats)
    df_stats.to_csv("chromosome_stats.tsv", sep='\\t', index=False)

    # Family distribution per chromosome
    family_stats = df_class.groupby(['chromosome', 'monomer_family']).size().reset_index(name='count')
    family_stats.to_csv("family_by_chromosome.tsv", sep='\\t', index=False)

    # HOR distribution per chromosome
    if len(df_hors) > 0:
        hor_stats = df_hors.groupby(['chromosome', 'hor_unit']).size().reset_index(name='count')
    else:
        hor_stats = pd.DataFrame(columns=['chromosome', 'hor_unit', 'count'])

    hor_stats.to_csv("hor_by_chromosome.tsv", sep='\\t', index=False)

    print("Chromosome statistics completed", file=sys.stderr)
    """
}
