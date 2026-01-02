/*
========================================================================================
    HOR_PLOTS: Generate HOR-specific visualizations
========================================================================================
    Creates comprehensive plots for HOR analysis:
    - Genome-wide HOR distribution across chromosomes
    - HOR architecture schematics (homHOR and hetHOR)
    - Large duplication visualizations
========================================================================================
*/

process HOR_PLOTS {
    tag "HOR plots"
    publishDir "${params.outdir}/05_plots/hors", mode: 'copy'

    input:
    path hors
    path large_duplications
    path chromosome_stats

    output:
    path "*.png", emit: plots, optional: true
    path "hor_plots.log", emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== HOR Visualization ===" > hor_plots.log
    echo "HORs file: ${hors}" >> hor_plots.log
    echo "Large duplications: ${large_duplications}" >> hor_plots.log

    # Count HORs
    n_hors=\$(tail -n +2 ${hors} | wc -l)
    echo "Total HORs: \$n_hors" >> hor_plots.log

    if [ \$n_hors -eq 0 ]; then
        echo "WARNING: No HORs found - skipping plots" >> hor_plots.log
        exit 0
    fi

    # Create expected filename symlinks for scripts
    ln -sf ${hors} reference_genome_hors_MONOMER_LEVEL.tsv
    cp ${large_duplications} large_duplications_input.tsv
    mv large_duplications_input.tsv large_duplications.tsv

    # 1. Genome-wide HOR distribution
    if [ -f "${projectDir}/bin/plot_monomer_level_genome_wide.py" ]; then
        echo "Generating genome-wide HOR distribution plot..." >> hor_plots.log
        python3 "${projectDir}/bin/plot_monomer_level_genome_wide.py" \
            2>&1 | tee -a hor_plots.log || echo "Genome-wide plot failed" >> hor_plots.log
    fi

    # 2. HOR schematics (homHOR and hetHOR architecture)
    if [ -f "${projectDir}/bin/plot_monomer_level_schematics.py" ]; then
        echo "Generating HOR architecture schematics..." >> hor_plots.log
        python3 "${projectDir}/bin/plot_monomer_level_schematics.py" \
            2>&1 | tee -a hor_plots.log || echo "Schematic plots failed" >> hor_plots.log
    fi

    # 3. Large duplication overview
    n_large_dups=\$(tail -n +2 ${large_duplications} | wc -l)
    echo "Large duplications: \$n_large_dups" >> hor_plots.log

    if [ \$n_large_dups -gt 0 ] && [ -f "${projectDir}/bin/plot_large_duplications_overview.py" ]; then
        echo "Generating large duplication overview..." >> hor_plots.log
        python3 "${projectDir}/bin/plot_large_duplications_overview.py" \
            2>&1 | tee -a hor_plots.log || echo "Duplication overview failed" >> hor_plots.log
    fi

    # 4. Large duplication details (for top duplications)
    if [ \$n_large_dups -gt 0 ] && [ -f "${projectDir}/bin/plot_large_duplications_detail.py" ]; then
        echo "Generating large duplication detail plots..." >> hor_plots.log
        python3 "${projectDir}/bin/plot_large_duplications_detail.py" \
            2>&1 | tee -a hor_plots.log || echo "Duplication detail failed" >> hor_plots.log
    fi

    # Count outputs
    n_plots=\$(ls -1 *.png 2>/dev/null | wc -l)
    echo "Generated \$n_plots HOR plots" >> hor_plots.log

    echo "✅ HOR visualization complete" >> hor_plots.log
    """
}
