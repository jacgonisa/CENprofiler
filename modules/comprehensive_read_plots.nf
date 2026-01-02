/*
========================================================================================
    COMPREHENSIVE_READ_PLOTS: All read-mode visualizations from results_v2
========================================================================================
    Generates complete visualization suite:
    1. Family distribution bar chart
    2. Family transition heatmap
    3. Top arrays combined view
    4. Individual array detail plots (top 5)
    5. Summary statistics
========================================================================================
*/

process COMPREHENSIVE_READ_PLOTS {
    tag "comprehensive read plots"
    publishDir "${params.outdir}/05_plots/comprehensive", mode: 'copy'

    input:
    path classifications
    path monomer_info

    output:
    path "family_summary.png", emit: family_summary, optional: true
    path "top_arrays_combined.png", emit: top_arrays, optional: true
    path "array_*.png", emit: array_plots, optional: true
    path "ARRAY_SUMMARY.txt", emit: summary, optional: true
    path "comprehensive_plots.log", emit: log

    when:
    params.alignment != null  // Run in read mode with alignment

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== Comprehensive Read Visualization ===" > comprehensive_plots.log
    echo "Classifications: ${classifications}" >> comprehensive_plots.log

    n_monomers=\$(tail -n +2 ${classifications} | wc -l)
    echo "Total monomers: \$n_monomers" >> comprehensive_plots.log

    if [ \$n_monomers -eq 0 ]; then
        echo "WARNING: No monomers found - skipping plots" >> comprehensive_plots.log
        exit 0
    fi

    # Generate comprehensive plots
    python3 ${projectDir}/bin/create_comprehensive_read_plots.py \\
        --classifications ${classifications} \\
        --output-dir . \\
        --n-arrays 5 \\
        2>&1 | tee -a comprehensive_plots.log

    # Count outputs
    n_plots=\$(ls -1 *.png 2>/dev/null | wc -l)
    echo "Generated \$n_plots comprehensive plots" >> comprehensive_plots.log

    echo "✅ Comprehensive visualization complete" >> comprehensive_plots.log
    """
}
