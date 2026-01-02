/*
========================================================================================
    RIBBON_PLOTS: Single molecule vs reference satellite visualization
========================================================================================
    Creates ribbon plots showing how satellites are remodelled in individual reads
    compared to the reference genome. Shows:
    - Reference track (top): monomers from reference
    - Read track (bottom): monomers from sequenced molecule
    - Gray ribbons: aligned regions
    - Gaps: large indels showing satellite remodelling
    - Colors: CEN178 families
========================================================================================
*/

process RIBBON_PLOTS {
    tag "ribbon plots"
    publishDir "${params.outdir}/05_plots/ribbon_plots", mode: 'copy'

    input:
    path bam_file
    path bai_file
    path indel_catalog
    path monomer_classifications
    path monomer_info

    output:
    path "ribbon_*.png", emit: plots, optional: true
    path "ribbon_plots.log", emit: log

    when:
    params.mode == 'reads' && params.generate_ribbon_plots

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== Ribbon Plot Generation ===" > ribbon_plots.log
    echo "BAM: ${bam_file}" >> ribbon_plots.log
    echo "Indel catalog: ${indel_catalog}" >> ribbon_plots.log
    echo "Monomers: ${monomer_classifications}" >> ribbon_plots.log

    # Prepare monomer data with read_id column
    python3 ${projectDir}/bin/prepare_monomer_data_for_viz.py \\
        ${monomer_classifications} \\
        monomers_with_readid.tsv

    # Find top reads with most large indels
    awk -F'\t' 'NR>1 && \$6>=100 {count[\$1]++} END {for (id in count) print count[id], id}' \\
        ${indel_catalog} | sort -rn | head -10 | awk '{print \$2}' > top_indel_reads.txt

    n_reads=\$(wc -l < top_indel_reads.txt)
    echo "Found \$n_reads reads with large indels" >> ribbon_plots.log

    if [ \$n_reads -eq 0 ]; then
        echo "No reads with large indels found - skipping ribbon plots" >> ribbon_plots.log
        exit 0
    fi

    # Limit to top 5 to avoid excessive runtime
    head -5 top_indel_reads.txt > selected_reads.txt

    # Generate ribbon plots
    read_ids=\$(tr '\\n' ' ' < selected_reads.txt)
    echo "Generating ribbon plots for: \$read_ids" >> ribbon_plots.log

    python3 ${projectDir}/bin/visualize_indel_families_v2.py \\
        --bam ${bam_file} \\
        --sv-info ${indel_catalog} \\
        --monomers monomers_with_readid.tsv \\
        --read-ids \$read_ids \\
        --output-dir . \\
        2>&1 | tee -a ribbon_plots.log || true

    # Rename outputs to match expected pattern
    for f in indel_families_*.png; do
        if [ -f "\$f" ]; then
            new_name=\$(echo "\$f" | sed 's/indel_families/ribbon/')
            mv "\$f" "\$new_name"
            echo "  Generated: \$new_name" >> ribbon_plots.log
        fi
    done

    n_plots=\$(ls -1 ribbon_*.png 2>/dev/null | wc -l)
    echo "Generated \$n_plots ribbon plots" >> ribbon_plots.log

    echo "✅ Ribbon plot generation complete" >> ribbon_plots.log
    """
}
