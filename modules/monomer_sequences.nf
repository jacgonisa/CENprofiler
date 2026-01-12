#!/usr/bin/env nextflow

/*
 * Nextflow module for monomer sequence extraction and analysis
 */

process EXTRACT_MONOMER_SEQUENCES {
    publishDir "${params.outdir}/08_monomer_sequences", mode: 'copy'

    input:
    path monomer_classifications
    path monomers_fasta

    output:
    path "family_fastas/*.fa"           , emit: family_fastas
    path "family_diversity.tsv"         , emit: diversity
    path "family_diversity_report.txt"  , emit: diversity_report
    path "sequence_summary.txt"         , emit: summary
    path "consensus/*.fa"               , emit: consensus, optional: true
    path "*.log"                        , emit: log, optional: true

    script:
    """
    python ${projectDir}/bin/extract_monomer_sequences.py \\
        ${monomer_classifications} \\
        ${monomers_fasta} \\
        . \\
        > sequence_extraction.log 2>&1
    """
}
