#!/usr/bin/env nextflow

/*
 * Nextflow module for comprehensive monomer statistics analysis
 */

process MONOMER_STATISTICS {
    publishDir "${params.outdir}/07_monomer_statistics", mode: 'copy'

    input:
    path monomer_classifications

    output:
    path "*.png"                        , emit: plots
    path "*.tsv"                        , emit: tables
    path "monomer_statistics.txt"       , emit: report
    path "monomer_statistics.json"      , emit: json
    path "*.log"                        , emit: log, optional: true

    script:
    """
    python ${projectDir}/bin/analyze_monomer_statistics.py \\
        ${monomer_classifications} \\
        . \\
        > monomer_statistics.log 2>&1
    """
}
