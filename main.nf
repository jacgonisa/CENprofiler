#!/usr/bin/env nextflow

/*
========================================================================================
    CENprofiler: Centromeric Satellite & HOR Analysis Pipeline
========================================================================================
    Github : https://github.com/yourusername/CENprofiler

    Main workflow entry point

    Modes:
        - genome: Analyze reference genome for satellites and HORs
        - reads:  Analyze long reads for satellite composition and variants
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    PRINT HELP
========================================================================================
*/

def helpMessage() {
    log.info"""

    ========================================
     CENprofiler v1.0
    ========================================

    Usage:

      nextflow run main.nf --mode <genome|reads> [options]

    Required Arguments:
      --mode              Analysis mode: 'genome' or 'reads'
      --input             Input file (FASTA for genome mode, FASTA/BAM for read mode)
      --reference_monomers Reference monomer FASTA for classification
      --family_assignments Family assignment file (TSV: monomer_id<tab>family_id)
      --outdir            Output directory

    Optional Arguments:
      --period_min        Minimum tandem repeat period (default: 160)
      --period_max        Maximum tandem repeat period (default: 200)
      --min_identity      Minimum alignment identity for classification (default: 70)
      --fastan_threads    Threads for FasTAN (default: 8)

    Genome Mode Options:
      --min_copies        Minimum HOR copies (default: 3)
      --min_monomers      Minimum monomers per HOR unit (default: 3)
      --max_gap           Maximum gap between consecutive monomers (default: 500)
      --large_dup_threshold Large duplication threshold in kb (default: 40)

    Read Mode Options:
      --analyze_indels    Perform indel analysis (default: true)
      --min_array_size    Minimum array size to analyze (default: 5 monomers)

    Examples:

      # Genome mode
      nextflow run main.nf \\
        --mode genome \\
        --input genome.fasta \\
        --reference_monomers representatives.fasta \\
        --family_assignments families.txt \\
        --outdir results/

      # Read mode
      nextflow run main.nf \\
        --mode reads \\
        --input reads.fasta \\
        --reference_monomers representatives.fasta \\
        --family_assignments families.txt \\
        --outdir results/

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

if (!params.mode) {
    error "Error: --mode is required (genome or reads)"
}

if (!params.input) {
    error "Error: --input is required"
}

if (!params.reference_monomers) {
    error "Error: --reference_monomers is required"
}

if (!params.family_assignments) {
    error "Error: --family_assignments is required"
}

if (!params.outdir) {
    error "Error: --outdir is required"
}

// Validate mode
if (params.mode != 'genome' && params.mode != 'reads') {
    error "Error: --mode must be 'genome' or 'reads'"
}

// Check input files exist
input_file = file(params.input)
if (!input_file.exists()) {
    error "Error: Input file does not exist: ${params.input}"
}

reference_file = file(params.reference_monomers)
if (!reference_file.exists()) {
    error "Error: Reference monomers file does not exist: ${params.reference_monomers}"
}

families_file = file(params.family_assignments)
if (!families_file.exists()) {
    error "Error: Family assignments file does not exist: ${params.family_assignments}"
}

/*
========================================================================================
    NAMED WORKFLOWS
========================================================================================
*/

include { GENOME_MODE } from './workflows/genome_mode'
include { READ_MODE } from './workflows/read_mode'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    log.info """
    ========================================
     CENprofiler Pipeline
    ========================================
    Mode              : ${params.mode}
    Input             : ${params.input}
    Reference Monomers: ${params.reference_monomers}
    Family Assignments: ${params.family_assignments}
    Output Directory  : ${params.outdir}
    ========================================
    """.stripIndent()

    // Create input channel
    ch_input = Channel.fromPath(params.input)
    ch_reference = Channel.fromPath(params.reference_monomers)
    ch_families = Channel.fromPath(params.family_assignments)

    // Run appropriate workflow based on mode
    if (params.mode == 'genome') {
        GENOME_MODE(
            ch_input,
            ch_reference,
            ch_families
        )
    } else if (params.mode == 'reads') {
        READ_MODE(
            ch_input,
            ch_reference,
            ch_families
        )
    }
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================
    Pipeline completed!
    ========================================
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Work Dir  : ${workflow.workDir}
    Results   : ${params.outdir}
    Duration  : ${workflow.duration}
    ========================================
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed!"
    log.error "Error: ${workflow.errorMessage}"
}
