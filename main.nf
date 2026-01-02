#!/usr/bin/env nextflow

/*
========================================================================================
    CENprofiler: Centromeric Satellite & HOR Analysis Pipeline
========================================================================================
    Github : https://github.com/jacgonisa/CENprofiler

    Main workflow entry point with AUTO-DETECTION:

    Detection Logic:
        - If --alignment provided → READ MODE with indel analysis
        - Otherwise → GENOME MODE
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
     CENprofiler v2.0 (Auto-Detect Mode)
    ========================================

    Usage:

      # GENOME MODE (auto-detected)
      nextflow run main.nf --input genome.fasta [options]

      # READ MODE with indel analysis (auto-detected when --alignment provided)
      nextflow run main.nf --input reads.fa --alignment reads.bam --reference_genome genome.fa [options]

    Required Arguments:
      --input             Input FASTA file (genome or reads)
      --reference_monomers Reference monomer FASTA for classification
      --family_assignments Family assignment file (TSV: monomer_id<tab>family_id)
      --outdir            Output directory

    Read Mode with BAM Arguments (triggers indel analysis):
      --alignment         BAM/SAM alignment file (AUTO-TRIGGERS read mode with indels)
      --reference_genome  Reference genome FASTA (required with --alignment)
      --annotation_dir    Genomic annotations directory (centromeres, rDNA BED files)
      --sample_name       Sample name for output files (default: 'sample')
      --min_indel_size    Minimum indel size to extract reads (default: 100bp)

    FasTAN Parameters:
      --period_min        Minimum tandem repeat period (default: 160bp)
      --period_max        Maximum tandem repeat period (default: 200bp)
      --fastan_threads    Threads for FasTAN (default: 8)

    Classification Parameters:
      --min_identity      Minimum alignment identity % (default: 70)
      --minimap2_threads  Threads for minimap2 (default: 4)

    Genome Mode - HOR Detection:
      --min_copies        Minimum HOR copies (default: 3)
      --min_monomers      Minimum monomers per HOR unit (default: 3)
      --max_gap           Maximum gap between monomers in bp (default: 500)
      --large_dup_threshold Large duplication threshold in kb (default: 40)

    Examples:

      # 1. Analyze reference genome for HORs
      nextflow run main.nf \\
        --input TAIR12/GCA_028009825.2_Col-CC_genomic.fna \\
        --reference_monomers Col-CC-V2-CEN178-representative.fasta \\
        --family_assignments itol_manual_phylo_clusters.txt \\
        --outdir results_genome/

      # 2. Analyze reads with large indels from BAM
      nextflow run main.nf \\
        --input dummy.fa \\
        --alignment tair12_indel_comparison/results/mapping/Col_9day.bam \\
        --reference_genome TAIR12/GCA_028009825.2_Col-CC_genomic.fna \\
        --annotation_dir TAIR12/curated_anno \\
        --reference_monomers Col-CC-V2-CEN178-representative.fasta \\
        --family_assignments itol_manual_phylo_clusters.txt \\
        --sample_name Col_9day \\
        --min_indel_size 100 \\
        --outdir results_reads/

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

// Check required parameters
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
    error "Error: --outdir is required (use './results' for default)"
}

// Check reference files exist
reference_file = file(params.reference_monomers)
if (!reference_file.exists()) {
    error "Error: Reference monomers file does not exist: ${params.reference_monomers}"
}

families_file = file(params.family_assignments)
if (!families_file.exists()) {
    error "Error: Family assignments file does not exist: ${params.family_assignments}"
}

// Validate read mode with alignment
if (params.alignment) {
    // READ MODE with indel analysis
    alignment_file = file(params.alignment)
    if (!alignment_file.exists()) {
        error "Error: Alignment file does not exist: ${params.alignment}"
    }

    if (!params.reference_genome) {
        error "Error: --reference_genome is required when --alignment is provided"
    }

    genome_file = file(params.reference_genome)
    if (!genome_file.exists()) {
        error "Error: Reference genome file does not exist: ${params.reference_genome}"
    }

    if (!params.annotation_dir) {
        error "Error: --annotation_dir is required for read mode with alignment (e.g., TAIR12/curated_anno)"
    }

    anno_dir = file(params.annotation_dir)
    if (!anno_dir.exists() || !anno_dir.isDirectory()) {
        error "Error: Annotation directory does not exist or is not a directory: ${params.annotation_dir}"
    }
}

/*
========================================================================================
    NAMED WORKFLOWS
========================================================================================
*/

include { GENOME_MODE } from './workflows/genome_mode'
include { READ_MODE_WITH_INDELS } from './workflows/read_mode_with_indels'

/*
========================================================================================
    AUTO-DETECT MODE & RUN WORKFLOW
========================================================================================
*/

workflow {

    // AUTO-DETECTION LOGIC
    def detected_mode = params.alignment ? "READ MODE with INDEL ANALYSIS" : "GENOME MODE"

    log.info """
    ========================================
     CENprofiler Pipeline v2.0
    ========================================
    Mode (AUTO-DETECTED): ${detected_mode}
    Input             : ${params.input}
    Reference Monomers: ${params.reference_monomers}
    Family Assignments: ${params.family_assignments}
    Output Directory  : ${params.outdir}
    ========================================
    """.stripIndent()

    // Create common input channels
    ch_reference = Channel.fromPath(params.reference_monomers)
    ch_families = Channel.fromPath(params.family_assignments)

    // Branch based on detected mode
    if (params.alignment) {
        //
        // READ MODE WITH INDEL ANALYSIS (BAM provided)
        //
        log.info """
        BAM Alignment     : ${params.alignment}
        Reference Genome  : ${params.reference_genome}
        Annotation Dir    : ${params.annotation_dir}
        Sample Name       : ${params.sample_name}
        Min Indel Size    : ${params.min_indel_size} bp
        ========================================
        """.stripIndent()

        // Create channels for read mode
        ch_bam = Channel.fromPath(params.alignment)
        ch_bai = Channel.fromPath("${params.alignment}.bai")
        ch_genome = Channel.fromPath(params.reference_genome)
        ch_anno_dir = Channel.fromPath(params.annotation_dir)
        ch_sample = Channel.value(params.sample_name)

        // Run READ MODE WITH INDELS
        READ_MODE_WITH_INDELS(
            ch_bam,
            ch_bai,
            ch_genome,
            ch_reference,
            ch_families,
            ch_anno_dir,
            ch_sample
        )

    } else {
        //
        // GENOME MODE (default)
        //
        log.info "Running genome-wide analysis for HORs and satellites\n========================================\n"

        // Check input exists
        input_file = file(params.input)
        if (!input_file.exists()) {
            error "Error: Input file does not exist: ${params.input}"
        }

        ch_input = Channel.fromPath(params.input)

        // Run GENOME MODE
        GENOME_MODE(
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
