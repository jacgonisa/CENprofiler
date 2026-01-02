/*
========================================================================================
    READ MODE WITH INDELS WORKFLOW
========================================================================================
    Analyzes reads with large indels from BAM alignment:
    - Extracts reads with large indels (≥ min_indel_size)
    - Detects tandem satellite arrays (FasTAN)
    - Classifies monomers by family
    - Analyzes indel-family relationships    - Generates visualizations
========================================================================================
*/

include { LOAD_GENOMIC_REGIONS } from '../modules/load_genomic_regions'
include { EXTRACT_READS_FROM_BAM } from '../modules/extract_reads_from_bam'
include { FASTAN } from '../modules/fastan'
include { TANBED } from '../modules/tanbed'
include { EXTRACT_MONOMERS } from '../modules/extract_monomers'
include { CLASSIFY_MONOMERS } from '../modules/classify_monomers'
include { READ_PLOTS } from '../modules/read_plots'
include { COMPREHENSIVE_READ_PLOTS } from '../modules/comprehensive_read_plots'
include { ANALYZE_DELETION_MONOMERS } from '../modules/analyze_deletion_monomers'
include { RIBBON_PLOTS } from '../modules/ribbon_plots'

workflow READ_MODE_WITH_INDELS {
    take:
    bam_file           // Channel: BAM alignment file
    bam_index          // Channel: BAM index (.bai)
    reference_genome   // Channel: Reference genome FASTA
    reference_monomers // Channel: Reference monomer FASTA
    family_assignments // Channel: Family assignment TSV
    annotation_dir     // Channel: Annotation directory path
    sample_name        // Channel: Sample name

    main:

    log.info """
    ========================================
     Running READ MODE with INDEL ANALYSIS
    ========================================
    Sample: ${sample_name}
    BAM file: ${bam_file}
    Reference genome: ${reference_genome}
    Min indel size: ${params.min_indel_size} bp
    ========================================
    """

    //
    // STEP 1: Load genomic region annotations
    //
    LOAD_GENOMIC_REGIONS(
        annotation_dir
    )

    //
    // STEP 2: Extract reads with large indels from BAM
    //
    EXTRACT_READS_FROM_BAM(
        bam_file,
        bam_index,
        LOAD_GENOMIC_REGIONS.out.regions,
        sample_name
    )

    //
    // STEP 3: Run FasTAN to detect tandem repeats in extracted reads
    //
    FASTAN(
        EXTRACT_READS_FROM_BAM.out.reads_fasta
    )

    //
    // STEP 4: Convert .1aln to BED format
    //
    TANBED(
        FASTAN.out.aln
    )

    //
    // STEP 5: Extract individual monomers from tandem arrays
    //
    EXTRACT_MONOMERS(
        EXTRACT_READS_FROM_BAM.out.reads_fasta,
        TANBED.out.bed,
        'reads'  // mode
    )

    //
    // STEP 6: Classify monomers using minimap2 + family assignments
    //
    CLASSIFY_MONOMERS(
        EXTRACT_MONOMERS.out.monomers_fasta,
        reference_monomers,
        family_assignments,
        EXTRACT_MONOMERS.out.monomer_info
    )

    //
    // STEP 7: Generate basic read-level visualizations
    //
    indel_stats_ch = EXTRACT_READS_FROM_BAM.out.indel_catalog

    READ_PLOTS(
        CLASSIFY_MONOMERS.out.classifications,
        EXTRACT_MONOMERS.out.monomer_info,
        indel_stats_ch
    )

    //
    // STEP 8: Generate comprehensive visualizations (results_v2 style)
    //
    COMPREHENSIVE_READ_PLOTS(
        CLASSIFY_MONOMERS.out.classifications,
        EXTRACT_MONOMERS.out.monomer_info
    )

    //
    // STEP 9: Analyze deletion monomers from reference genome
    //
    ANALYZE_DELETION_MONOMERS(
        bam_file,
        bam_index,
        reference_genome,
        reference_genome.map { file -> file.parent.resolve("${file.baseName}.fai") },
        EXTRACT_READS_FROM_BAM.out.indel_catalog,
        reference_monomers,
        family_assignments,
        sample_name
    )

    //
    // STEP 10: Generate ribbon plots showing satellite remodelling
    //
    RIBBON_PLOTS(
        bam_file,
        bam_index,
        EXTRACT_READS_FROM_BAM.out.indel_catalog,
        CLASSIFY_MONOMERS.out.classifications,
        EXTRACT_MONOMERS.out.monomer_info
    )

    emit:
    regions               = LOAD_GENOMIC_REGIONS.out.regions
    reads_fasta           = EXTRACT_READS_FROM_BAM.out.reads_fasta
    indel_catalog         = EXTRACT_READS_FROM_BAM.out.indel_catalog
    extraction_stats      = EXTRACT_READS_FROM_BAM.out.stats
    fastan_bed            = TANBED.out.bed
    monomer_fasta         = EXTRACT_MONOMERS.out.monomers_fasta
    classifications       = CLASSIFY_MONOMERS.out.classifications
    plots                 = READ_PLOTS.out.plots
    comprehensive_plots   = COMPREHENSIVE_READ_PLOTS.out.family_summary
    deletion_monomers     = ANALYZE_DELETION_MONOMERS.out.deletion_monomers
    ribbon_plots          = RIBBON_PLOTS.out.plots
}
