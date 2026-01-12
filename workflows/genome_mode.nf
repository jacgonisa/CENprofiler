/*
========================================================================================
    GENOME MODE WORKFLOW
========================================================================================
    Analyzes reference genome for:
    - Tandem satellite arrays (FasTAN)
    - Monomer classification
    - Higher-Order Repeats (HORs)
    - Chromosome-level statistics
    - Genome-wide visualizations
========================================================================================
*/

include { FASTAN } from '../modules/fastan'
include { TANBED } from '../modules/tanbed'
include { EXTRACT_MONOMERS } from '../modules/extract_monomers'
include { CLASSIFY_MONOMERS } from '../modules/classify_monomers'
include { DETECT_HORS_REFINED as DETECT_HORS } from '../modules/detect_hors_refined'
include { CHROMOSOME_STATS } from '../modules/chromosome_stats'
include { SATELLITE_PLOTS } from '../modules/satellite_plots'
include { HOR_PLOTS } from '../modules/hor_plots'

workflow GENOME_MODE {
    take:
    genome_fasta      // Channel: genome FASTA file
    reference_monomers // Channel: reference monomer FASTA
    family_assignments // Channel: family assignment TSV

    main:

    log.info """
    ========================================
     Running GENOME MODE
    ========================================
    Analyzing genome-wide satellite arrays
    and Higher-Order Repeats (HORs)
    ========================================
    """

    //
    // STEP 1: Run FasTAN to detect tandem repeats
    //
    FASTAN(
        genome_fasta
    )

    //
    // STEP 2: Convert .1aln to BED format
    //
    TANBED(
        FASTAN.out.aln
    )

    //
    // STEP 3: Extract individual monomers from tandem arrays
    //
    EXTRACT_MONOMERS(
        genome_fasta,
        TANBED.out.bed,
        'genome'  // mode
    )

    //
    // STEP 4: Classify monomers using minimap2 + family assignments
    //
    CLASSIFY_MONOMERS(
        EXTRACT_MONOMERS.out.monomers_fasta,
        reference_monomers,
        family_assignments,
        EXTRACT_MONOMERS.out.monomer_info
    )

    //
    // STEP 5: Detect Higher-Order Repeats (HORs)
    //
    DETECT_HORS(
        CLASSIFY_MONOMERS.out.classifications,
        EXTRACT_MONOMERS.out.monomer_info
    )

    //
    // STEP 6: Generate chromosome-level statistics
    //
    CHROMOSOME_STATS(
        CLASSIFY_MONOMERS.out.classifications,
        DETECT_HORS.out.hors
    )

    //
    // STEP 7: Generate satellite/monomer-level visualizations
    //
    SATELLITE_PLOTS(
        CLASSIFY_MONOMERS.out.classifications,
        DETECT_HORS.out.hors,
        CHROMOSOME_STATS.out.stats
    )

    //
    // STEP 8: Generate HOR-specific visualizations
    //
    HOR_PLOTS(
        DETECT_HORS.out.hors,
        DETECT_HORS.out.large_duplications,
        CHROMOSOME_STATS.out.stats
    )

    emit:
    fastan_bed        = TANBED.out.bed
    monomer_fasta     = EXTRACT_MONOMERS.out.monomers_fasta
    classifications   = CLASSIFY_MONOMERS.out.classifications
    hors              = DETECT_HORS.out.hors
    large_duplications = DETECT_HORS.out.large_duplications
    chromosome_stats  = CHROMOSOME_STATS.out.stats
    satellite_plots   = SATELLITE_PLOTS.out.plots
    hor_plots         = HOR_PLOTS.out.plots
}
