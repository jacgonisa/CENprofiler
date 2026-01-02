/*
========================================================================================
    READ MODE WORKFLOW
========================================================================================
    Analyzes long reads for:
    - Tandem satellite arrays (FasTAN)
    - Monomer classification per read
    - Family transitions within arrays
    - Indel/variant analysis
    - Per-read visualizations
========================================================================================
*/

include { FASTAN } from '../modules/fastan'
include { TANBED } from '../modules/tanbed'
include { EXTRACT_MONOMERS } from '../modules/extract_monomers'
include { CLASSIFY_MONOMERS } from '../modules/classify_monomers'
include { ANALYZE_INDELS } from '../modules/analyze_indels'
include { READ_PLOTS } from '../modules/read_plots'

workflow READ_MODE {
    take:
    reads_fasta        // Channel: reads FASTA file
    reference_monomers // Channel: reference monomer FASTA
    family_assignments // Channel: family assignment TSV

    main:

    log.info """
    ========================================
     Running READ MODE
    ========================================
    Analyzing satellite composition and
    variants in long reads
    ========================================
    """

    //
    // STEP 1: Run FasTAN to detect tandem repeats
    //
    FASTAN(
        reads_fasta
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
        reads_fasta,
        TANBED.out.bed,
        'reads'  // mode
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
    // STEP 5: Analyze indels and structural variants (optional)
    //
    if (params.analyze_indels) {
        ANALYZE_INDELS(
            CLASSIFY_MONOMERS.out.classifications,
            EXTRACT_MONOMERS.out.monomer_info
        )
    }

    //
    // STEP 6: Generate per-read visualizations
    //
    READ_PLOTS(
        CLASSIFY_MONOMERS.out.classifications,
        EXTRACT_MONOMERS.out.monomer_info,
        params.analyze_indels ? ANALYZE_INDELS.out.indel_stats : Channel.empty()
    )

    emit:
    fastan_bed        = TANBED.out.bed
    monomer_fasta     = EXTRACT_MONOMERS.out.monomers_fasta
    classifications   = CLASSIFY_MONOMERS.out.classifications
    indel_stats       = params.analyze_indels ? ANALYZE_INDELS.out.indel_stats : Channel.empty()
    plots             = READ_PLOTS.out.plots
}
