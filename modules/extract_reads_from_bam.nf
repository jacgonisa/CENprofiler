/*
========================================================================================
    EXTRACT_READS_FROM_BAM: Extract reads with large indels
========================================================================================
    Extracts reads containing insertions or deletions ≥ min_indel_size.
    Classifies indels by genomic region.
    Outputs: reads FASTA, indel catalog, statistics
========================================================================================
*/

process EXTRACT_READS_FROM_BAM {
    tag "${sample_name}"
    publishDir "${params.outdir}/01_extracted_reads", mode: 'copy'

    input:
    path bam_file
    path bam_index
    path regions_file
    val sample_name

    output:
    path "${sample_name}_reads.fa", emit: reads_fasta
    path "${sample_name}_indel_catalog.tsv", emit: indel_catalog
    path "${sample_name}_stats.txt", emit: stats
    path "extraction.log", emit: log

    script:
    """
    python3 ${projectDir}/bin/extract_reads_from_bam.py \\
        ${bam_file} \\
        ${regions_file} \\
        ${params.min_indel_size} \\
        ${sample_name} \\
        2> extraction.log

    echo "✅ Extracted reads with indels ≥${params.min_indel_size}bp" >> extraction.log
    """
}
