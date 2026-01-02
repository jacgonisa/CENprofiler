/*
========================================================================================
    LOAD_GENOMIC_REGIONS: Load genomic region annotations
========================================================================================
    Loads centromere, pericentromere, and rDNA regions from BED files.
    Creates unified regions file for indel classification.
========================================================================================
*/

process LOAD_GENOMIC_REGIONS {
    tag "genomic regions"
    publishDir "${params.outdir}/00_regions", mode: 'copy'

    input:
    path annotation_dir

    output:
    path "genomic_regions.tsv", emit: regions
    path "regions.log", emit: log

    script:
    """
    python3 ${projectDir}/bin/load_genomic_regions.py \\
        ${annotation_dir} \\
        --output genomic_regions.tsv \\
        2> regions.log

    echo "✅ Loaded genomic regions" >> regions.log
    """
}
