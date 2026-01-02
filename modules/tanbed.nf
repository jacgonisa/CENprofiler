/*
========================================================================================
    TANBED: Convert FasTAN .1aln to BED format
========================================================================================
    Uses alntools tanbed to convert FasTAN output to BED format
    https://github.com/richarddurbin/alntools
========================================================================================
*/

process TANBED {
    tag "${aln.baseName}"
    publishDir "${params.outdir}/01_fastan", mode: 'copy'

    input:
    path aln

    output:
    path "*.bed", emit: bed
    path "tanbed.log", emit: log

    script:
    def prefix = aln.baseName.replaceAll(/\.1aln$/, '')
    """
    # Convert .1aln to BED format
    ${params.tanbed_bin} \\
        ${aln} \\
        > ${prefix}.bed \\
        2> tanbed.log

    # Verify output was created
    if [ ! -f "${prefix}.bed" ]; then
        echo "ERROR: tanbed did not produce output file"
        exit 1
    fi

    # Count tandem arrays
    n_arrays=\$(wc -l < ${prefix}.bed)
    echo "tanbed completed successfully" >> tanbed.log
    echo "Detected \${n_arrays} tandem arrays" >> tanbed.log
    """
}
