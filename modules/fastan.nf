/*
========================================================================================
    FASTAN: Tandem Repeat Detection
========================================================================================
    Detects tandem repeats using FasTAN
    https://github.com/thegenemyers/FASTAN
========================================================================================
*/

process FASTAN {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/01_fastan", mode: 'copy'

    input:
    path fasta

    output:
    path "*.1aln", emit: aln
    path "fastan.log", emit: log

    script:
    def prefix = fasta.baseName
    """
    # Run FasTAN for tandem repeat detection
    ${params.fastan_bin} \\
        -v \\
        -T${task.cpus} \\
        ${fasta} \\
        ${prefix}.1aln \\
        2>&1 | tee fastan.log

    # Verify output was created
    if [ ! -f "${prefix}.1aln" ]; then
        echo "ERROR: FasTAN did not produce output file"
        exit 1
    fi

    echo "FasTAN completed successfully" >> fastan.log
    echo "Output: ${prefix}.1aln" >> fastan.log
    """
}
