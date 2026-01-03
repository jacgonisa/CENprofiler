/*
========================================================================================
    ANALYZE_DELETION_MONOMERS: Extract and classify monomers from deletions
========================================================================================
    For large deletions, the sequence exists in the reference but not in the read.
    This process extracts reference sequences from deletion regions and classifies
    the monomers using FASTAN + minimap2.

    Critical for understanding which satellite families are being lost!
========================================================================================
*/

process ANALYZE_DELETION_MONOMERS {
    tag "${sample_name}"
    publishDir "${params.outdir}/03_deletion_monomers", mode: 'copy'

    input:
    path bam_file
    path bai_file
    path reference_genome
    path reference_fai
    path indel_catalog
    path reference_monomers
    path family_assignments
    val sample_name

    output:
    path "deletion_monomers_*.tsv", emit: deletion_monomers, optional: true
    path "deletion_analysis.log", emit: log
    path "temp_deletions/*", optional: true

    when:
    params.alignment != null && params.analyze_deletions  // Run in read mode with alignment

    script:
    """
    #!/bin/bash
    set -euo pipefail

    echo "=== Deletion Monomer Analysis ===" > deletion_analysis.log
    echo "Sample: ${sample_name}" >> deletion_analysis.log
    echo "BAM: ${bam_file}" >> deletion_analysis.log
    echo "Reference: ${reference_genome}" >> deletion_analysis.log

    # Create temp directory
    mkdir -p temp_deletions

    # Count reads with large deletions
    n_reads_with_dels=\$(awk -F'\t' '\$5=="Deletion" && \$6>=100' ${indel_catalog} | cut -f1 | sort -u | wc -l)
    echo "Reads with large deletions (≥100bp): \$n_reads_with_dels" >> deletion_analysis.log

    if [ \$n_reads_with_dels -eq 0 ]; then
        echo "No large deletions found - skipping analysis" >> deletion_analysis.log
        exit 0
    fi

    # Get top reads with deletions (limit to avoid excessive runtime)
    awk -F'\t' '\$5=="Deletion" && \$6>=100' ${indel_catalog} | \\
        cut -f1 | sort | uniq -c | sort -rn | head -20 | \\
        awk '{print \$2}' > top_deletion_reads.txt

    n_analyze=\$(wc -l < top_deletion_reads.txt)
    echo "Analyzing top \$n_analyze reads with most deletions" >> deletion_analysis.log

    # Analyze each read's deletions
    # Disable exit-on-error for the entire loop to ensure all reads are processed
    set +e
    analyzed=0
    while IFS= read -r read_id; do
        echo "  Processing \$read_id..." >> deletion_analysis.log

        # Run analysis for this read
        python3 ${projectDir}/bin/analyze_deletion_monomers.py \\
            --bam ${bam_file} \\
            --ref-fasta ${reference_genome} \\
            --read-id "\$read_id" \\
            --ref-monomers ${reference_monomers} \\
            --cluster-file ${family_assignments} \\
            --output "deletion_monomers_\${read_id:0:8}.tsv" \\
            2>&1 | tee -a deletion_analysis.log

        if [ -f "deletion_monomers_\${read_id:0:8}.tsv" ]; then
            analyzed=\$((analyzed + 1))
            echo "    ✓ Analyzed \$read_id" >> deletion_analysis.log
        else
            echo "    ✗ No CEN178 monomers found for \$read_id" >> deletion_analysis.log
        fi
    done < top_deletion_reads.txt
    set -e

    echo "" >> deletion_analysis.log
    echo "Successfully analyzed deletions in \$analyzed reads" >> deletion_analysis.log

    # Combine all deletion monomer files if any exist
    if ls deletion_monomers_*.tsv 1> /dev/null 2>&1; then
        head -1 \$(ls deletion_monomers_*.tsv | head -1) > all_deletion_monomers.tsv
        for f in deletion_monomers_*.tsv; do
            tail -n +2 "\$f" >> all_deletion_monomers.tsv
        done

        n_del_monomers=\$(tail -n +2 all_deletion_monomers.tsv | wc -l)
        echo "Total deletion monomers classified: \$n_del_monomers" >> deletion_analysis.log
    else
        echo "No CEN178 monomers found in deletions" >> deletion_analysis.log
    fi

    echo "✅ Deletion monomer analysis complete" >> deletion_analysis.log
    """
}
