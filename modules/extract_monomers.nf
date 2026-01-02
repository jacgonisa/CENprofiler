/*
========================================================================================
    EXTRACT_MONOMERS: Extract individual monomers from tandem arrays
========================================================================================
    Extracts individual repeat units based on FasTAN detected arrays
========================================================================================
*/

process EXTRACT_MONOMERS {
    tag "${fasta.baseName}"
    publishDir "${params.outdir}/02_monomers", mode: 'copy'

    input:
    path fasta
    path bed
    val mode  // 'genome' or 'reads'

    output:
    path "monomers.fa", emit: monomers_fasta
    path "monomer_info.tsv", emit: monomer_info
    path "extraction.log", emit: log

    script:
    """
    #!/usr/bin/env python3
    import sys
    from Bio import SeqIO
    import pandas as pd

    # Load sequences
    print("Loading sequences...", file=sys.stderr)
    sequences = {}
    for record in SeqIO.parse("${fasta}", "fasta"):
        seq_id = record.id.split()[0]
        sequences[seq_id] = str(record.seq).upper()
    print(f"Loaded {len(sequences)} sequences", file=sys.stderr)

    # Load tandem arrays from BED
    print("Loading tandem arrays...", file=sys.stderr)
    arrays = pd.read_csv("${bed}", sep='\\t', header=None,
                         names=['seq_id', 'start', 'end', 'period', 'quality'])

    # Filter for CEN178-like period
    arrays = arrays[
        (arrays['period'] >= ${params.period_min}) &
        (arrays['period'] <= ${params.period_max})
    ]
    print(f"Found {len(arrays)} CEN178-like arrays (period ${params.period_min}-${params.period_max}bp)",
          file=sys.stderr)

    # Extract monomers
    monomer_records = []
    monomer_info = []
    monomer_count = 0

    for idx, row in arrays.iterrows():
        seq_id = row['seq_id']
        start = int(row['start'])
        end = int(row['end'])
        period = int(row['period'])
        quality = row['quality']

        if seq_id not in sequences:
            print(f"WARNING: Sequence {seq_id} not found", file=sys.stderr)
            continue

        # Extract array sequence
        array_seq = sequences[seq_id][start:end]
        array_length = end - start

        # Extract individual monomers
        pos = 0
        monomer_idx = 0

        while pos + period <= array_length:
            monomer_seq = array_seq[pos:pos + period]
            monomer_start = start + pos
            monomer_end = start + pos + period
            monomer_length = period

            # Filter by length
            if monomer_length < ${params.min_monomer_length}:
                pos += period
                continue
            if monomer_length > ${params.max_monomer_length}:
                pos += period
                continue

            # Create monomer ID
            monomer_id = f"{seq_id}_array{idx}_mon{monomer_idx}"

            # Store monomer
            monomer_records.append(f">{monomer_id}\\n{monomer_seq}")

            # Store info
            monomer_info.append({
                'monomer_id': monomer_id,
                'seq_id': seq_id,
                'array_idx': idx,
                'monomer_idx': monomer_idx,
                'monomer_start': monomer_start,
                'monomer_end': monomer_end,
                'monomer_length': monomer_length,
                'array_period': period,
                'array_quality': quality
            })

            monomer_idx += 1
            monomer_count += 1
            pos += period

    # Write monomers FASTA
    print(f"Writing {monomer_count} monomers to monomers.fa", file=sys.stderr)
    with open("monomers.fa", 'w') as f:
        f.write('\\n'.join(monomer_records) + '\\n')

    # Write monomer info
    df_info = pd.DataFrame(monomer_info)
    df_info.to_csv("monomer_info.tsv", sep='\\t', index=False)

    # Write log
    with open("extraction.log", 'w') as f:
        f.write(f"Extraction completed\\n")
        f.write(f"Mode: ${mode}\\n")
        f.write(f"Input sequences: {len(sequences)}\\n")
        f.write(f"Tandem arrays: {len(arrays)}\\n")
        f.write(f"Monomers extracted: {monomer_count}\\n")
        f.write(f"Period range: ${params.period_min}-${params.period_max}bp\\n")

    print("Extraction completed successfully", file=sys.stderr)
    """
}
