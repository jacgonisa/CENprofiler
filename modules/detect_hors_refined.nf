/*
========================================================================================
    DETECT_HORS_REFINED: Detect HORs with Quality Metrics
========================================================================================
    Refined gap-aware monomer-level HOR detection with quality scoring

    Improvements over standard detection:
    - Purity scoring (0-1): Measures pattern perfection
    - Quality scoring (0-100): Overall HOR quality
    - Gap metrics: max, mean, std deviation
    - Better overlap resolution: Prefers high-quality HORs
    - Flexible filtering: min_purity and min_score parameters

    Requirements:
        - min_copies ≥ 3 (pattern repeats at least 3 times)
        - min_monomers ≥ 3 (pattern comprises at least 3 monomers)
        - Gaps > max_gap break HORs
        - Purity >= min_purity (default 0.9)
        - Quality >= min_score (default 50)
========================================================================================
*/

process DETECT_HORS_REFINED {
    tag "HOR detection (refined)"
    publishDir "${params.outdir}/03_hors", mode: 'copy'

    input:
    path classifications
    path monomer_info

    output:
    path "hors_detected.tsv", emit: hors
    path "large_duplications.tsv", emit: large_duplications
    path "hor_detection.log", emit: log

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import sys
    from collections import defaultdict

    # Load refined HOR detection functions
    exec(open("${projectDir}/bin/detect_hors_refined.py").read())

    # Load data
    print("Loading monomer classifications...", file=sys.stderr)
    df = pd.read_csv("${classifications}", sep='\\t')
    print(f"Loaded {len(df)} monomers", file=sys.stderr)

    # Filter for classified monomers
    df_classified = df[df['monomer_family'].notna()].copy()
    print(f"Classified: {len(df_classified)} monomers", file=sys.stderr)

    if len(df_classified) == 0:
        print("ERROR: No classified monomers found", file=sys.stderr)
        pd.DataFrame().to_csv("hors_detected.tsv", sep='\\t', index=False)
        pd.DataFrame().to_csv("large_duplications.tsv", sep='\\t', index=False)
        with open("hor_detection.log", 'w') as f:
            f.write("No classified monomers - no HORs detected\\n")
        sys.exit(0)

    # Group by sequence/array
    df_classified['array_id'] = df_classified['seq_id'] + '_array' + df_classified['array_idx'].astype(str)

    # Detect HORs per array using refined algorithm
    all_hors = []

    for array_id, group in df_classified.groupby('array_id'):
        group = group.sort_values('monomer_idx').reset_index(drop=True)

        # Use refined detection with quality metrics
        hors_df = analyze_centromere_array_refined(
            group,
            min_pattern_length=${params.min_monomers},
            max_pattern_length=${params.max_pattern_length},
            min_copies=${params.min_copies},
            max_gap=${params.max_gap},
            min_purity=${params.hor_min_purity ?: 0.9},
            min_score=${params.hor_min_score ?: 50}
        )

        if len(hors_df) > 0:
            # Add seq_id and array_idx
            hors_df['seq_id'] = group.iloc[0]['seq_id']
            hors_df['array_idx'] = group.iloc[0]['array_idx']

            # Reorder columns
            hors_df = hors_df[['seq_id', 'array_idx', 'hor_start', 'hor_end', 'hor_unit',
                              'hor_unit_length', 'hor_copies', 'total_monomers', 'hor_type',
                              'purity', 'quality_score', 'max_gap', 'mean_gap', 'gap_std']]

            all_hors.append(hors_df)

    # Combine all HORs
    if all_hors:
        df_hors = pd.concat(all_hors, ignore_index=True)

        # Calculate bp and kb lengths
        df_hors['hor_length_bp'] = df_hors['hor_end'] - df_hors['hor_start']
        df_hors['hor_length_kb'] = df_hors['hor_length_bp'] / 1000

        # Find large duplications
        large_threshold = ${params.large_dup_threshold} * 1000  # Convert kb to bp
        df_large = df_hors[df_hors['hor_length_bp'] >= large_threshold].copy()
        df_large = df_large.sort_values('hor_length_bp', ascending=False)
    else:
        df_hors = pd.DataFrame()
        df_large = pd.DataFrame()

    # Save results
    df_hors.to_csv("hors_detected.tsv", sep='\\t', index=False)
    df_large.to_csv("large_duplications.tsv", sep='\\t', index=False)

    # Write log
    with open("hor_detection.log", 'w') as f:
        f.write("Refined HOR Detection Completed\\n")
        f.write(f"Total HORs detected: {len(df_hors)}\\n")

        if len(df_hors) > 0:
            n_homhor = (df_hors['hor_type'] == 'homHOR').sum()
            n_hethor = (df_hors['hor_type'] == 'hetHOR').sum()
            f.write(f"homHORs: {n_homhor}\\n")
            f.write(f"hetHORs: {n_hethor}\\n")
            f.write(f"Large duplications (≥${params.large_dup_threshold}kb): {len(df_large)}\\n")
            f.write(f"\\nQuality Metrics:\\n")
            f.write(f"  Mean purity: {df_hors['purity'].mean():.3f}\\n")
            f.write(f"  Mean quality score: {df_hors['quality_score'].mean():.1f}\\n")
            f.write(f"  Mean max gap: {df_hors['max_gap'].mean():.1f} bp\\n")

            if len(df_large) > 0:
                f.write(f"\\nLargest HOR: {df_large.iloc[0]['hor_length_kb']:.1f} kb\\n")
                f.write(f"  Unit: {df_large.iloc[0]['hor_unit']}\\n")
                f.write(f"  Quality: {df_large.iloc[0]['quality_score']:.1f}\\n")

    print(f"Detected {len(df_hors)} HORs total", file=sys.stderr)
    if len(df_hors) > 0:
        print(f"Mean quality score: {df_hors['quality_score'].mean():.1f}", file=sys.stderr)
    """
}
