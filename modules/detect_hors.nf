/*
========================================================================================
    DETECT_HORS: Detect Higher-Order Repeats (HORs)
========================================================================================
    Gap-aware monomer-level HOR detection
    Requirements:
        - min_copies ≥ 3 (pattern repeats at least 3 times)
        - min_monomers ≥ 3 (pattern comprises at least 3 monomers)
        - Gaps > max_gap break HORs
========================================================================================
*/

process DETECT_HORS {
    tag "HOR detection"
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

    # Load data
    print("Loading monomer classifications...", file=sys.stderr)
    df = pd.read_csv("${classifications}", sep='\\t')
    print(f"Loaded {len(df)} monomers", file=sys.stderr)

    # Filter for classified monomers
    df_classified = df[df['monomer_family'].notna()].copy()
    print(f"Classified: {len(df_classified)} monomers", file=sys.stderr)

    if len(df_classified) == 0:
        print("ERROR: No classified monomers found", file=sys.stderr)
        # Create empty output files
        pd.DataFrame().to_csv("hors_detected.tsv", sep='\\t', index=False)
        pd.DataFrame().to_csv("large_duplications.tsv", sep='\\t', index=False)
        with open("hor_detection.log", 'w') as f:
            f.write("No classified monomers - no HORs detected\\n")
        sys.exit(0)

    # Group by sequence/array
    df_classified['array_id'] = df_classified['seq_id'] + '_array' + df_classified['array_idx'].astype(str)

    # Gap-aware HOR detection function
    def find_repeating_patterns_monomer_level(family_sequence, monomers_df,
                                             min_pattern_length=${params.min_monomers},
                                             max_pattern_length=${params.max_pattern_length},
                                             min_copies=${params.min_copies},
                                             max_gap=${params.max_gap}):
        """
        Find repeating patterns in monomer family sequence with gap checking
        """
        if len(family_sequence) < min_pattern_length * min_copies:
            return []

        detected_hors = []
        occupied = set()

        # Try different pattern lengths (prefer shorter)
        for pattern_len in range(min_pattern_length,
                                 min(max_pattern_length + 1, len(family_sequence) // min_copies + 1)):

            for start_pos in range(len(family_sequence) - pattern_len * min_copies + 1):
                # Skip if already part of a HOR
                if any(i in occupied for i in range(start_pos, start_pos + pattern_len)):
                    continue

                # Extract candidate pattern
                pattern = tuple(family_sequence[start_pos:start_pos + pattern_len])

                # Count consecutive copies with gap checking
                copies = 1
                current_pos = start_pos + pattern_len

                while current_pos + pattern_len <= len(family_sequence):
                    next_segment = tuple(family_sequence[current_pos:current_pos + pattern_len])

                    if next_segment == pattern:
                        # Check for gaps
                        prev_mono_idx = current_pos - 1
                        curr_mono_idx = current_pos

                        prev_mono_end = monomers_df.iloc[prev_mono_idx]['monomer_end']
                        curr_mono_start = monomers_df.iloc[curr_mono_idx]['monomer_start']
                        gap = curr_mono_start - prev_mono_end

                        if gap > max_gap:
                            break  # Gap too large

                        copies += 1
                        current_pos += pattern_len
                    else:
                        break

                # Check if valid HOR
                if copies >= min_copies:
                    # Final gap verification
                    has_large_gap = False
                    for i in range(start_pos, start_pos + (copies * pattern_len) - 1):
                        mono_end = monomers_df.iloc[i]['monomer_end']
                        next_mono_start = monomers_df.iloc[i + 1]['monomer_start']
                        gap = next_mono_start - mono_end

                        if gap > max_gap:
                            has_large_gap = True
                            break

                    if not has_large_gap:
                        hor_end = start_pos + (copies * pattern_len)

                        detected_hors.append({
                            'pattern': pattern,
                            'pattern_length': pattern_len,
                            'copies': copies,
                            'start_monomer_idx': start_pos,
                            'end_monomer_idx': hor_end,
                            'total_monomers': copies * pattern_len
                        })

                        # Mark positions as occupied
                        for i in range(start_pos, hor_end):
                            occupied.add(i)

        return detected_hors

    # Detect HORs per array
    all_hors = []

    for array_id, group in df_classified.groupby('array_id'):
        group = group.sort_values('monomer_idx').reset_index(drop=True)

        family_sequence = group['monomer_family'].astype(int).tolist()

        hors = find_repeating_patterns_monomer_level(family_sequence, group)

        for hor in hors:
            # Get genomic coordinates
            start_idx = hor['start_monomer_idx']
            end_idx = hor['end_monomer_idx']

            hor_start = group.iloc[start_idx]['monomer_start']
            hor_end = group.iloc[end_idx - 1]['monomer_end']
            hor_length_bp = hor_end - hor_start

            # Format HOR unit
            pattern_str = '-'.join([f"1F{fam}" for fam in hor['pattern']])

            # Classify as homHOR or hetHOR
            hor_type = 'homHOR' if len(set(hor['pattern'])) == 1 else 'hetHOR'

            all_hors.append({
                'seq_id': group.iloc[0]['seq_id'],
                'array_idx': group.iloc[0]['array_idx'],
                'hor_start': hor_start,
                'hor_end': hor_end,
                'hor_unit': pattern_str,
                'hor_unit_length': hor['pattern_length'],
                'hor_copies': hor['copies'],
                'total_monomers': hor['total_monomers'],
                'hor_type': hor_type,
                'hor_length_bp': hor_length_bp
            })

    # Convert to DataFrame
    df_hors = pd.DataFrame(all_hors)

    # Calculate kb
    if len(df_hors) > 0:
        df_hors['hor_length_kb'] = df_hors['hor_length_bp'] / 1000

        # Find large duplications
        large_threshold = ${params.large_dup_threshold} * 1000  # Convert kb to bp
        df_large = df_hors[df_hors['hor_length_bp'] >= large_threshold].copy()
        df_large = df_large.sort_values('hor_length_bp', ascending=False)
    else:
        df_hors['hor_length_kb'] = []
        df_large = pd.DataFrame()

    # Save results
    df_hors.to_csv("hors_detected.tsv", sep='\\t', index=False)
    df_large.to_csv("large_duplications.tsv", sep='\\t', index=False)

    # Write log
    with open("hor_detection.log", 'w') as f:
        f.write("HOR Detection Completed\\n")
        f.write(f"Total HORs detected: {len(df_hors)}\\n")

        if len(df_hors) > 0:
            n_homhor = (df_hors['hor_type'] == 'homHOR').sum()
            n_hethor = (df_hors['hor_type'] == 'hetHOR').sum()
            f.write(f"homHORs: {n_homhor}\\n")
            f.write(f"hetHORs: {n_hethor}\\n")
            f.write(f"Large duplications (≥${params.large_dup_threshold}kb): {len(df_large)}\\n")

            if len(df_large) > 0:
                f.write(f"\\nLargest HOR: {df_large.iloc[0]['hor_length_kb']:.1f} kb\\n")

    print(f"Detected {len(df_hors)} HORs total", file=sys.stderr)
    """
}
