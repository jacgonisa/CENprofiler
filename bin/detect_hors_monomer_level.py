#!/usr/bin/env python3
"""
HOR detection at MONOMER level (not RLE level)

This correctly detects patterns like:
- F3-F3-F3 repeated 356 times → 3-homHOR
- F4-F5-F7 repeated 10 times → 3-hetHOR

Works on the original monomer family sequence before RLE compression.
"""

import pandas as pd
import numpy as np
from collections import defaultdict
import json

def find_repeating_patterns_monomer_level(family_sequence, monomers_df, min_pattern_length=3,
                                          max_pattern_length=20, min_copies=3, max_gap=500):
    """
    Find repeating patterns in monomer family sequence.
    IMPORTANT: Checks for gaps - HORs must be consecutive with no large gaps!

    Args:
        family_sequence: List of families [3, 3, 3, 4, 5, 3, 3, 3, 4, 5, ...]
        monomers_df: DataFrame with monomer positions (to check for gaps)
        min_pattern_length: Minimum monomers in pattern (≥3 for exact matching)
        max_pattern_length: Maximum monomers in pattern
        min_copies: Minimum repetitions (≥3 for exact matching)
        max_gap: Maximum allowed gap between consecutive monomers (bp)

    Returns:
        List of detected HORs
    """
    if len(family_sequence) < min_pattern_length * min_copies:
        return []

    detected_hors = []

    # Try different pattern lengths
    for pattern_len in range(min_pattern_length,
                             min(max_pattern_length + 1, len(family_sequence) // min_copies + 1)):

        # Slide window through sequence
        for start_pos in range(len(family_sequence) - pattern_len * min_copies + 1):

            # Extract candidate pattern
            pattern = tuple(family_sequence[start_pos:start_pos + pattern_len])

            # Count consecutive copies (checking for gaps!)
            copies = 1
            current_pos = start_pos + pattern_len

            while current_pos + pattern_len <= len(family_sequence):
                # Check pattern match
                next_segment = tuple(family_sequence[current_pos:current_pos + pattern_len])

                if next_segment == pattern:
                    # Pattern matches - but check for gaps!
                    # Check gap between last monomer of previous copy and first of this copy
                    prev_mono_idx = current_pos - 1
                    curr_mono_idx = current_pos

                    prev_mono_end = monomers_df.iloc[prev_mono_idx]['monomer_end']
                    curr_mono_start = monomers_df.iloc[curr_mono_idx]['monomer_start']
                    gap = curr_mono_start - prev_mono_end

                    if gap > max_gap:
                        # Gap too large - BREAK the HOR!
                        break

                    copies += 1
                    current_pos += pattern_len
                else:
                    break

            # Check if we found a valid HOR
            if copies >= min_copies:
                # Final verification: check for gaps within this HOR
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

    return detected_hors

def curate_overlapping_hors(hors):
    """
    Resolve overlapping HORs by preferring:
    1. SHORTER pattern length (simplest repeating unit)
    2. More copies (if pattern length equal)

    This prefers detecting 3F3 × 356 over 20F3 × 53
    """
    if len(hors) == 0:
        return []

    # Sort by start position
    hors_sorted = sorted(hors, key=lambda x: x['start_monomer_idx'])

    curated = []
    covered = set()

    # Group overlapping HORs
    for hor in hors_sorted:
        hor_range = set(range(hor['start_monomer_idx'], hor['end_monomer_idx']))

        # Check if this HOR overlaps with already covered regions
        overlap = hor_range & covered

        if len(overlap) == 0:
            # No overlap - add it
            curated.append(hor)
            covered.update(hor_range)
        else:
            # Overlaps - need to decide which one to keep
            # Find overlapping HORs in curated list
            overlapping_hors = []
            for i, existing in enumerate(curated):
                existing_range = set(range(existing['start_monomer_idx'],
                                          existing['end_monomer_idx']))
                if existing_range & hor_range:
                    overlapping_hors.append((i, existing))

            # Compare with best overlapping HOR
            # Prefer SHORTER patterns (simpler units), then more copies
            if overlapping_hors:
                best_idx, best_hor = min(overlapping_hors,
                                        key=lambda x: (x[1]['pattern_length'],
                                                      -x[1]['copies']))

                # Should we replace it?
                if (hor['pattern_length'] < best_hor['pattern_length'] or
                    (hor['pattern_length'] == best_hor['pattern_length'] and
                     hor['copies'] > best_hor['copies'])):
                    # Remove old ones and add new one
                    for idx, _ in sorted(overlapping_hors, reverse=True):
                        old_hor = curated.pop(idx)
                        old_range = set(range(old_hor['start_monomer_idx'],
                                             old_hor['end_monomer_idx']))
                        covered -= old_range

                    curated.append(hor)
                    covered.update(hor_range)

    return curated

def analyze_centromere_array(monomers_df, min_pattern_length=3, max_pattern_length=20,
                             min_copies=3):
    """
    Analyze a centromere array for HORs at monomer level.

    Args:
        monomers_df: DataFrame with columns: monomer_family, monomer_start, monomer_end
        min_pattern_length: Minimum monomers in pattern
        max_pattern_length: Maximum monomers in pattern
        min_copies: Minimum repetitions

    Returns:
        DataFrame of detected HORs
    """
    # Filter out NaN families and extract sequence
    monomers_valid = monomers_df[monomers_df['monomer_family'].notna()].copy()
    monomers_valid = monomers_valid.reset_index(drop=True)  # Reset index for proper iloc access
    family_sequence = monomers_valid['monomer_family'].astype(int).tolist()

    # Detect HORs (with gap checking!)
    raw_hors = find_repeating_patterns_monomer_level(
        family_sequence,
        monomers_valid,  # Pass monomers DataFrame for gap checking
        min_pattern_length=min_pattern_length,
        max_pattern_length=max_pattern_length,
        min_copies=min_copies
    )

    # Curate overlaps
    curated_hors = curate_overlapping_hors(raw_hors)

    # Convert to DataFrame with genomic coordinates
    hor_records = []
    for hor in curated_hors:
        # Get genomic coordinates from monomer indices
        start_mono = monomers_valid.iloc[hor['start_monomer_idx']]
        end_mono = monomers_valid.iloc[hor['end_monomer_idx'] - 1]

        hor_start = start_mono['monomer_start']
        hor_end = end_mono['monomer_end']

        # Build HOR unit string (like "3F3" or "1F4-1F5-1F7")
        pattern = hor['pattern']

        # Check if homHOR or hetHOR
        if len(set(pattern)) == 1:
            # homHOR
            hor_unit = f"{len(pattern)}F{pattern[0]}"
            hor_type = 'homHOR'
        else:
            # hetHOR - build RLE representation
            elements = []
            current_fam = pattern[0]
            current_count = 1

            for fam in pattern[1:]:
                if fam == current_fam:
                    current_count += 1
                else:
                    elements.append(f"{current_count}F{current_fam}")
                    current_fam = fam
                    current_count = 1
            elements.append(f"{current_count}F{current_fam}")

            hor_unit = '-'.join(elements)
            hor_type = 'hetHOR'

        hor_records.append({
            'hor_start': hor_start,
            'hor_end': hor_end,
            'hor_unit': hor_unit,
            'hor_unit_length': hor['pattern_length'],
            'hor_copies': hor['copies'],
            'total_monomers': hor['total_monomers'],
            'hor_type': hor_type,
            'pattern_tuple': str(pattern),
            'start_monomer_idx': hor['start_monomer_idx'],
            'end_monomer_idx': hor['end_monomer_idx']
        })

    return pd.DataFrame(hor_records)

# Test on a known case
if __name__ == '__main__':
    print("=== TESTING MONOMER-LEVEL HOR DETECTION (WITH GAP CHECKING) ===\n")

    # Test case 1: 1068 consecutive F3 monomers (NO GAPS)
    print("Test 1: 1068 consecutive F3 monomers (NO GAPS)")
    print("Expected: Should detect 3F3 pattern\n")

    test_monomers = pd.DataFrame({
        'monomer_family': [3] * 1068,
        'monomer_start': [i * 178 for i in range(1068)],
        'monomer_end': [(i + 1) * 178 for i in range(1068)]
    })

    hors_detected = analyze_centromere_array(test_monomers,
                                             min_pattern_length=3,
                                             max_pattern_length=20,
                                             min_copies=3)

    print(f"Detected {len(hors_detected)} HORs:")
    if len(hors_detected) > 0:
        print(hors_detected[['hor_unit', 'hor_copies', 'total_monomers',
                            'hor_type']].to_string(index=False))

    print("\n" + "="*60 + "\n")

    # Test case 2: F3 monomers WITH A LARGE GAP in the middle
    print("Test 2: 500 F3 + GAP + 500 F3")
    print("Expected: Should detect TWO separate 3F3 HORs (gap breaks them)\n")

    test_monomers_gap = pd.DataFrame({
        'monomer_family': [3] * 500 + [3] * 500,
        'monomer_start': [i * 178 for i in range(500)] +
                        [i * 178 + 100000 for i in range(500, 1000)],  # 100kb gap!
        'monomer_end': [(i + 1) * 178 for i in range(500)] +
                      [(i + 1) * 178 + 100000 for i in range(500, 1000)]
    })

    hors_detected_gap = analyze_centromere_array(test_monomers_gap,
                                                 min_pattern_length=3,
                                                 max_pattern_length=20,
                                                 min_copies=3)

    print(f"Detected {len(hors_detected_gap)} HORs:")
    if len(hors_detected_gap) > 0:
        print(hors_detected_gap[['hor_unit', 'hor_copies', 'total_monomers',
                                'hor_type', 'hor_start', 'hor_end']].to_string(index=False))

    print("\n" + "="*60 + "\n")

    # Test case 3: Repeating hetHOR pattern
    print("Test 3: F4-F5-F7 repeated 10 times")
    print("Expected: Should detect 1F4-1F5-1F7 × 10 copies\n")

    test_pattern = [4, 5, 7] * 10
    test_monomers2 = pd.DataFrame({
        'monomer_family': test_pattern,
        'monomer_start': [i * 178 for i in range(len(test_pattern))],
        'monomer_end': [(i + 1) * 178 for i in range(len(test_pattern))]
    })

    hors_detected2 = analyze_centromere_array(test_monomers2,
                                              min_pattern_length=3,
                                              max_pattern_length=10,
                                              min_copies=3)

    print(f"Detected {len(hors_detected2)} HORs:")
    if len(hors_detected2) > 0:
        print(hors_detected2[['hor_unit', 'hor_copies', 'total_monomers',
                             'hor_type']].to_string(index=False))

    print("\n✅ Monomer-level HOR detection ready (with gap checking)!")
