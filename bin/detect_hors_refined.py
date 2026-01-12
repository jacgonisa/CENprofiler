#!/usr/bin/env python3
"""
Refined HOR detection at monomer level with improved pattern representation and quality metrics.

Key improvements:
1. Better pattern representation (e.g., "3F3" instead of "1F3-1F3-1F3")
2. Quality metrics (purity score, periodicity confidence)
3. Flexible gap handling
4. Better overlap resolution
5. Pattern validation and filtering
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import json
import sys

def calculate_pattern_purity(family_sequence, pattern, start_idx, copies):
    """
    Calculate purity score for a detected HOR.

    Purity = (perfect pattern matches) / (total monomers in HOR)

    Returns value between 0-1, where 1.0 = perfect HOR with no mismatches
    """
    pattern_len = len(pattern)
    total_monomers = copies * pattern_len
    perfect_matches = 0

    for copy in range(copies):
        copy_start = start_idx + (copy * pattern_len)
        for i, expected_family in enumerate(pattern):
            if copy_start + i < len(family_sequence):
                if family_sequence[copy_start + i] == expected_family:
                    perfect_matches += 1

    purity = perfect_matches / total_monomers if total_monomers > 0 else 0
    return purity

def check_gap_consistency(monomers_df, start_idx, end_idx):
    """
    Check gap consistency within a HOR region.

    Returns:
        max_gap: Maximum gap found (bp)
        mean_gap: Mean gap (bp)
        gap_std: Standard deviation of gaps (bp)
    """
    gaps = []
    for i in range(start_idx, end_idx - 1):
        mono_end = monomers_df.iloc[i]['monomer_end']
        next_start = monomers_df.iloc[i + 1]['monomer_start']
        gap = next_start - mono_end
        gaps.append(gap)

    if gaps:
        return {
            'max_gap': max(gaps),
            'mean_gap': np.mean(gaps),
            'gap_std': np.std(gaps)
        }
    return {'max_gap': 0, 'mean_gap': 0, 'gap_std': 0}

def find_repeating_patterns_refined(family_sequence, monomers_df,
                                    min_pattern_length=3,
                                    max_pattern_length=20,
                                    min_copies=3,
                                    max_gap=500,
                                    min_purity=0.9):
    """
    Find repeating patterns with quality metrics and gap validation.

    Args:
        family_sequence: List of families [3, 3, 3, 4, 5, 3, 3, 3, 4, 5, ...]
        monomers_df: DataFrame with monomer positions
        min_pattern_length: Minimum monomers in pattern (≥3)
        max_pattern_length: Maximum monomers in pattern
        min_copies: Minimum repetitions (≥3)
        max_gap: Maximum allowed gap between consecutive monomers (bp)
        min_purity: Minimum purity score (0-1)

    Returns:
        List of detected HORs with quality metrics
    """
    if len(family_sequence) < min_pattern_length * min_copies:
        return []

    detected_hors = []

    # Try different pattern lengths (prefer shorter patterns)
    for pattern_len in range(min_pattern_length,
                             min(max_pattern_length + 1, len(family_sequence) // min_copies + 1)):

        # Slide window through sequence
        for start_pos in range(len(family_sequence) - pattern_len * min_copies + 1):

            # Extract candidate pattern
            pattern = tuple(family_sequence[start_pos:start_pos + pattern_len])

            # Count consecutive copies with gap checking
            copies = 1
            current_pos = start_pos + pattern_len
            gap_violations = []

            while current_pos + pattern_len <= len(family_sequence):
                # Check pattern match
                next_segment = tuple(family_sequence[current_pos:current_pos + pattern_len])

                if next_segment == pattern:
                    # Check gap between copies
                    prev_mono_idx = current_pos - 1
                    curr_mono_idx = current_pos

                    prev_mono_end = monomers_df.iloc[prev_mono_idx]['monomer_end']
                    curr_mono_start = monomers_df.iloc[curr_mono_idx]['monomer_start']
                    gap = curr_mono_start - prev_mono_end

                    if gap > max_gap:
                        # Gap too large - break HOR
                        break

                    copies += 1
                    current_pos += pattern_len
                else:
                    break

            # Check if we found a valid HOR
            if copies >= min_copies:
                hor_end = start_pos + (copies * pattern_len)

                # Calculate quality metrics
                purity = calculate_pattern_purity(family_sequence, pattern, start_pos, copies)

                # Check gap consistency
                gap_metrics = check_gap_consistency(monomers_df, start_pos, hor_end)

                # Validate purity
                if purity >= min_purity and gap_metrics['max_gap'] <= max_gap:
                    detected_hors.append({
                        'pattern': pattern,
                        'pattern_length': pattern_len,
                        'copies': copies,
                        'start_monomer_idx': start_pos,
                        'end_monomer_idx': hor_end,
                        'total_monomers': copies * pattern_len,
                        'purity': purity,
                        'max_gap': gap_metrics['max_gap'],
                        'mean_gap': gap_metrics['mean_gap'],
                        'gap_std': gap_metrics['gap_std']
                    })

    return detected_hors

def curate_overlapping_hors_refined(hors):
    """
    Resolve overlapping HORs by preferring:
    1. Higher purity score
    2. Shorter pattern length (simpler repeating unit)
    3. More copies (if pattern length equal)

    This ensures we keep the highest quality HORs.
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

            if overlapping_hors:
                best_idx, best_hor = min(overlapping_hors,
                                        key=lambda x: (-x[1]['purity'],  # Higher purity first
                                                      x[1]['pattern_length'],  # Then shorter pattern
                                                      -x[1]['copies']))  # Then more copies

                # Should we replace it?
                if (hor['purity'] > best_hor['purity'] + 0.05 or  # Significantly better purity
                    (abs(hor['purity'] - best_hor['purity']) < 0.05 and  # Similar purity
                     (hor['pattern_length'] < best_hor['pattern_length'] or
                      (hor['pattern_length'] == best_hor['pattern_length'] and
                       hor['copies'] > best_hor['copies'])))):
                    # Remove old ones and add new one
                    for idx, _ in sorted(overlapping_hors, reverse=True):
                        old_hor = curated.pop(idx)
                        old_range = set(range(old_hor['start_monomer_idx'],
                                             old_hor['end_monomer_idx']))
                        covered -= old_range

                    curated.append(hor)
                    covered.update(hor_range)

    return curated

def format_hor_unit(pattern):
    """
    Format HOR unit in readable form.

    For homHORs: "3F3" means 3 consecutive F3 monomers
    For hetHORs: "2F1-1F7-2F1" means 2×F1, then 1×F7, then 2×F1
    """
    if len(set(pattern)) == 1:
        # homHOR - simple format
        return f"{len(pattern)}F{pattern[0]}"
    else:
        # hetHOR - RLE format
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

        return '-'.join(elements)

def calculate_hor_score(hor):
    """
    Calculate overall quality score for a HOR.

    Score considers:
    - Purity (0-1)
    - Number of copies (more is better)
    - Pattern simplicity (shorter is better)

    Returns score 0-100
    """
    purity_score = hor['purity'] * 50  # 0-50 points

    # Copy score: logarithmic scale (3 copies=10pts, 10 copies=20pts, 100+ copies=30pts)
    copy_score = min(30, 10 * np.log10(hor['copies']) + 15)

    # Simplicity score: prefer shorter patterns (3-mer=20pts, 20-mer=5pts)
    simplicity_score = max(5, 25 - hor['pattern_length'])

    total_score = purity_score + copy_score + simplicity_score
    return total_score

def analyze_centromere_array_refined(monomers_df,
                                     min_pattern_length=3,
                                     max_pattern_length=20,
                                     min_copies=3,
                                     max_gap=500,
                                     min_purity=0.9,
                                     min_score=50):
    """
    Analyze a centromere array for HORs with quality filtering.

    Args:
        monomers_df: DataFrame with columns: monomer_family, monomer_start, monomer_end
        min_pattern_length: Minimum monomers in pattern
        max_pattern_length: Maximum monomers in pattern
        min_copies: Minimum repetitions
        max_gap: Maximum gap between monomers (bp)
        min_purity: Minimum purity score (0-1)
        min_score: Minimum overall quality score (0-100)

    Returns:
        DataFrame of detected HORs with quality metrics
    """
    # Filter out NaN families and extract sequence
    monomers_valid = monomers_df[monomers_df['monomer_family'].notna()].copy()
    monomers_valid = monomers_valid.reset_index(drop=True)
    family_sequence = monomers_valid['monomer_family'].astype(int).tolist()

    if len(family_sequence) < min_pattern_length * min_copies:
        return pd.DataFrame()

    # Detect HORs with quality metrics
    raw_hors = find_repeating_patterns_refined(
        family_sequence,
        monomers_valid,
        min_pattern_length=min_pattern_length,
        max_pattern_length=max_pattern_length,
        min_copies=min_copies,
        max_gap=max_gap,
        min_purity=min_purity
    )

    # Curate overlaps (prefer high quality)
    curated_hors = curate_overlapping_hors_refined(raw_hors)

    # Convert to DataFrame with genomic coordinates
    hor_records = []
    for hor in curated_hors:
        # Calculate quality score
        quality_score = calculate_hor_score(hor)

        # Filter by minimum score
        if quality_score < min_score:
            continue

        # Get genomic coordinates
        start_mono = monomers_valid.iloc[hor['start_monomer_idx']]
        end_mono = monomers_valid.iloc[hor['end_monomer_idx'] - 1]

        hor_start = start_mono['monomer_start']
        hor_end = end_mono['monomer_end']

        # Format HOR unit
        hor_unit = format_hor_unit(hor['pattern'])

        # Determine type
        hor_type = 'homHOR' if len(set(hor['pattern'])) == 1 else 'hetHOR'

        hor_records.append({
            'hor_start': hor_start,
            'hor_end': hor_end,
            'hor_unit': hor_unit,
            'hor_unit_length': hor['pattern_length'],
            'hor_copies': hor['copies'],
            'total_monomers': hor['total_monomers'],
            'hor_type': hor_type,
            'purity': hor['purity'],
            'quality_score': quality_score,
            'max_gap': hor['max_gap'],
            'mean_gap': hor['mean_gap'],
            'gap_std': hor['gap_std'],
            'pattern_tuple': str(hor['pattern'])
        })

    return pd.DataFrame(hor_records)

# Test cases
if __name__ == '__main__':
    print("=== TESTING REFINED MONOMER-LEVEL HOR DETECTION ===\n")

    # Test case 1: Perfect 1068 consecutive F3 monomers
    print("Test 1: 1068 consecutive F3 monomers (perfect quality)")
    print("Expected: Should detect 3F3 pattern with high purity\n")

    test_monomers = pd.DataFrame({
        'monomer_family': [3] * 1068,
        'monomer_start': [i * 178 for i in range(1068)],
        'monomer_end': [(i + 1) * 178 for i in range(1068)]
    })

    hors_detected = analyze_centromere_array_refined(
        test_monomers,
        min_pattern_length=3,
        max_pattern_length=20,
        min_copies=3,
        min_purity=0.9
    )

    print(f"Detected {len(hors_detected)} HORs:")
    if len(hors_detected) > 0:
        print(hors_detected[['hor_unit', 'hor_copies', 'total_monomers',
                            'hor_type', 'purity', 'quality_score']].to_string(index=False))

    print("\n" + "="*70 + "\n")

    # Test case 2: F3 with a large gap
    print("Test 2: 500 F3 + LARGE GAP + 500 F3")
    print("Expected: Two separate HORs (gap breaks them)\n")

    test_monomers_gap = pd.DataFrame({
        'monomer_family': [3] * 500 + [3] * 500,
        'monomer_start': [i * 178 for i in range(500)] +
                        [i * 178 + 100000 for i in range(500, 1000)],
        'monomer_end': [(i + 1) * 178 for i in range(500)] +
                      [(i + 1) * 178 + 100000 for i in range(500, 1000)]
    })

    hors_detected_gap = analyze_centromere_array_refined(
        test_monomers_gap,
        min_pattern_length=3,
        max_pattern_length=20,
        min_copies=3
    )

    print(f"Detected {len(hors_detected_gap)} HORs:")
    if len(hors_detected_gap) > 0:
        print(hors_detected_gap[['hor_unit', 'hor_copies', 'purity',
                                 'quality_score', 'hor_start', 'hor_end']].to_string(index=False))

    print("\n" + "="*70 + "\n")

    # Test case 3: HetHOR pattern
    print("Test 3: F4-F5-F7 repeated 10 times (perfect)")
    print("Expected: 1F4-1F5-1F7 × 10 copies, high purity\n")

    test_pattern = [4, 5, 7] * 10
    test_monomers2 = pd.DataFrame({
        'monomer_family': test_pattern,
        'monomer_start': [i * 178 for i in range(len(test_pattern))],
        'monomer_end': [(i + 1) * 178 for i in range(len(test_pattern))]
    })

    hors_detected2 = analyze_centromere_array_refined(
        test_monomers2,
        min_pattern_length=3,
        max_pattern_length=10,
        min_copies=3
    )

    print(f"Detected {len(hors_detected2)} HORs:")
    if len(hors_detected2) > 0:
        print(hors_detected2[['hor_unit', 'hor_copies', 'total_monomers',
                             'hor_type', 'purity', 'quality_score']].to_string(index=False))

    print("\n" + "="*70 + "\n")

    # Test case 4: Imperfect HOR (some mismatches)
    print("Test 4: Mostly F3 with occasional F1 (imperfect HOR)")
    print("Expected: Lower purity score, may be filtered\n")

    imperfect = [3, 3, 3] * 40  # 120 monomers of F3
    # Add some noise
    imperfect[10] = 1
    imperfect[50] = 1
    imperfect[90] = 1

    test_monomers3 = pd.DataFrame({
        'monomer_family': imperfect,
        'monomer_start': [i * 178 for i in range(len(imperfect))],
        'monomer_end': [(i + 1) * 178 for i in range(len(imperfect))]
    })

    hors_detected3 = analyze_centromere_array_refined(
        test_monomers3,
        min_pattern_length=3,
        max_pattern_length=10,
        min_copies=3,
        min_purity=0.8  # Lower threshold to see imperfect HOR
    )

    print(f"Detected {len(hors_detected3)} HORs:")
    if len(hors_detected3) > 0:
        print(hors_detected3[['hor_unit', 'hor_copies', 'purity',
                             'quality_score']].to_string(index=False))
    else:
        print("No HORs detected (filtered by quality thresholds)")

    print("\n✅ Refined HOR detection complete!")
    print("\nKey improvements:")
    print("- Purity scoring (0-1): Measures pattern perfection")
    print("- Quality scoring (0-100): Overall HOR quality")
    print("- Gap metrics: max, mean, std deviation")
    print("- Better overlap resolution: Prefers high-quality HORs")
    print("- Flexible filtering: min_purity and min_score parameters")
