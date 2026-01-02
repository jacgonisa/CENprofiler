# Ribbon Plot Visualizations: Satellite Remodelling

**Date**: 2026-01-02
**Sample**: Col_9day Chr3 centromere test subset
**Purpose**: Visualize structural variation and satellite remodelling in single molecules

---

## Overview

These ribbon plots visualize how centromeric satellites are remodelled in individual DNA molecules compared to the reference genome. Each plot shows:

1. **Reference track (top)**: Satellite monomers from the reference genome
2. **Read track (bottom)**: Satellite monomers from the sequenced DNA molecule
3. **Gray ribbons**: Aligned regions connecting reference to read
4. **Gaps in ribbons**: Large indels (insertions/deletions ≥100bp)
5. **Color coding**: Each CEN178 family has a distinct color

---

## Visualization Details

### Color Scheme

- **Family 1**: Red
- **Family 2**: Blue
- **Family 3**: Green
- **Family 4**: Purple
- **Family 5**: Dark Gray
- **Family 6**: Light Gray
- **Family 7**: Orange
- **Family 8**: Yellow
- **Family 10**: Light Gray
- **Family 11**: Brown
- **Family 13**: Pink
- **Family 16**: Light Gray
- **Family 18**: Cyan
- **Unclassified**: Light gray

### Visual Elements

- **Insertions**: Highlighted with blue markers on read track
- **Deletions**: Highlighted with red markers on reference track
- **Individual monomers**: ~178bp units shown as vertical bars
- **Black borders**: Separate individual monomers within arrays

---

## Generated Plots

### 1. indel_families_1_a5a8102b.png (428 KB)

**Read**: a5a8102b-b1a4-47c6-8fb9-fe8cf7c48b74
**Region**: Chr3:13,673,187-13,744,633 bp
**Read length**: 66,348 bp
**Large indels**: 8 insertions, 1 deletion

**Key observations**:
- Multiple large insertions visible as gaps in ribbons
- Different satellite family composition in read vs reference
- Clear family transitions visible within arrays
- Shows extensive satellite array expansion

### 2. indel_families_2_d9d95b8c.png (188 KB)

**Read**: d9d95b8c-f334-42af-9e4d-cce45e7b83db
**Region**: Chr3 centromere
**Large indels**: 6 total

**Key observations**:
- Shorter alignment region with focused structural variation
- Clear satellite family patterns
- Mix of insertions and deletions

### 3. indel_families_3_172620b0.png (395 KB)

**Read**: 172620b0-df36-4b65-af36-610a46d8b2a2
**Region**: Chr3 centromere
**Large indels**: 6 total

**Key observations**:
- Long continuous satellite arrays
- Complex family composition
- Multiple structural variants

---

## Biological Interpretation

### Satellite Remodelling Patterns

1. **Copy number variation**: Ribbon gaps show satellite array expansion/contraction
2. **Family composition changes**: Color patterns differ between reference and read
3. **Structural plasticity**: Centromeric satellites are highly dynamic
4. **Tandem array dynamics**: Individual monomers can be gained or lost

### What the Ribbons Show

**Continuous ribbons** = Regions where read matches reference exactly
**Ribbon gaps on reference side** = Read has INSERTION (extra sequence)
**Ribbon gaps on read side** = Read has DELETION (missing sequence relative to reference)

### Satellite Family Dynamics

By comparing the color patterns (family composition) between reference and read:
- **Family enrichment**: Some families appear more frequently in reads
- **Family depletion**: Some families are lost or reduced
- **New arrangements**: Different family orders in read vs reference
- **Homogenization**: Conversion of one family to another

---

## Technical Details

### Input Data

- **BAM file**: test_data/bam/Col_9day_test_Chr3cen.bam
- **Indel catalog**: 93 large indels (≥100bp)
- **Monomers classified**: 13,236 CEN178 monomers
- **Family assignments**: From phylogenetic clustering

### Visualization Method

The plots use CIGAR string parsing to:
1. Extract alignment blocks (matches, insertions, deletions)
2. Map read coordinates to reference coordinates
3. Project monomer positions onto both tracks
4. Draw ribbons connecting aligned regions
5. Highlight indels with colored markers

### Plot Resolution

- **DPI**: 150
- **Figure size**: 20 × 5 inches
- **Format**: PNG
- **Monomer borders**: Visible to show individual ~178bp units

---

## Usage

These visualizations demonstrate:
- ✅ Satellite structural variation detection
- ✅ Family-level satellite composition analysis
- ✅ Single-molecule resolution of centromere dynamics
- ✅ Validation of indel calling from CIGAR strings

---

## Next Steps

To generate more ribbon plots:

```bash
python3 bin/visualize_indel_families_v2.py \
  --bam <bam_file> \
  --sv-info <indel_catalog.tsv> \
  --monomers <monomer_classifications.tsv> \
  --read-ids <read_id_1> <read_id_2> ... \
  --output-dir <output_directory>
```

To find reads with most indels:

```bash
awk -F'\t' 'NR>1 {count[$1]++} END {for (id in count) print count[id], id}' \
  indel_catalog.tsv | sort -rn | head -20
```

---

## Comparison with Reference Genome Analysis

**Genome-wide plots** (from bin/plot_monomer_level_genome_wide.py):
- Show entire chromosomes
- Display HOR patterns
- Identify large duplications

**Single-molecule ribbon plots** (this visualization):
- Show individual reads
- Compare read vs reference directly
- Visualize structural variation at monomer resolution
- Demonstrate satellite remodelling events

Both approaches are complementary and provide different insights into centromeric satellite organization.

---

**Generated**: 2026-01-02
**Pipeline**: CENprofiler v2.0
**Visualization**: visualize_indel_families_v2.py
