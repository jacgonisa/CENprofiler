# CENprofiler Workflow: Generating Ribbon Plots

## Overview

This document explains in detail how CENprofiler was used to generate the **ribbon plots** (e.g., `indel_families_2_d9d95b8c.png`) that visualize indel patterns across centromeric tandem repeat arrays. The workflow consists of multiple steps from raw sequencing data to final visualization.

## Table of Contents

1. [Input Data Requirements](#input-data-requirements)
2. [Complete Workflow](#complete-workflow)
3. [Testing on Single Molecules](#testing-on-single-molecules)
4. [Detailed Step-by-Step Guide](#detailed-step-by-step-guide)
5. [Output Files Explanation](#output-files-explanation)
6. [Troubleshooting](#troubleshooting)

---

## Input Data Requirements

### Required Files

1. **Query sequence(s)**: FASTA file containing centromeric region(s)
   - Can be single molecule (e.g., long-read HiFi sequencing)
   - Can be multiple molecules or genome assembly
   - Example: `centromeric_region.fa` or `single_molecule.fa`

2. **Reference monomers**: FASTA file with representative centromeric monomers
   - Example: `Col-CC-V2-CEN178-representative.fasta`
   - Contains ~178bp centromeric repeat units
   - One sequence per monomer family

3. **Family assignments**: TSV file mapping monomer IDs to family classifications
   - Example: `itol_manual_phylo_clusters.txt`
   - Format: `monomer_id<TAB>family_number`

### Optional Files

- **Indel family classification**: TSV file with indel→family assignments (auto-generated if not provided)

---

## Complete Workflow

### High-Level Pipeline Steps

```
Input FASTA (centromeric regions)
    ↓
[1] MONOMER DETECTION (FASTaN)
    ↓
[2] MONOMER CLASSIFICATION (minimap2 + family assignment)
    ↓
[3] ARRAY SEGMENTATION (identify tandem repeat arrays)
    ↓
[4] INDEL DETECTION (alignment-based)
    ↓
[5] INDEL FAMILY CLASSIFICATION (phylogenetic grouping)
    ↓
[6] RIBBON PLOT GENERATION (visualization)
```

### Execution Command

**For Nextflow pipeline (recommended):**

```bash
nextflow run main.nf \
    --mode genome \
    --input centromeric_regions.fa \
    --reference_monomers assets/Col-CC-V2-CEN178-representative.fasta \
    --family_assignments assets/itol_manual_phylo_clusters.txt \
    --outdir results/ \
    --fastan_threads 8 \
    --minimap2_threads 4
```

**For single molecule analysis:**

```bash
nextflow run main.nf \
    --mode genome \
    --input single_molecule.fa \
    --reference_monomers assets/Col-CC-V2-CEN178-representative.fasta \
    --family_assignments assets/itol_manual_phylo_clusters.txt \
    --outdir results_single_molecule/ \
    --fastan_threads 4 \
    --minimap2_threads 2
```

---

## Testing on Single Molecules

### Quick Start for Single Molecule Analysis

If you have a **single long-read sequencing read** spanning a centromeric region:

#### Step 1: Prepare Your Data

```bash
# Create a FASTA file with your single molecule
cat > single_molecule.fa << 'EOF'
>molecule_1
ACGTACGT...  # Your sequence here
EOF
```

#### Step 2: Run CENprofiler

```bash
# Run the pipeline
nextflow run main.nf \
    --mode genome \
    --input single_molecule.fa \
    --reference_monomers assets/Col-CC-V2-CEN178-representative.fasta \
    --family_assignments assets/itol_manual_phylo_clusters.txt \
    --outdir results_molecule1/
```

#### Step 3: Generate Ribbon Plots

```bash
# After pipeline completes, run visualization
python bin/visualize_indel_families_v2.py \
    results_molecule1/02_monomers/monomer_classifications.tsv \
    results_molecule1/03_indels/indel_classifications.tsv \
    results_molecule1/04_visualizations/
```

### Expected Runtime

- **Single molecule (1-5 Mb)**: 5-15 minutes
- **Multiple molecules (10-20 Mb)**: 30-60 minutes
- **Whole genome centromeres (50-100 Mb)**: 2-4 hours

---

## Detailed Step-by-Step Guide

### Step 1: Monomer Detection (FASTaN)

**What it does:** Identifies individual ~178bp centromeric repeat units (monomers) in the input sequence.

**Tool:** FASTaN (tandem repeat finder optimized for centromeres)

**Parameters:**
```bash
fastan \
    --input centromeric_region.fa \
    --output monomers.fa \
    --threads 8 \
    --min-length 160 \
    --max-length 200
```

**Output:**
- `monomers.fa`: FASTA file with detected monomers
- Each sequence header contains position info: `>seq_id:start-end`

### Step 2: Monomer Classification

**What it does:** Assigns each monomer to a phylogenetic family (F1-F20) based on sequence similarity.

**Tool:** minimap2 + custom classification script

**Process:**
1. Align each monomer to reference monomers
2. Find best match based on alignment identity
3. Assign family based on reference monomer's classification

**Command:**
```bash
# Align to reference
minimap2 -x map-ont \
    assets/Col-CC-V2-CEN178-representative.fasta \
    monomers.fa > alignments.paf

# Classify based on alignments
python bin/classify_monomers.py \
    --alignments alignments.paf \
    --family-assignments assets/itol_manual_phylo_clusters.txt \
    --output monomer_classifications.tsv
```

**Output:** `monomer_classifications.tsv`
```
seq_id              monomer_idx  monomer_start  monomer_end  monomer_family  array_idx
Chr4_centromere     0            1000           1178         3               0
Chr4_centromere     1            1178           1356         3               0
Chr4_centromere     2            1356           1534         4               0
```

### Step 3: Array Segmentation

**What it does:** Groups consecutive monomers into tandem repeat arrays.

**Logic:**
- Consecutive monomers with gaps <1000bp → same array
- Large gaps (≥1000bp) → new array starts
- Assigns unique `array_idx` to each array

**Output:** Updated `monomer_classifications.tsv` with `array_idx` column

### Step 4: Indel Detection

**What it does:** Identifies insertions and deletions within monomers by comparing to reference.

**Algorithm:**
1. For each monomer, align to its family reference
2. Parse CIGAR string to find indels
3. Classify indel location within monomer
4. Record indel sequence and surrounding context

**Command:**
```bash
python bin/detect_indels.py \
    --monomers monomers.fa \
    --classifications monomer_classifications.tsv \
    --reference assets/Col-CC-V2-CEN178-representative.fasta \
    --output indels_raw.tsv
```

**Output:** `indels_raw.tsv`
```
monomer_id           indel_type  indel_start  indel_end  indel_length  indel_seq      context_before  context_after
Chr4_centromere:1000 deletion    50           55         5             -----          ACGTAC          GTACGT
Chr4_centromere:1178 insertion   100          100        3             ATG            GCTAGC          CTAGCT
```

### Step 5: Indel Family Classification

**What it does:** Groups similar indels into families based on:
- Indel type (insertion/deletion)
- Length
- Sequence similarity
- Position within monomer

**Algorithm:**
```python
# Pseudocode
for each indel:
    - Extract indel features (type, length, sequence, position)
    - Compare to existing indel families
    - If similarity > 80%: assign to existing family
    - Else: create new family
```

**Command:**
```bash
python bin/classify_indel_families.py \
    --indels indels_raw.tsv \
    --monomers monomer_classifications.tsv \
    --output indel_classifications.tsv \
    --min-similarity 0.8
```

**Output:** `indel_classifications.tsv`
```
monomer_id           monomer_family  indel_family  indel_type  indel_start  indel_end  indel_length
Chr4_centromere:1000 3               DEL_F3_1      deletion    50           55         5
Chr4_centromere:1178 3               INS_F3_2      insertion   100          100        3
```

### Step 6: Ribbon Plot Generation

**What it does:** Creates the final visualization showing:
- Monomer positions and family colors
- Indel positions marked as vertical lines
- Different colors for indel families
- Statistical summaries

**Command:**
```bash
python bin/visualize_indel_families_v2.py \
    monomer_classifications.tsv \
    indel_classifications.tsv \
    output_dir/ \
    --array-idx 2  # Optional: visualize specific array
```

**Visualization Components:**

1. **Top panel**: Monomer ribbon
   - Each rectangle = one monomer
   - Color = monomer family (F1=red, F2=blue, F3=green, etc.)
   - Width = monomer length

2. **Indel markers**: Vertical lines on monomers
   - Color = indel family
   - Position = location within monomer
   - Height indicates multiple indels at same position

3. **Bottom panel**: Position scale (kb)

4. **Statistics box**:
   - Total indels detected
   - Indel families identified
   - Most common indel types

**Output:** `indel_families_2_d9d95b8c.png` (example filename)

---

## Output Files Explanation

### Directory Structure

```
results/
├── 01_tandem_repeats/
│   ├── monomers.fa                    # Detected monomers (FASTaN output)
│   └── monomer_positions.tsv          # Monomer coordinates
│
├── 02_monomers/
│   ├── alignments.paf                 # Minimap2 alignments
│   └── monomer_classifications.tsv    # Classified monomers with families
│
├── 03_indels/
│   ├── indels_raw.tsv                 # All detected indels
│   ├── indel_classifications.tsv      # Indels grouped into families
│   └── indel_statistics.txt           # Summary statistics
│
└── 04_visualizations/
    ├── indel_families_0_*.png         # Ribbon plot for array 0
    ├── indel_families_1_*.png         # Ribbon plot for array 1
    ├── indel_families_2_*.png         # Ribbon plot for array 2
    └── ...
```

### Key Output Files

**`monomer_classifications.tsv`**
- Contains all monomers with their positions and family assignments
- Used as input for downstream analyses
- Columns: seq_id, monomer_idx, monomer_start, monomer_end, monomer_family, array_idx

**`indel_classifications.tsv`**
- Contains all indels grouped into families
- Links indels to their parent monomers
- Columns: monomer_id, monomer_family, indel_family, indel_type, indel_start, indel_end, indel_length

**`indel_families_*.png`**
- Visual representation of one tandem repeat array
- Shows monomer organization and indel patterns
- File naming: `indel_families_{array_idx}_{random_id}.png`

---

## Troubleshooting

### Issue 1: No monomers detected

**Symptoms:** FASTaN output is empty or very few monomers

**Solutions:**
```bash
# Check input sequence length
seqkit stats centromeric_region.fa

# Adjust FASTaN parameters
fastan --min-length 150 --max-length 210  # More lenient

# Verify you're using centromeric regions (not whole genome)
```

### Issue 2: Poor monomer classification

**Symptoms:** Many monomers classified as "unknown" or wrong families

**Solutions:**
```bash
# Check reference monomer quality
seqkit stats assets/Col-CC-V2-CEN178-representative.fasta

# Verify family assignments file format
head assets/itol_manual_phylo_clusters.txt

# Lower stringency in minimap2
minimap2 -x asm20  # More sensitive mode
```

### Issue 3: No indels detected

**Symptoms:** `indels_raw.tsv` is empty

**Solutions:**
```bash
# Check alignment quality
samtools view -h alignments.bam | head

# Adjust indel detection sensitivity
python bin/detect_indels.py --min-indel-length 1  # Detect smaller indels

# Verify reference sequences match query species
```

### Issue 4: Ribbon plots not generated

**Symptoms:** Visualization script fails or produces blank images

**Solutions:**
```bash
# Check input files exist and have data
wc -l monomer_classifications.tsv
wc -l indel_classifications.tsv

# Run with debug mode
python bin/visualize_indel_families_v2.py ... --debug

# Check matplotlib backend
python -c "import matplotlib; print(matplotlib.get_backend())"
```

---

## Example: Complete Single-Molecule Workflow

Here's a complete example from start to finish:

```bash
#!/bin/bash
# complete_workflow.sh - Process a single centromeric molecule

# Configuration
MOLECULE="single_molecule.fa"
REFERENCE="assets/Col-CC-V2-CEN178-representative.fasta"
FAMILIES="assets/itol_manual_phylo_clusters.txt"
OUTDIR="results_molecule1"

# Step 1: Run CENprofiler pipeline
echo "Running CENprofiler pipeline..."
nextflow run main.nf \
    --mode genome \
    --input $MOLECULE \
    --reference_monomers $REFERENCE \
    --family_assignments $FAMILIES \
    --outdir $OUTDIR \
    --fastan_threads 4 \
    --minimap2_threads 2

# Step 2: Generate ribbon plots
echo "Generating ribbon plots..."
python bin/visualize_indel_families_v2.py \
    $OUTDIR/02_monomers/monomer_classifications.tsv \
    $OUTDIR/03_indels/indel_classifications.tsv \
    $OUTDIR/04_visualizations/

# Step 3: Generate summary statistics
echo "Generating statistics..."
python bin/analyze_indel_patterns.py \
    $OUTDIR/03_indels/indel_classifications.tsv \
    $OUTDIR/05_statistics/

echo "Analysis complete! Results in: $OUTDIR"
echo "Ribbon plots: $OUTDIR/04_visualizations/"
```

---

## Citation

If you use CENprofiler in your research, please cite:

```
CENprofiler: A comprehensive pipeline for centromeric tandem repeat analysis
[Your citation information here]
```

---

## Contact & Support

For questions, issues, or feature requests:
- GitHub Issues: https://github.com/jacgonisa/CENprofiler/issues
- Email: [Your contact email]

Last updated: 2026-01-14
