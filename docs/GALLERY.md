# CENprofiler Output Gallery

This document showcases example outputs from the CENprofiler pipeline.

---

## Test Dataset

**Genome**: Arabidopsis thaliana Col-CC centromeric regions
**Input**: `/mnt/ssd-8tb/atrx_china/CEN178profiler/v3_genome_analysis_CORRECT/centromeric_regions.fa`
**Reference**: 5,597 CEN178 representative monomers (20 phylogenetic families)
**Results**: `test_data/results_genome/`

---

## Pipeline Statistics

### Monomer-Level Analysis
- **Total monomers detected**: 63,233
- **Tandem arrays analyzed**: Multiple arrays across Chr1-Chr5
- **Period range**: 160-200 bp (CEN178 satellite)

### HOR Detection
(Gap-aware detection with min_copies ≥ 3 AND monomers_per_unit ≥ 3)

- **Total HORs detected**: 355
- **homHORs**: 184 (51.8%) - Single family repeated
- **hetHORs**: 171 (48.2%) - Multiple families in pattern

---

##  HOR Architecture Schematics

### 1. homHOR Patterns (Top 10)

![homHOR Schematics](../test_data/results_genome/05_plots/hors/hor_schematics_homHOR_monomer_level.png)

**Top homHOR patterns detected:**

| Rank | Pattern | Occurrences | Description |
|------|---------|-------------|-------------|
| 1 | **3F3** (1F3-1F3-1F3) | 90 | Family 3 repeated 3 times |
| 2 | **3F1** (1F1-1F1-1F1) | 31 | Family 1 repeated 3 times |
| 3 | **3F4** (1F4-1F4-1F4) | 27 | Family 4 repeated 3 times |
| 4 | **3F5** (1F5-1F5-1F5) | 14 | Family 5 repeated 3 times |
| 5 | **3F6** (1F6-1F6-1F6) | 7 | Family 6 repeated 3 times |
| 6 | **3F8** (1F8-1F8-1F8) | 7 | Family 8 repeated 3 times |
| 7 | **3F7** (1F7-1F7-1F7) | 3 | Family 7 repeated 3 times |

**Key Observations:**
- Family 3 (F3) forms the most abundant homHOR pattern (90 occurrences)
- homHORs comprise simple tandem repetitions of a single monomer family
- Most homHORs have exactly 3 copies (minimum threshold)

---

### 2. hetHOR Patterns (Top 10)

![hetHOR Schematics](../test_data/results_genome/05_plots/hors/hor_schematics_hetHOR_monomer_level.png)

**Top hetHOR patterns detected:**

| Rank | Pattern | Occurrences | Description |
|------|---------|-------------|-------------|
| 1 | **1F6-1F2-1F2** | 14 | F6 followed by 2×F2 |
| 2 | **1F4-1F5-1F4** | 6 | F4-F5-F4 unit |
| 3 | **1F1-1F7-1F1** | 5 | F1-F7-F1 unit |
| 4 | **1F5-1F5-1F4** | 4 | 2×F5 followed by F4 |
| 5 | **1F4-1F5-1F5** | 4 | F4 followed by 2×F5 |
| 6 | **1F2-1F6-1F2** | 4 | F2-F6-F2 unit |

**Key Observations:**
- hetHORs show diverse family combinations
- F6-F2-F2 is the most common hetHOR pattern (14 occurrences)
- Many hetHORs display palindromic-like structures (e.g., F1-F7-F1, F2-F6-F2)
- Some patterns involve 4+ monomers per unit (e.g., 1F4-1F5-1F4-1F5)

---

## Output Directory Structure

```
test_data/results_genome/
├── 01_tandem_detection/
│   ├── arrays.bed              # Tandem array coordinates
│   └── arrays.1aln             # FasTAN alignment output
├── 02_monomers/
│   ├── monomers.fa             # Extracted individual monomers
│   └── monomer_info.tsv        # Monomer metadata
├── 03_classification/
│   └── monomer_classifications.tsv   # Family assignments
├── 03_hors/
│   ├── hors_detected.tsv       # All HORs with coordinates
│   └── large_duplications.tsv  # HORs ≥10kb
├── 04_stats/
│   └── chromosome_stats.tsv    # Per-chromosome statistics
└── 05_plots/
    ├── hors/
    │   ├── hor_schematics_homHOR_monomer_level.png
    │   ├── hor_schematics_hetHOR_monomer_level.png
    │   └── hor_plots.log
    └── satellites/
        └── satellite_plots.log
```

---

## Data Files

### HORs Detected (`hors_detected.tsv`)

Contains all detected HORs with the following columns:
- `seq_id`: Sequence/chromosome identifier
- `array_idx`: Array index within sequence
- `hor_start`, `hor_end`: Genomic coordinates
- `hor_unit`: Pattern notation (e.g., "1F3-1F3-1F3")
- `hor_unit_length`: Number of monomers per repeat unit
- `hor_copies`: Number of times the unit is repeated
- `total_monomers`: Total monomers in the HOR
- `hor_type`: homHOR or hetHOR
- `hor_length_bp`: Length in basepairs
- `hor_length_kb`: Length in kilobases

### Chromosome Statistics (`chromosome_stats.tsv`)

Per-chromosome summary statistics including:
- Total monomers and classified monomers
- Family distribution
- HOR counts (total, homHOR, hetHOR)
- Family-specific HOR patterns

---

## Visualization Types

### Currently Implemented

1. **HOR Architecture Schematics** ✅
   - Visual representation of top homHOR and hetHOR patterns
   - Shows monomer composition with family colors
   - Displays occurrence frequency

### To Be Integrated

2. **Genome-Wide HOR Distribution**
   - Chromosome tracks showing HOR locations
   - Coverage density plots
   - Requires column name fix: `seq_id` → `read_id` compatibility

3. **Satellite Monomer Plots**
   - Family enrichment analysis
   - Monomer distribution across chromosomes
   - Requires script path fixes for integration

4. **Large Duplication Visualizations**
   - Overview of duplications ≥10kb
   - Detailed zooms of specific regions

---

## Future Enhancements

### Genome Mode
- [ ] Fix genome-wide HOR distribution plot (column name compatibility)
- [ ] Integrate satellite enrichment analysis
- [ ] Add large duplication visualizations
- [ ] Chromosome-level comparison plots

### Read Mode
- [ ] Requires architecture redesign (see below)
- [ ] Per-read monomer ribbon diagrams
- [ ] Indel/SV family composition plots
- [ ] Family transition heatmaps

---

## Notes

### Genome Mode vs Read Mode

The current pipeline distinguishes between two modes, but **architecture improvements are needed**:

**Current Approach:**
- `--mode genome`: Analyze reference genome
- `--mode reads`: Analyze long reads

**Proposed Approach:**
The read mode needs additional inputs for proper indel/SV analysis:
- If input is **genome FASTA** → analyze genome directly
- If input is **reads FASTA** → also require:
  - Reference genome (what reads were aligned to)
  - Alignment file (BAM/SAM)

This allows proper indel detection and family-specific variant analysis.

---

## Example Command

```bash
# Genome mode (current test)
nextflow run main.nf \
    --mode genome \
    --input centromeric_regions.fa \
    --reference_monomers Col-CC-V2-CEN178-representative.fasta \
    --family_assignments itol_manual_phylo_clusters.txt \
    --outdir results/
```

---

**Generated**: 2026-01-02
**Pipeline Version**: v1.0 (initial release)
**Test Dataset**: Arabidopsis thaliana Col-CC centromeres
