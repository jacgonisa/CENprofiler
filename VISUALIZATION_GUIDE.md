# CENprofiler Visualization Guide

**Date**: 2026-01-02
**Pipeline Version**: v2.0

---

## Overview

CENprofiler generates multiple types of visualizations to understand satellite organization and remodelling at different scales:

1. **Ribbon plots**: Single molecule vs reference - shows how satellites are remodelled
2. **Summary plots**: Family distributions and statistics
3. **Array plots**: Individual tandem array structures
4. **Per-read plots**: Classification rates and indel distributions

---

## 1. Ribbon Plots - Single Molecule vs Reference

### Location

**New test results** (CENprofiler pipeline):
```
test_data/results_test/06_ribbon_plots/
├── indel_families_1_a5a8102b.png  (428 KB) - 8 INS + 1 DEL
├── indel_families_2_d9d95b8c.png  (188 KB) - 6 indels
├── indel_families_3_172620b0.png  (395 KB) - 6 indels
└── VISUALIZATION_SUMMARY.md
```

**Previous analysis** (CEN178profiler - more comprehensive):
```
../CEN178profiler/results_v2_test/
├── ribbon_monomer_1_5b96a790.png  (315 KB) ⭐ EXCELLENT
├── ribbon_monomer_2_6da52b44.png  (491 KB) ⭐ EXCELLENT
├── indel_families_1_5b96a790.png  (183 KB)
├── indel_families_2_6da52b44.png  (408 KB)
├── indel_families_3_215a2258.png  (370 KB)
└── README.md  (comprehensive documentation)
```

### What Ribbon Plots Show

**Two tracks:**
- **Reference track (top)**: Satellite monomers from reference genome
- **Read track (bottom)**: Satellite monomers from sequenced DNA molecule

**Visual elements:**
- **Gray ribbons**: Connect aligned regions (read matches reference)
- **Ribbon gaps**: Show large indels (insertions/deletions ≥100bp)
- **Color bars**: Individual ~178bp monomers, colored by CEN178 family
- **Black borders**: Separate individual monomers within arrays

**Insertions**: Extra sequence in read (gap on reference side)
**Deletions**: Missing sequence in read (gap on read side, highlighted in red/pink)

### Biological Interpretation

✅ **Satellite copy number variation**: Ribbon gaps show array expansion/contraction
✅ **Family composition changes**: Different color patterns between reference and read
✅ **Structural plasticity**: Centromeric satellites are highly dynamic
✅ **Tandem array remodelling**: Individual monomers gained or lost

---

## 2. Summary Plots

### Family Distribution

**Location**: `test_data/results_test/05_plots/reads/family_distribution.png`

**Shows**: Bar chart of CEN178 family counts across all reads

**Previous analysis**: `../CEN178profiler/results_v2_test/family_summary.png`
- Two panels: family distribution + transition heatmap
- Reveals which families transition into each other

### Indel Distribution

**Location**: `test_data/results_test/05_plots/reads/indel_distribution.png`

**Shows**: Histogram of indel sizes (insertions vs deletions)

### Read Statistics

**Location**: `test_data/results_test/05_plots/reads/read_statistics.png`

**Shows**:
- Array size distribution (monomers per array)
- Classification rate per read (% monomers classified)

---

## 3. Array Structure Plots

**Location**: `../CEN178profiler/results_v2_test/`

### Top Arrays Combined

**File**: `top_arrays_combined.png` (244 KB)

**Shows**: Top 3 largest satellite arrays side-by-side
- Per-monomer family coloring
- Reveals family heterogeneity within single arrays

### Individual Arrays

- `array_1_idx70.png` (112 KB) - 99 monomers, 17.6 kb
- `array_2_idx69.png` (109 KB) - 82 monomers, 14.5 kb
- `array_3_idx5.png` (76 KB) - 8 monomers
- `array_4_idx42.png` (81 KB) - 6 monomers
- `array_5_idx0.png` (75 KB) - 5 monomers

**Key finding**: Arrays contain 4-8 different families with frequent transitions!

---

## 4. Color Scheme - CEN178 Families

All visualizations use consistent family colors:

| Family | Color | Description |
|--------|-------|-------------|
| 1 | Red | Dominant family (58.9% in test) |
| 2 | Blue | Minor |
| 3 | Green | Minor |
| 4 | Purple | Minor |
| 7 | Orange/Brown | Common (16.1% in test) |
| 8 | Yellow | Rare |
| 11 | Brown/Orange | Common (12.2% in test) |
| 13 | Pink | Rare |
| 18 | Cyan | Common (7.2% in test) |
| Unclassified | Gray | Not assigned |

---

## 5. Key Results from CEN178profiler v2 Analysis

### Summary Stats

```
Total monomers:     209
Classified:         180 (86.1%)
Families detected:  9/20
Largest array:      99 monomers (17.6 kb)
```

### Major Discovery: Family Heterogeneity

**Previous expectation**: Tandem arrays should be homogeneous
**Actual observation**: Arrays contain 4-8 different families with frequent transitions

**Example** (Array 69, first 10 monomers):
```
Position    Family  Identity
0-178       Fam 7   97.6%
178-356     Fam 1   98.8%
356-534     Fam 18  100%   ← Perfect match
534-712     Fam 18  100%
712-890     Fam 7   96.8%
890-1068    Fam 1   100%
1068-1246   Fam 11  100%
1246-1424   Fam 11  100%
1424-1602   Fam 7   98.9%
1602-1780   Fam 4   88.9%
```

**Biological interpretation**: Frequent recombination, gene conversion, or recent expansion events mixing multiple evolutionary lineages.

---

## 6. How to Generate Visualizations

### Ribbon Plots

```bash
python3 bin/visualize_indel_families_v2.py \
  --bam <bam_file> \
  --sv-info <indel_catalog.tsv> \
  --monomers <monomer_classifications.tsv> \
  --read-ids <read_id_1> <read_id_2> ... \
  --output-dir <output_directory>
```

**Find reads with most indels:**
```bash
awk -F'\t' 'NR>1 && $6>=100 {count[$1]++} END {for (id in count) print count[id], id}' \
  indel_catalog.tsv | sort -rn | head -20
```

### Full Pipeline

```bash
nextflow run main.nf \
  --alignment <bam_file> \
  --reference_genome <genome.fna> \
  --annotation_dir <annotation_dir> \
  --reference_monomers <CEN178_representatives.fasta> \
  --family_assignments <phylo_clusters.txt> \
  --sample_name <sample> \
  --min_indel_size 100 \
  --outdir results/
```

---

## 7. File Formats

### Indel Catalog

**File**: `01_extracted_reads/*_indel_catalog.tsv`

**Columns**:
- read_id, chromosome, ref_pos, read_pos
- type (Insertion/Deletion), size
- region (centromere/pericentromere/etc)
- mapping_quality
- multiple_of_178, closest_178_multiple, distance_to_178_multiple

### Monomer Classifications

**File**: `02_monomers/monomer_classifications.tsv`

**Columns**:
- monomer_id, seq_id, array_idx, monomer_idx
- monomer_start, monomer_end, monomer_length
- array_period, array_quality
- best_match, alignment_identity, mapq
- monomer_family (1-20 or NA)

---

## 8. Comparison: Multiple Visualization Approaches

### CEN178profiler Results (Previous - Very Detailed)

**Location**: `../CEN178profiler/results_v2_test/`

**Strengths**:
- ✅ Comprehensive README with biological interpretation
- ✅ Multiple ribbon plot styles
- ✅ Family transition heatmaps
- ✅ Per-array detailed plots
- ✅ Combined top arrays view
- ✅ Extensive documentation

**Coverage**: 5 reads, 209 monomers, 86.1% classified

### CENprofiler v2.0 Test (Current - Production Pipeline)

**Location**: `test_data/results_test/`

**Strengths**:
- ✅ Full automated pipeline
- ✅ Auto-detection of analysis mode
- ✅ BAM file support
- ✅ Integrated visualization generation
- ✅ Production-ready workflow

**Coverage**: 55 reads, 13,236 monomers, high classification rate

### Recommendation

For **publication-quality figures** → Use **CEN178profiler results**
- Better documented
- More plot types
- Detailed biological interpretation

For **large-scale analysis** → Use **CENprofiler v2.0 pipeline**
- Automated end-to-end
- Handles full BAM files
- Scalable to whole genomes

---

## 9. Next Steps

### Immediate Visualizations Needed

1. ✅ Ribbon plots showing satellite remodelling (DONE)
2. ✅ Family distributions (DONE)
3. [ ] Compare Col_9day vs atxr56_9day family distributions
4. [ ] Chromosome-specific ribbon plots
5. [ ] Statistical tests for family enrichment/depletion

### Full-Scale Analysis

1. [ ] Run full Col_9day.bam through pipeline
2. [ ] Run full atxr56_9day.bam through pipeline
3. [ ] Generate comparative visualizations
4. [ ] Statistical analysis of differences
5. [ ] Create publication-ready figure gallery

### Advanced Visualizations

1. [ ] Interactive HTML plots (plotly/bokeh)
2. [ ] Sequence logos for each family
3. [ ] Phylogenetic tree with read positions mapped
4. [ ] Network diagrams of family transitions
5. [ ] 3D visualization of family clustering

---

## 10. Key Visualizations for Understanding Satellite Remodelling

### Must-See Plots

1. **ribbon_monomer_1_5b96a790.png** ⭐⭐⭐
   - Shows deletion clearly highlighted in red/pink
   - Perfect example of satellite loss
   - Individual monomers visible

2. **top_arrays_combined.png** ⭐⭐⭐
   - Shows family heterogeneity
   - Demonstrates frequent family transitions
   - Clear monomer-level resolution

3. **family_summary.png** ⭐⭐
   - Family distribution bar chart
   - Transition heatmap shows which families neighbor each other

4. **test_data/results_test/06_ribbon_plots/indel_families_1_a5a8102b.png** ⭐⭐⭐
   - Large read with 8 insertions + 1 deletion
   - Shows satellite array expansion
   - Clear family patterns

---

## Contact

For questions about visualizations:
- Email: jg2070@cam.ac.uk
- Check `README.md` files in each results directory
- See `VISUALIZATION_SUMMARY.md` in ribbon_plots directories

---

**CENprofiler v2.0** - Production Ready
**Status**: All visualization types validated ✅
**Generated**: 2026-01-02
