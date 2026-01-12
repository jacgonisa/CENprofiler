# Monomer-Level Analysis Guide

## Overview

CENprofiler v2.0 provides comprehensive monomer-level analysis capabilities focused on understanding centromeric satellite organization, family composition, and spatial patterns. This guide covers all monomer-level analyses available in the pipeline.

## Table of Contents

1. [Automated Pipeline Analyses](#automated-pipeline-analyses)
2. [Manual Analysis Scripts](#manual-analysis-scripts)
3. [Output Directories](#output-directories)
4. [Analysis Workflows](#analysis-workflows)
5. [Interpretation Guide](#interpretation-guide)

---

## Automated Pipeline Analyses

These analyses run automatically when you execute the pipeline in READ_MODE_WITH_INDELS:

### 1. Basic Monomer Classification
**Output:** `02_monomers/monomer_classifications.tsv`

Main data file containing all classified monomers with:
- Family assignments (1-20)
- Alignment identity scores
- Array membership
- Genomic coordinates
- Quality metrics

### 2. Read-Level Plots
**Output:** `05_plots/reads/`

Basic visualizations:
- Family distribution bar charts
- Indel size distributions
- Read statistics summaries

### 3. Comprehensive Read Plots
**Output:** `05_plots/comprehensive/`

Advanced visualizations:
- **Family summary** with transition heatmap (KEY)
- Top arrays combined view
- Individual array details (top 5)
- Array summary statistics

**Family Transition Heatmap** shows which families neighbor each other sequentially, revealing:
- Recombination hotspots
- Preferred family transitions
- Array heterogeneity patterns

### 4. Ribbon Plots
**Output:** `05_plots/ribbon_plots/`

Single-molecule visualizations showing:
- Reference genome track (top)
- Individual read track (bottom)
- Gray ribbons connecting aligned regions
- Gaps showing large indels (≥100bp)
- Monomers colored by family
- Deletion monomers with red borders

### 5. Deletion Monomer Analysis
**Output:** `03_deletion_monomers/`

Analysis of satellite families LOST in deletions:
- Extracted reference sequences from deleted regions
- Classified monomers from deletions
- Integrated into ribbon plots
- Shows what was present in reference but absent in reads

### 6. Monomer Statistics ⭐ NEW
**Output:** `07_monomer_statistics/`

Comprehensive statistical analysis including:

**Plots:**
- `length_distribution.png` - Overall and per-family length distributions
- `identity_distribution.png` - Alignment identity distributions
- `family_composition.png` - Bar chart with percentages
- `transition_matrix.png` - Family transition probability heatmap
- `heterogeneity_metrics.png` - 4-panel diversity metrics
- `array_size_vs_diversity.png` - Scatter plots

**Reports:**
- `monomer_statistics.txt` - Human-readable comprehensive report
- `monomer_statistics.json` - Machine-readable JSON format
- `family_statistics.tsv` - Per-family summary table
- `array_heterogeneity.tsv` - Per-array diversity metrics

**Metrics Calculated:**
- Classification rates
- Length/identity statistics per family
- Transition frequencies and probabilities
- Shannon entropy (array complexity)
- Simpson diversity (1 - sum of squared proportions)
- Run length analysis (consecutive same-family stretches)
- Number of families per array
- Dominant family fractions

### 7. Monomer Sequence Extraction ⭐ NEW
**Output:** `08_monomer_sequences/`

Organizes sequences by family for downstream analysis:

**Outputs:**
- `family_fastas/family_*.fa` - One FASTA per family
- `family_diversity.tsv` - Within-family pairwise identity
- `family_diversity_report.txt` - Summary statistics
- `sequence_summary.txt` - Overall summary
- `consensus/all_consensus.fa` - Consensus sequences (if MUSCLE installed)

**Analysis:**
- Pairwise sequence identity within families
- Mean/std length per family
- Diversity metrics (heterogeneity of sequences)
- Consensus generation (requires MUSCLE aligner)

---

## Manual Analysis Scripts

These scripts can be run manually on classification outputs for additional insights:

### 1. Sample Comparison ⭐ NEW

```bash
python bin/compare_samples.py \\
    results_Col0/02_monomers/monomer_classifications.tsv \\
    results_mutant/02_monomers/monomer_classifications.tsv \\
    "Col-0" \\
    "atxr5/6" \\
    comparison_output/
```

**What it does:**
- Compares family composition between samples (chi-square test)
- Tests transition pattern differences
- Compares heterogeneity metrics (t-tests)
- Identifies enriched/depleted families
- Statistical significance testing throughout

**Outputs:**
- `family_comparison.png` - Side-by-side bar charts + difference plot
- `heterogeneity_comparison.png` - Boxplots for diversity metrics
- `transition_comparison.png` - Top transition differences
- `comparison_report.txt` - Detailed statistical report
- `comparison_report.json` - JSON format
- `family_comparison.tsv` - Table with fold changes
- `transition_comparison.tsv` - Transition differences

**Use cases:**
- Col-0 vs mutant comparisons
- Developmental stage comparisons
- Treatment effect analysis
- Any pairwise sample comparison

### 2. Positional/Spatial Analysis ⭐ NEW

```bash
python bin/analyze_monomer_positions.py \\
    results/02_monomers/monomer_classifications.tsv \\
    position_analysis_output/
```

**What it analyzes:**
- **Positional preferences**: Do families prefer array starts, centers, or ends?
- **Boundary enrichment**: Which families occur at array boundaries?
- **Clustering**: Are same-family monomers grouped together?
- **Co-occurrence**: Which families are found near each other?

**Outputs:**
- `positional_preferences.png` - Distribution histograms + mean positions
- `clustering_scores.png` - Clustering index per family
- `cooccurrence_heatmap.png` - Family pair enrichment matrix
- `position_analysis.txt` - Comprehensive report
- `position_analysis.json` - JSON format
- `boundary_enrichment.tsv` - Start/end/middle frequencies
- `cooccurrence.tsv` - Observed vs expected co-occurrence

**Metrics:**
- **Relative position**: 0.0 = array start, 0.5 = center, 1.0 = end
- **Clustering index**: >1 = clustered, <1 = dispersed, ~1 = random
- **Boundary enrichment**: % at start vs end vs middle
- **Co-occurrence enrichment**: observed/expected within 5bp

**Discoveries enabled:**
- Spatial organization patterns
- Non-random family arrangements
- Specific family associations
- Array boundary specialization

### 3. Run Standalone Statistics

If you want to re-run statistics without the full pipeline:

```bash
python bin/analyze_monomer_statistics.py \\
    results/02_monomers/monomer_classifications.tsv \\
    custom_stats_output/
```

### 4. Extract Sequences Separately

If you want just the sequence extraction:

```bash
python bin/extract_monomer_sequences.py \\
    results/02_monomers/monomer_classifications.tsv \\
    results/02_monomers/monomers.fa \\
    sequence_output/
```

---

## Output Directories

Full pipeline READ_MODE_WITH_INDELS produces:

```
results/
├── 00_regions/              # Genomic region annotations
├── 01_extracted_reads/      # Reads with large indels
├── 01_fastan/              # Tandem repeat detection
├── 02_monomers/            # Monomer classification (MAIN)
├── 03_deletion_monomers/   # Deleted satellite analysis
├── 05_plots/
│   ├── reads/              # Basic visualizations
│   ├── comprehensive/      # Advanced read plots
│   └── ribbon_plots/       # Single-molecule views
├── 07_monomer_statistics/  # Statistical analysis ⭐ NEW
└── 08_monomer_sequences/   # Sequence extraction ⭐ NEW
```

---

## Analysis Workflows

### Workflow 1: Single Sample Analysis

For analyzing one sample (e.g., Col-0):

1. **Run pipeline:**
   ```bash
   nextflow run main.nf \\
       --input dummy.fa \\
       --alignment Col_9day.bam \\
       --reference_genome genome.fa \\
       --annotation_dir annotations/ \\
       --reference_monomers monomers.fa \\
       --family_assignments families.txt \\
       --sample_name Col_9day \\
       --outdir results_Col/
   ```

2. **Automated outputs:**
   - All visualizations in `05_plots/`
   - Statistics in `07_monomer_statistics/`
   - Sequences in `08_monomer_sequences/`

3. **Optional manual analyses:**
   ```bash
   # Spatial patterns
   python bin/analyze_monomer_positions.py \\
       results_Col/02_monomers/monomer_classifications.tsv \\
       results_Col/position_analysis/
   ```

### Workflow 2: Two-Sample Comparison

For comparing Col-0 vs mutant:

1. **Run pipeline for both samples:**
   ```bash
   # Sample 1
   nextflow run main.nf --sample_name Col_9day --outdir results_Col/ [...]

   # Sample 2
   nextflow run main.nf --sample_name atxr56_9day --outdir results_atxr56/ [...]
   ```

2. **Compare samples:**
   ```bash
   python bin/compare_samples.py \\
       results_Col/02_monomers/monomer_classifications.tsv \\
       results_atxr56/02_monomers/monomer_classifications.tsv \\
       "Col-0" \\
       "atxr5/6" \\
       comparison_Col_vs_atxr56/
   ```

3. **Interpret results:**
   - Check `comparison_report.txt` for statistical significance
   - Review `family_comparison.png` for composition changes
   - Examine `transition_comparison.png` for pattern shifts
   - Analyze `heterogeneity_comparison.png` for diversity changes

### Workflow 3: Multi-Sample Time Series

For developmental or treatment time series:

1. **Run pipeline for all timepoints**

2. **Pairwise comparisons:**
   ```bash
   # T0 vs T1
   python bin/compare_samples.py t0.tsv t1.tsv "T0" "T1" compare_T0_T1/

   # T1 vs T2
   python bin/compare_samples.py t1.tsv t2.tsv "T1" "T2" compare_T1_T2/

   # etc.
   ```

3. **Aggregate results:**
   - Collect chi-square p-values across comparisons
   - Track family frequency changes over time
   - Monitor heterogeneity trends

### Workflow 4: Publication-Quality Figures

For generating publication figures:

1. **Run full pipeline** to get all automated plots

2. **Key figures to include:**
   - `05_plots/comprehensive/family_summary.png` - Shows overall composition + transitions
   - `05_plots/ribbon_plots/ribbon_*.png` - Single-molecule examples
   - `07_monomer_statistics/family_composition.png` - Clean bar chart
   - `07_monomer_statistics/transition_matrix.png` - Transition probabilities
   - `07_monomer_statistics/heterogeneity_metrics.png` - Diversity analysis

3. **For comparisons:**
   - `comparison_*/family_comparison.png` - Sample differences
   - `comparison_*/heterogeneity_comparison.png` - Diversity comparison

4. **For spatial analysis:**
   - `position_analysis/positional_preferences.png`
   - `position_analysis/cooccurrence_heatmap.png`

---

## Interpretation Guide

### Family Composition

**What to look for:**
- **Dominant families**: Typically 1-3 families account for 60-80% of monomers
- **Rare families**: Families <1% may be ancestral or under selection
- **Missing families**: Not all 20 families may be present in every sample

**Example interpretation:**
```
Family 1: 40.3% - Dominant, high abundance
Family 7: 27.6% - Second most abundant
Family 11: 7.3% - Minor but significant
```

This suggests a hierarchical organization with two major families (1 and 7) and multiple minor families.

### Transition Patterns

**What to look for:**
- **Homotypic runs** (F1→F1): Long stretches of same family
- **Heterotypic transitions** (F1→F7): Frequent switching
- **Asymmetric transitions**: F1→F7 ≠ F7→F1 suggests directional bias
- **Rare transitions**: Some family pairs never co-occur

**Example interpretation:**
```
F1→F1: 1,445 (45%) - Long homotypic runs
F1→F7: 995 (31%) - Frequent heterotypic
F7→F1: 211 (7%) - Less frequent reverse
```

This indicates:
1. Family 1 forms long homogeneous blocks
2. Transitions to Family 7 are common
3. Reverse transitions are less frequent → asymmetric recombination?

### Heterogeneity Metrics

**Shannon Entropy:**
- **Low (0-1)**: Homogeneous arrays, dominated by one family
- **Medium (1-2)**: Moderate diversity, 2-3 families
- **High (>2)**: Highly diverse arrays, many families

**Simpson Diversity:**
- **Low (0-0.3)**: Low diversity
- **Medium (0.3-0.7)**: Moderate diversity
- **High (>0.7)**: High diversity

**Example interpretation:**
```
Mean Shannon entropy: 1.85
Mean Simpson diversity: 0.68
Mean families per array: 5.2
```

This indicates moderate heterogeneity - arrays contain multiple families but are not completely mixed.

### Clustering Analysis

**Clustering Index Interpretation:**
- **>1.5**: Strongly clustered (non-random aggregation)
- **1.0-1.5**: Weakly clustered
- **0.8-1.0**: Random distribution
- **<0.8**: Dispersed (actively separated)

**Example interpretation:**
```
Family 1: Index = 2.3 (Clustered)
Family 7: Index = 1.8 (Clustered)
Family 11: Index = 0.6 (Dispersed)
```

Families 1 and 7 form blocks, while Family 11 is scattered - suggests different organizational principles.

### Positional Preferences

**Mean Position Interpretation:**
- **<0.3**: Start preference
- **0.4-0.6**: No preference (uniform)
- **>0.7**: End preference

**Example interpretation:**
```
Family 3: Mean = 0.15 (Start enriched)
Family 7: Mean = 0.52 (Uniform)
Family 13: Mean = 0.82 (End enriched)
```

This suggests:
- Family 3 may initiate arrays
- Family 7 is distributed throughout
- Family 13 may terminate arrays

### Co-occurrence Enrichment

**Enrichment Interpretation:**
- **>2.0**: Strongly associated (2x more than expected)
- **1.2-2.0**: Moderately associated
- **0.8-1.2**: Independent
- **<0.8**: Negative association (avoid each other)

**Example interpretation:**
```
F1-F7: Enrichment = 3.5 (Strong association)
F1-F11: Enrichment = 0.4 (Avoidance)
```

Families 1 and 7 frequently co-occur, while 1 and 11 are rarely adjacent - suggests functional incompatibility?

### Comparison Statistics

**Chi-square test (composition):**
- **p < 0.001**: Very strong evidence of difference
- **p < 0.05**: Significant difference
- **p > 0.05**: No significant difference

**t-test (heterogeneity):**
- Same interpretation as chi-square
- Tests if diversity metrics differ between samples

**Example interpretation:**
```
Chi-square p-value: 1.2e-15 → Highly significant
Family 1 difference: +8.3% → Enriched in mutant
Shannon entropy: p = 0.002 → Significantly different
```

This indicates the mutant has significantly altered family composition and increased diversity.

---

## Common Questions

### Q: Why are some monomers unclassified?

**A:** Monomers may be unclassified due to:
1. Low alignment identity (<70% threshold)
2. Highly degraded sequences
3. Non-satellite sequences
4. Chimeric or rearranged sequences

Typical classification rates: 60-85%

### Q: What's a "good" Shannon entropy value?

**A:** It depends on your biological question:
- **Expected homogeneous arrays**: Low entropy (0-1)
- **Expected mixed arrays**: High entropy (>2)
- **Natural centromeres**: Typically 1-2 (moderate mixing)

### Q: How many monomers needed for robust statistics?

**A:**
- **Basic stats**: >1,000 monomers
- **Family comparisons**: >5,000 monomers
- **Spatial analysis**: >10,000 monomers
- **Rare family analysis**: >50,000 monomers

### Q: Can I use these scripts on non-centromeric satellites?

**A:** Yes! The scripts are general-purpose for any tandem repeat with family classifications:
- Telomeric repeats
- Ribosomal DNA arrays
- Other satellite DNAs
- Any classified tandem repeats

Just ensure you have:
1. FasTAN-detected arrays
2. Minimap2 classifications
3. Family assignments

### Q: How do I cite this work?

**A:** CENprofiler manuscript in preparation. For now:
```
CENprofiler v2.0: Comprehensive monomer-level analysis
of centromeric satellites. https://github.com/jacgonisa/CENprofiler
```

---

## Troubleshooting

### Issue: MUSCLE not found for consensus

**Solution:** Install MUSCLE v5:
```bash
# Linux
wget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64
mv muscle5.1.linux_intel64 muscle
chmod +x muscle
sudo mv muscle /usr/local/bin/

# Or use conda
conda install -c bioconda muscle
```

### Issue: Memory errors on large datasets

**Solution:** For very large datasets (>100k monomers):
1. Increase memory limits in nextflow.config
2. Run manual scripts with data subsets
3. Use filtering (e.g., top arrays only)

### Issue: Few families detected

**Possible causes:**
1. Wrong reference monomers
2. Wrong family assignment file
3. Too stringent identity threshold (try --min_identity 60)
4. Limited diversity in sample

### Issue: Plots look cluttered

**Solution:**
- Transition heatmaps: Filter to top N families
- Position plots: Show only families with >50 occurrences
- Co-occurrence: Use top 30 pairs

---

## Advanced Tips

### 1. Filtering for high-quality monomers

```python
import pandas as pd

df = pd.read_csv('monomer_classifications.tsv', sep='\\t')

# Filter by identity
high_quality = df[df['alignment_identity'] > 8000]  # Top identity

# Filter by array quality
good_arrays = df[df['array_quality'] > 800]

# Save filtered
high_quality.to_csv('filtered_classifications.tsv', sep='\\t', index=False)
```

### 2. Custom family grouping

Group related families for analysis:

```python
# Define family groups
family_groups = {
    'Group_A': [1, 7, 11],
    'Group_B': [2, 4, 6],
    'Group_C': [3, 5, 8]
}

# Add group column
def assign_group(family):
    for group, members in family_groups.items():
        if family in members:
            return group
    return 'Other'

df['family_group'] = df['monomer_family'].apply(assign_group)

# Analyze by group instead of individual families
```

### 3. Time-series heatmap

Create a heatmap showing family frequency changes over time:

```python
import matplotlib.pyplot as plt
import seaborn as sns

# Collect data from multiple timepoints
timepoints = ['T0', 'T1', 'T2', 'T3']
data = {}

for tp in timepoints:
    df = pd.read_csv(f'results_{tp}/family_statistics.tsv', sep='\\t')
    data[tp] = df.set_index('family')['percentage']

# Combine
combined = pd.DataFrame(data)

# Plot
plt.figure(figsize=(10, 8))
sns.heatmap(combined, cmap='YlOrRd', annot=True, fmt='.1f')
plt.xlabel('Timepoint')
plt.ylabel('Family')
plt.title('Family Composition Over Time (%)')
plt.tight_layout()
plt.savefig('timeseries_heatmap.png', dpi=300)
```

### 4. Export for Cytoscape network analysis

Export co-occurrence data for network visualization:

```python
df_cooccur = pd.read_csv('cooccurrence.tsv', sep='\\t')

# Filter significant enrichments
sig = df_cooccur[df_cooccur['enrichment'] > 1.5]

# Export edge list
sig[['family1', 'family2', 'enrichment']].to_csv(
    'network_edges.txt', sep='\\t', index=False, header=False
)

# Import into Cytoscape for visualization
```

---

## Future Enhancements

Planned additions:
- [ ] Automated multi-sample comparison workflows
- [ ] Interactive HTML reports with plotly
- [ ] Integration with IGV for visualization
- [ ] Phylogenetic tree integration
- [ ] Sequence motif discovery within families
- [ ] Machine learning classification
- [ ] Long-read phasing analysis

---

## Contact & Support

For questions, issues, or feature requests:
- GitHub Issues: https://github.com/jacgonisa/CENprofiler/issues
- Email: [maintainer email]

---

**Last updated:** 2026-01-12
**Pipeline version:** CENprofiler v2.0
