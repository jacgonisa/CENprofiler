# Quick Start: Single Molecule Analysis

This guide will help you quickly analyze a single centromeric molecule using CENprofiler.

## Prerequisites

- Nextflow installed
- Python 3.8+
- Required Python packages: pandas, matplotlib, numpy, biopython

## 3-Step Workflow

### Step 1: Prepare Your Data

Create a FASTA file with your centromeric sequence:

```bash
# Your sequence can come from:
# - Long-read HiFi sequencing (PacBio)
# - ONT (Oxford Nanopore) reads
# - Assembled centromeric contigs
# - PCR amplicons

cat > my_centromere.fa << 'EOF'
>my_centromere_chr4
ATCGATCG...your_sequence_here...GCTAGCTA
EOF
```

### Step 2: Run the Pipeline

```bash
nextflow run main.nf \
    --mode genome \
    --input my_centromere.fa \
    --reference_monomers assets/Col-CC-V2-CEN178-representative.fasta \
    --family_assignments assets/itol_manual_phylo_clusters.txt \
    --outdir results_my_centromere/ \
    --fastan_threads 4
```

**Expected runtime:** 5-15 minutes for a single molecule (1-5 Mb)

### Step 3: View Results

```bash
# Ribbon plots showing indel patterns
ls results_my_centromere/04_visualizations/indel_families_*.png

# Summary statistics
cat results_my_centromere/03_indels/indel_statistics.txt

# Monomer classifications
head results_my_centromere/02_monomers/monomer_classifications.tsv
```

## Understanding the Output

### Ribbon Plot Components

The ribbon plot (e.g., `indel_families_0_abc123.png`) shows:

1. **Horizontal colored bars**: Individual monomers
   - Color indicates monomer family (F1-F20)
   - Width indicates monomer length (~178bp)

2. **Vertical lines**: Indel positions
   - Different colors = different indel families
   - Line position = location within monomer
   - Multiple lines = multiple indels at that position

3. **Bottom axis**: Genomic position in kilobases (kb)

4. **Statistics box**: Summary of detected patterns

### Key Output Files

| File | Description |
|------|-------------|
| `monomer_classifications.tsv` | All detected monomers with family assignments |
| `indel_classifications.tsv` | All indels grouped into families |
| `indel_families_*.png` | Visual ribbon plots (one per tandem array) |
| `indel_statistics.txt` | Summary statistics |

## Example Output Interpretation

```
Array 0: 245 monomers
  - 180 F3 monomers (73%)
  - 45 F4 monomers (18%)
  - 20 F1 monomers (8%)

Indels detected: 67
  - 45 deletions (67%)
  - 22 insertions (33%)

Indel families: 12
  - DEL_F3_1: 23 occurrences (most common deletion in F3)
  - INS_F3_1: 15 occurrences (most common insertion in F3)
```

This tells you:
- The array is dominated by F3 monomers (73%)
- Deletions are more common than insertions
- Certain indel patterns repeat frequently (potential functional significance)

## Common Use Cases

### 1. Compare Wild-Type vs. Mutant

```bash
# Process wild-type
nextflow run main.nf --input wildtype.fa --outdir results_wt/

# Process mutant
nextflow run main.nf --input mutant.fa --outdir results_mutant/

# Compare indel patterns
python bin/compare_indel_patterns.py \
    results_wt/03_indels/indel_classifications.tsv \
    results_mutant/03_indels/indel_classifications.tsv \
    comparison_output/
```

### 2. Analyze Multiple Replicates

```bash
# Create a batch script
for replicate in rep1.fa rep2.fa rep3.fa; do
    nextflow run main.nf \
        --input $replicate \
        --outdir results_$(basename $replicate .fa)/
done

# Aggregate results
python bin/aggregate_replicates.py results_rep*/
```

### 3. Focus on Specific Monomer Family

```bash
# Generate ribbon plot for F3 monomers only
python bin/visualize_indel_families_v2.py \
    results/02_monomers/monomer_classifications.tsv \
    results/03_indels/indel_classifications.tsv \
    results/04_visualizations/ \
    --filter-family 3
```

## Troubleshooting

### Problem: No monomers detected

**Solution:**
```bash
# Check if your sequence is centromeric
# Centromeric regions should have ~178bp repeats
# Try BLAST against known centromeric sequences first

# Or adjust detection parameters
nextflow run main.nf ... --fastan_min_length 150 --fastan_max_length 200
```

### Problem: Many "unknown" families

**Solution:**
```bash
# Your sequence may be from a different species
# You need species-specific reference monomers

# Check reference compatibility
head assets/Col-CC-V2-CEN178-representative.fasta

# Consider creating custom reference for your species
```

### Problem: Pipeline fails with memory error

**Solution:**
```bash
# Reduce thread counts
nextflow run main.nf ... --fastan_threads 2 --minimap2_threads 1

# Or split large input into smaller chunks
seqkit split -s 5000000 large_centromere.fa  # Split into 5Mb chunks
```

## Next Steps

- See [WORKFLOW_RIBBON_PLOTS.md](WORKFLOW_RIBBON_PLOTS.md) for detailed explanation
- See [README.md](../README.md) for full pipeline documentation
- See [EXAMPLES.md](EXAMPLES.md) for more analysis examples

## Questions?

Open an issue at: https://github.com/jacgonisa/CENprofiler/issues
