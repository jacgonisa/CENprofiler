# CENprofiler: Centromeric Satellite & HOR Analysis Pipeline

**Version:** 2.0.0
**Status:** Production Ready - Auto-detection and BAM support

---

## Overview

CENprofiler is a Nextflow pipeline for comprehensive analysis of centromeric satellites and Higher-Order Repeats (HORs) in both reference genomes and long-read sequencing data.

### Key Features

✅ **Two Analysis Modes:**
- **Genome Mode**: Analyze reference genomes for satellite arrays and HORs
- **Read Mode**: Analyze long reads for satellite composition and structural variants

✅ **Modular Architecture:**
- Easy to modify and extend
- Clear separation of concerns
- Flexible parameter configuration

✅ **Gap-Aware HOR Detection:**
- Detects HORs with both criteria: min_copies ≥3 AND monomers_per_unit ≥3
- Respects gaps between consecutive monomers (max_gap threshold)
- Identifies large duplications

✅ **Chromosome-Aware Statistics:**
- Track family and HOR prevalence per chromosome
- Compare satellite landscapes across chromosomes

---

## Quick Start

### Prerequisites

- **Nextflow** (≥21.10.0)
- **FasTAN** (installed at `/home/jg2070/bin/FasTAN`)
- **tanbed** (from alntools, at `/home/jg2070/alntools/tanbed`)
- **minimap2**
- **Python 3** with pandas, BioPython

### Installation

```bash
# Clone repository
git clone https://github.com/yourusername/CENprofiler.git
cd CENprofiler

# Test installation
nextflow run main.nf --help
```

---

## Usage

### Genome Mode

Analyze a reference genome for satellites and HORs:

```bash
nextflow run main.nf \\
  --mode genome \\
  --input genome.fasta \\
  --reference_monomers /path/to/Col-CC-V2-CEN178-representative.fasta \\
  --family_assignments /path/to/itol_manual_phylo_clusters.txt \\
  --outdir results/genome_analysis
```

### Read Mode

Analyze long reads for satellite composition:

```bash
nextflow run main.nf \\
  --mode reads \\
  --input reads.fasta \\
  --reference_monomers /path/to/Col-CC-V2-CEN178-representative.fasta \\
  --family_assignments /path/to/itol_manual_phylo_clusters.txt \\
  --outdir results/read_analysis \\
  --analyze_indels true
```

---

## Pipeline Workflow

### Genome Mode

```
genome.fasta
    ↓
[1] FasTAN (tandem repeat detection)
    ↓
[2] tanbed (convert to BED)
    ↓
[3] Extract Monomers
    ↓
[4] Classify Monomers (minimap2 + family assignment)
    ↓
[5] Detect HORs (gap-aware, dual criteria)
    ↓
[6] Chromosome Statistics
    ↓
[7] Generate Plots
```

### Read Mode

```
reads.fasta
    ↓
[1] FasTAN (tandem repeat detection)
    ↓
[2] tanbed (convert to BED)
    ↓
[3] Extract Monomers per Read
    ↓
[4] Classify Monomers (minimap2 + family assignment)
    ↓
[5] Analyze Indels (optional)
    ↓
[6] Generate Plots
```

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--mode` | Analysis mode: 'genome' or 'reads' |
| `--input` | Input FASTA file |
| `--reference_monomers` | Reference monomer FASTA for classification |
| `--family_assignments` | Family assignment TSV (monomer_id<tab>family_id) |
| `--outdir` | Output directory |

### FasTAN Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--period_min` | 160 | Minimum tandem repeat period (bp) |
| `--period_max` | 200 | Maximum tandem repeat period (bp) |
| `--fastan_threads` | 8 | Threads for FasTAN |

### Classification Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_identity` | 70 | Minimum alignment identity (%) |
| `--minimap2_threads` | 4 | Threads for minimap2 |
| `--min_monomer_length` | 150 | Minimum monomer length (bp) |
| `--max_monomer_length` | 210 | Maximum monomer length (bp) |

### HOR Detection (Genome Mode)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_copies` | 3 | Minimum HOR copies (pattern repeats ≥3×) |
| `--min_monomers` | 3 | Minimum monomers per HOR unit |
| `--max_pattern_length` | 20 | Maximum monomers in HOR pattern |
| `--max_gap` | 500 | Maximum gap between monomers (bp) |
| `--large_dup_threshold` | 40 | Large duplication threshold (kb) |

### Read Mode

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--analyze_indels` | true | Perform indel analysis |
| `--min_array_size` | 5 | Minimum array size (monomers) |

---

## Output Structure

```
results/
├── 01_fastan/
│   ├── *.1aln              # FasTAN alignment output
│   ├── *.bed               # Tandem arrays in BED format
│   └── fastan.log
├── 02_monomers/
│   ├── monomers.fa         # Extracted monomer sequences
│   ├── monomer_info.tsv    # Monomer positions and metadata
│   ├── monomer_classifications.tsv  # Main output!
│   └── classification.log
├── 03_hors/  (genome mode)
│   ├── hors_detected.tsv   # All detected HORs
│   ├── large_duplications.tsv
│   └── hor_detection.log
├── 04_stats/  (genome mode)
│   ├── chromosome_stats.tsv
│   ├── family_by_chromosome.tsv
│   └── hor_by_chromosome.tsv
├── 04_indels/  (read mode)
│   ├── indel_stats.tsv
│   └── deletion_monomers.tsv
├── 05_plots/
│   └── plots/
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    └── execution_trace.txt
```

---

## Output Files

### Monomer Classifications (`monomer_classifications.tsv`)

Primary output with per-monomer information:

| Column | Description |
|--------|-------------|
| `monomer_id` | Unique monomer identifier |
| `seq_id` | Source sequence (chromosome or read) |
| `array_idx` | Array index within sequence |
| `monomer_idx` | Monomer index within array |
| `monomer_start` | Start position (bp) |
| `monomer_end` | End position (bp) |
| `monomer_length` | Length (bp) |
| `array_period` | Tandem array period from FasTAN |
| `array_quality` | FasTAN quality score |
| `best_match` | Best matching reference monomer |
| `alignment_identity` | Alignment identity (%) |
| `mapq` | Mapping quality |
| `monomer_family` | Assigned family (1-20 or NA) |

### HORs Detected (`hors_detected.tsv`)

| Column | Description |
|--------|-------------|
| `seq_id` | Source sequence |
| `hor_start` | HOR start position (bp) |
| `hor_end` | HOR end position (bp) |
| `hor_unit` | HOR pattern (e.g., "1F3-1F3-1F3") |
| `hor_unit_length` | Monomers per unit |
| `hor_copies` | Number of repetitions |
| `total_monomers` | Total monomers in HOR |
| `hor_type` | homHOR or hetHOR |
| `hor_length_bp` | Length in base pairs |
| `hor_length_kb` | Length in kilobases |

---

## Reference Files

The pipeline requires two reference files for monomer classification:

### 1. Representative Monomers FASTA

Example: `Col-CC-V2-CEN178-representative.fasta`

```
>M1000_Chr5_9
ATCGATCGATCG...
>M1001_Chr5_9
ATCGATCGATCG...
```

### 2. Family Assignments TSV

Example: `itol_manual_phylo_clusters.txt`

```
# sequence_name	cluster_id
M1000_Chr5_9	2
M1001_Chr5_9	2
M1002_Chr5_9	8
...
```

**Note:** This file defines the phylogenetic classification of monomers into families (1-20).
The pipeline performs **assignment** (not classification) by mapping query monomers to these pre-classified references.

---

## Examples

### Example 1: Analyze Arabidopsis Genome

```bash
nextflow run main.nf \\
  --mode genome \\
  --input TAIR12/GCA_028009825.2_Col-CC_genomic.fna \\
  --reference_monomers kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta \\
  --family_assignments kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt \\
  --outdir results/arabidopsis_genome \\
  --min_copies 3 \\
  --min_monomers 3 \\
  --max_gap 500
```

### Example 2: Analyze ONT Reads

```bash
nextflow run main.nf \\
  --mode reads \\
  --input reads/col-sorted_reads.fasta \\
  --reference_monomers kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta \\
  --family_assignments kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt \\
  --outdir results/col-sorted_reads \\
  --analyze_indels true
```

---

## Development Status

### ✅ Completed

- Modular Nextflow backbone
- FasTAN integration
- Monomer extraction and classification
- Gap-aware HOR detection
- Chromosome-aware statistics
- Dual-mode architecture (genome/reads)

### 🚧 To Be Implemented

- Visualization modules (genome plots, HOR schematics, etc.)
- Indel analysis (integrate existing scripts)
- Read-level HOR detection
- MultiQC-style HTML reports
- Comprehensive testing suite

### 📝 Integration Needed

The following existing scripts should be integrated:

**Genome Mode Plots:**
- `plot_monomer_level_genome_wide.py`
- `plot_large_duplications_detail.py`
- `plot_large_duplications_overview.py`
- `plot_monomer_level_schematics.py`
- `analyze_monomer_enrichment_monomer_level.py`

**Read Mode Analysis:**
- `analyze_deletion_monomers.py`
- `large_scale_indel_analysis.py`
- `visualize_indel_families_v2.py`

---

## Troubleshooting

### FasTAN not found

```
Error: FasTAN not found at /home/jg2070/bin/FasTAN
```

Update the path in `nextflow.config`:
```groovy
params.fastan_bin = '/path/to/FasTAN'
```

### No monomers classified

Check:
1. Minimum identity threshold (--min_identity, default 70%)
2. Reference monomers file exists and is properly formatted
3. Family assignments file matches monomer IDs

### No HORs detected

Check:
1. Classification rate (need classified monomers)
2. HOR parameters (min_copies, min_monomers)
3. Gap threshold (max_gap) - may be too strict

---

## Citation

If you use CENprofiler, please cite:

- **FasTAN**: Myers et al. (TBD)
- **alntools**: Durbin et al. (TBD)
- **minimap2**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.

---

## Contact

For questions or issues:
- GitHub Issues: https://github.com/yourusername/CENprofiler/issues
- Email: jg2070@cam.ac.uk

---

## License

MIT License (or specify your preferred license)

---

**CENprofiler** - Comprehensive Centromeric Satellite Analysis
Version 1.0.0 | Development Version
