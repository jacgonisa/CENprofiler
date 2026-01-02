# CENprofiler v2.0 Architecture

## Major Redesign: Auto-Detection Mode

**Date**: 2026-01-02
**Version**: 2.0
**Status**: Implemented

---

## Overview

CENprofiler v2.0 introduces **automatic mode detection** based on input parameters, eliminating the artificial genome/reads mode split and enabling proper indel analysis with BAM files.

---

## Auto-Detection Logic

```
IF --alignment provided:
    → READ MODE with INDEL ANALYSIS
    Required inputs:
      - --alignment (BAM/SAM file)
      - --reference_genome (reference FASTA)
      - --annotation_dir (genomic regions)
      - --reference_monomers
      - --family_assignments

ELSE:
    → GENOME MODE
    Required inputs:
      - --input (genome FASTA)
      - --reference_monomers
      - --family_assignments
```

---

## Workflows

### 1. GENOME MODE (unchanged)

**Triggered when**: `--alignment` is NOT provided

**Input**: Reference genome FASTA

**Steps**:
1. FasTAN → detect tandem arrays
2. TANBED → convert to BED
3. EXTRACT_MONOMERS → individual monomers
4. CLASSIFY_MONOMERS → family assignment
5. DETECT_HORS → gap-aware HOR detection
6. CHROMOSOME_STATS → per-chromosome statistics
7. SATELLITE_PLOTS → monomer visualizations
8. HOR_PLOTS → HOR schematics

**Outputs**:
- Monomer classifications
- HOR catalog with coordinates
- Large duplications (≥40kb)
- Chromosome statistics
- Visualization plots

---

### 2. READ MODE with INDEL ANALYSIS (NEW!)

**Triggered when**: `--alignment` is provided

**Input**: BAM alignment file

**Steps**:
1. **LOAD_GENOMIC_REGIONS** → Load annotations
   - Centromere coordinates
   - Pericentromere coordinates
   - 5S and 45S rDNA regions
   - Handle chromosome name mapping (CP → Chr)

2. **EXTRACT_READS_FROM_BAM** → Extract reads with large indels
   - Parse CIGAR strings for insertions/deletions ≥ min_indel_size
   - Classify indels by genomic region
   - Calculate CEN178 multiple information
   - Output: reads FASTA + indel catalog TSV

3. **FASTAN** → Detect tandem arrays in extracted reads

4. **TANBED** → Convert to BED format

5. **EXTRACT_MONOMERS** → Individual monomers from arrays

6. **CLASSIFY_MONOMERS** → Assign families using minimap2

7. **READ_PLOTS** → Generate visualizations
   - Family distribution
   - Per-read statistics
   - Indel size distributions

**Outputs**:
- Extracted reads FASTA
- Indel catalog with region classifications
- Monomer classifications
- Family-specific indel statistics
- Visualization plots

---

## New Modules

### `LOAD_GENOMIC_REGIONS`

**Purpose**: Load genomic region annotations from BED files

**Script**: `bin/load_genomic_regions.py`

**Input**:
- Annotation directory containing:
  - `centromere.bed`
  - `pericentromere_clean.bed`
  - `5s_rdna_regions.bed`
  - `45s_rdna_regions.bed`

**Output**:
- `genomic_regions.tsv` - Unified regions file

**Features**:
- Chromosome name mapping (CP accessions → Chr format)
- Priority-based region classification

---

### `EXTRACT_READS_FROM_BAM`

**Purpose**: Extract reads containing large indels from BAM alignment

**Script**: `bin/extract_reads_from_bam.py`

**Input**:
- BAM file with index
- Genomic regions TSV
- Minimum indel size (default: 100bp)

**Output**:
- `{sample}_reads.fa` - FASTA of extracted reads
- `{sample}_indel_catalog.tsv` - Complete indel catalog
- `{sample}_stats.txt` - Extraction statistics

**Process**:
1. Scan BAM chromosome by chromosome
2. Parse CIGAR strings for insertions (op=1) and deletions (op=2)
3. Filter indels ≥ min_indel_size
4. Classify indel position by region (priority: 5s_rdna > 45s_rdna > centromere > pericentromere > arms)
5. Calculate CEN178 multiple metrics:
   - `multiple_of_178`: Boolean (length % 178 == 0)
   - `closest_178_multiple`: Nearest 178bp multiple
   - `distance_to_178_multiple`: Distance in bp
6. Extract unique reads with large indels

**Performance**: Memory-efficient streaming, processes millions of reads

---

## New Parameters

### Read Mode with BAM

```groovy
params {
    // Triggers read mode with indel analysis
    alignment               = null  // BAM/SAM file path
    reference_genome        = null  // Reference genome FASTA
    annotation_dir          = null  // Annotations directory
    sample_name             = 'sample'  // Sample identifier

    // Indel extraction
    min_indel_size          = 100   // Minimum indel size (bp)
}
```

---

## Indel Catalog Format

```
read_id  chromosome  ref_pos  read_pos  type  size  region  mapping_quality  multiple_of_178  closest_178_multiple  distance_to_178_multiple
```

**Columns**:
- `read_id`: Read identifier
- `chromosome`: Chromosome name (e.g., Chr1)
- `ref_pos`: Reference position
- `read_pos`: Position within read
- `type`: Insertion or Deletion
- `size`: Indel size in bp
- `region`: Genomic region (centromere, pericentromere, arms, 5s_rdna, 45s_rdna, other)
- `mapping_quality`: MAPQ score
- `multiple_of_178`: Boolean, exact CEN178 multiple
- `closest_178_multiple`: Ratio (size / 178)
- `distance_to_178_multiple`: Distance to nearest 178bp multiple

---

## Region Classification

**Priority order** (highest to lowest):
1. **5s_rdna** - 5S ribosomal DNA arrays
2. **45s_rdna** - 45S ribosomal DNA arrays
3. **centromere** - Core centromeric regions
4. **pericentromere** - Flanking regions (±500kb from centromere)
5. **arms** - Chromosome arms (everything else)
6. **other** - Non-standard chromosomes

---

## Example Commands

### 1. Genome Mode (Auto-Detected)

```bash
nextflow run main.nf \
    --input /mnt/ssd-8tb/atrx_china/TAIR12/GCA_028009825.2_Col-CC_genomic.fna \
    --reference_monomers /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta \
    --family_assignments /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt \
    --outdir results_genome/
```

### 2. Read Mode with Indels (Auto-Detected)

```bash
nextflow run main.nf \
    --input dummy.fa \
    --alignment /mnt/ssd-8tb/atrx_china/tair12_indel_comparison/results/mapping/Col_9day.bam \
    --reference_genome /mnt/ssd-8tb/atrx_china/TAIR12/GCA_028009825.2_Col-CC_genomic.fna \
    --annotation_dir /mnt/ssd-8tb/atrx_china/TAIR12/curated_anno \
    --reference_monomers /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta \
    --family_assignments /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt \
    --sample_name Col_9day \
    --min_indel_size 100 \
    --outdir results_Col_9day/
```

---

## Migration from v1.0

### What Changed

**Removed**:
- `--mode genome|reads` parameter (auto-detected now)
- Simple read mode without indel analysis

**Added**:
- `--alignment` parameter (triggers read mode with indels)
- `--reference_genome` parameter (required with --alignment)
- `--annotation_dir` parameter (genomic regions)
- `--sample_name` parameter
- `--min_indel_size` parameter (default: 100bp)

**Preserved**:
- All genome mode functionality
- All FasTAN parameters
- All classification parameters
- All HOR detection parameters

### Backward Compatibility

**Genome mode**: Fully compatible - same command structure, just remove `--mode genome`

**Read mode**: Requires migration to new BAM-based workflow

---

## Performance Considerations

### Memory Requirements

- **GENOME MODE**: 4-32GB (depends on genome size)
- **READ MODE**: 32GB+ for large BAM files (54-62GB BAMs tested)
- **EXTRACT_READS_FROM_BAM**: Memory-efficient streaming (doesn't load all reads)

### Execution Time

- **Genome mode** (Arabidopsis centromeres, ~5Mb): ~5-10 minutes
- **Read mode** (54GB BAM, ~millions of reads): Variable
  - BAM scanning: Hours (depends on BAM size)
  - FasTAN on extracted reads: Minutes to hours
  - Classification: Hours (depends on read count)

---

## Future Enhancements

### Planned Features
- [ ] Match monomers to specific indels (monomer-indel mapping)
- [ ] Family-specific indel visualizations
- [ ] Per-read ribbon plots with indel annotations
- [ ] Statistical analysis of family enrichment in indels
- [ ] Comparison between samples (Col vs atxr56)

### Under Consideration
- [ ] Containerization (Docker/Singularity)
- [ ] Parallelization of BAM extraction
- [ ] Streaming indel analysis (avoid writing intermediate FASTA)
- [ ] Interactive HTML reports

---

## Technical Notes

### Chromosome Name Mapping

The pipeline handles mismatches between GenBank accessions and standard chromosome names:

```python
CHROM_NAME_MAP = {
    'CP116280.1': 'Chr1',
    'CP116281.2': 'Chr2',
    'CP116282.1': 'Chr3',
    'CP116283.2': 'Chr4',
    'CP116284.1': 'Chr5',
}
```

### CIGAR String Parsing

Indels are extracted from CIGAR tuples:
- **Match (op=0)**: Advance both ref_pos and read_pos
- **Insertion (op=1)**: Advance read_pos only
- **Deletion (op=2)**: Advance ref_pos only
- **Soft clip (op=4)**: Advance read_pos only

### CEN178 Multiple Detection

```python
multiple_of_178 = (length % 178 == 0)
closest_178 = round(length / 178, 2)
closest_multiple = round(length / 178) * 178
distance = abs(length - closest_multiple)
```

This helps identify indels that are exact multiples of the CEN178 monomer size (~178bp), suggesting repeat insertion/deletion events.

---

## Credits

**Architecture Design**: Claude Code (Anthropic)
**Biological Context**: J. Gonzalez-Sanchez
**Testing Data**: Arabidopsis thaliana Col-CC and atxr56 mutant

---

**Generated**: 2026-01-02
**Version**: 2.0
**Status**: Production-ready (pending testing)
