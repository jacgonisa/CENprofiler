# CENprofiler - Quick Start Guide

## 🚀 Get Started in 5 Minutes

### Step 1: Test the Pipeline

```bash
cd /mnt/ssd-8tb/atrx_china/CENprofiler

# Run automated tests
./test_pipeline.sh
```

This will test both genome and read modes using existing validated data.

---

### Step 2: Run Your First Analysis

#### Option A: Genome Mode

```bash
nextflow run main.nf \
  --mode genome \
  --input /path/to/your/genome.fasta \
  --reference_monomers /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta \
  --family_assignments /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt \
  --outdir results/my_genome
```

#### Option B: Read Mode

```bash
nextflow run main.nf \
  --mode reads \
  --input /path/to/your/reads.fasta \
  --reference_monomers /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta \
  --family_assignments /mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt \
  --outdir results/my_reads
```

---

### Step 3: Check Results

```bash
# Main output: monomer classifications
head results/my_genome/02_monomers/monomer_classifications.tsv

# HORs detected (genome mode)
head results/my_genome/03_hors/hors_detected.tsv

# Statistics summary
cat results/my_genome/02_monomers/classification.log
cat results/my_genome/03_hors/hor_detection.log
```

---

## 📊 Key Output Files

### 1. Monomer Classifications
**Location**: `results/02_monomers/monomer_classifications.tsv`

Contains per-monomer information:
- Position coordinates
- Family assignment
- Alignment identity
- Quality metrics

### 2. HORs Detected (Genome Mode)
**Location**: `results/03_hors/hors_detected.tsv`

Contains detected HORs:
- HOR pattern (e.g., "3F3")
- Number of copies
- Genomic coordinates
- homHOR vs hetHOR classification

### 3. Chromosome Statistics (Genome Mode)
**Location**: `results/04_stats/chromosome_stats.tsv`

Per-chromosome:
- Total monomers
- Family distribution
- HOR counts

---

## ⚙️ Common Parameters

### Adjust HOR Detection

```bash
nextflow run main.nf \
  --mode genome \
  --input genome.fasta \
  ... \
  --min_copies 5 \         # Require 5+ copies (default: 3)
  --min_monomers 3 \       # Require 3+ monomers per unit (default: 3)
  --max_gap 1000           # Allow 1kb gaps (default: 500bp)
```

### Adjust Classification

```bash
nextflow run main.nf \
  --mode genome \
  --input genome.fasta \
  ... \
  --min_identity 80 \      # Require 80% identity (default: 70%)
  --period_min 170 \       # Min period 170bp (default: 160)
  --period_max 190         # Max period 190bp (default: 200)
```

---

## 🔧 Troubleshooting

### No Monomers Classified

**Check:**
1. Classification threshold: `--min_identity` (try lowering to 60%)
2. Reference files exist and are readable
3. FasTAN detected arrays (check `01_fastan/*.log`)

### No HORs Detected

**Check:**
1. Monomers were classified (check `02_monomers/*.log`)
2. HOR parameters not too strict:
   - `--min_copies` (try 2)
   - `--min_monomers` (try 2)
   - `--max_gap` (try 1000)

### Pipeline Crashes

**Check:**
1. Nextflow version: `nextflow -version` (need ≥21.10.0)
2. FasTAN installed: `ls -lh /home/jg2070/bin/FasTAN`
3. Disk space: `df -h`
4. Memory: `free -h`

---

## 📁 Where is Everything?

```
CENprofiler/
├── main.nf                 # Main pipeline
├── nextflow.config         # Configuration
├── README.md               # Full documentation
├── QUICK_START.md          # This file
├── PIPELINE_SUMMARY.md     # Development summary
├── test_pipeline.sh        # Test suite
├── workflows/              # Genome/Read workflows
├── modules/                # Individual processes
└── bin/                    # Helper scripts
```

---

## 📚 Next Steps

1. **Read full documentation**: `README.md`
2. **Review pipeline design**: `PIPELINE_SUMMARY.md`
3. **Run tests**: `./test_pipeline.sh`
4. **Analyze your data**: Use examples above
5. **Adjust parameters**: Based on results

---

## 🆘 Need Help?

- **Documentation**: See `README.md`
- **Test data**: Run `./test_pipeline.sh`
- **Parameters**: `nextflow run main.nf --help`
- **Logs**: Check `results/*/pipeline_info/`

---

**Happy profiling! 🧬**
