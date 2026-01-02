# CENprofiler Pipeline - Development Summary

**Created:** January 1, 2026
**Status:** Modular backbone completed ✅
**Next Steps:** Integration of plotting and indel analysis scripts

---

## 🎯 What Has Been Created

### Core Pipeline Components

✅ **Main Workflow** (`main.nf`)
- Entry point for the pipeline
- Parameter validation
- Help message system
- Workflow routing (genome/read modes)

✅ **Configuration** (`nextflow.config`)
- All parameters with sensible defaults
- Process-specific resource requirements
- Execution profiles (local, cluster, debug)
- Pipeline manifest and reporting

✅ **Workflow Definitions**
- `workflows/genome_mode.nf` - Reference genome analysis
- `workflows/read_mode.nf` - Long read analysis

✅ **Process Modules** (9 modules)
1. `fastan.nf` - Tandem repeat detection
2. `tanbed.nf` - Convert FasTAN output to BED
3. `extract_monomers.nf` - Extract individual monomers
4. `classify_monomers.nf` - Minimap2 + family assignment
5. `detect_hors.nf` - Gap-aware HOR detection
6. `chromosome_stats.nf` - Chromosome-level statistics
7. `analyze_indels.nf` - Indel analysis (placeholder)
8. `genome_plots.nf` - Genome visualizations (placeholder)
9. `read_plots.nf` - Read visualizations (placeholder)

✅ **Testing Framework**
- `test_pipeline.sh` - Automated test script
- Tests both genome and read modes
- Uses existing test data
- Provides detailed output summary

✅ **Documentation**
- `README.md` - Comprehensive user guide
- Usage examples
- Parameter descriptions
- Output file specifications

---

## 🏗️ Architecture

### Modular Design

```
CENprofiler/
├── main.nf                    # Entry point
├── nextflow.config            # Configuration
├── workflows/                 # Mode-specific workflows
│   ├── genome_mode.nf
│   └── read_mode.nf
├── modules/                   # Individual processes
│   ├── fastan.nf
│   ├── tanbed.nf
│   ├── extract_monomers.nf
│   ├── classify_monomers.nf
│   ├── detect_hors.nf
│   ├── chromosome_stats.nf
│   ├── analyze_indels.nf
│   ├── genome_plots.nf
│   └── read_plots.nf
├── bin/                       # Helper scripts
│   └── detect_hors_monomer_level.py
└── test_pipeline.sh           # Test suite
```

### Key Design Principles

1. **Modularity**: Each step is a separate module - easy to modify
2. **Flexibility**: Clear parameter system - easy to configure
3. **Extensibility**: Placeholder modules for future additions
4. **Testability**: Test script using real data

---

## ✅ Implemented Features

### Core Functionality

- [x] Dual-mode operation (genome/reads)
- [x] FasTAN integration
- [x] Monomer extraction and filtering
- [x] Minimap2-based classification
- [x] Family assignment from phylogenetic clusters
- [x] Gap-aware HOR detection
- [x] Chromosome-aware statistics
- [x] Configurable parameters

### HOR Detection Algorithm

- [x] Monomer-level detection (not RLE)
- [x] Dual criteria: min_copies ≥3 AND monomers_per_unit ≥3
- [x] Gap checking (max_gap threshold)
- [x] Prefers shorter patterns
- [x] Large duplication identification

---

## 🚧 To Be Implemented

### High Priority

1. **Visualization Modules**
   - Integrate existing plotting scripts:
     - `plot_monomer_level_genome_wide.py`
     - `plot_large_duplications_detail.py`
     - `plot_large_duplications_overview.py`
     - `plot_monomer_level_schematics.py`
     - `analyze_monomer_enrichment_monomer_level.py`

2. **Indel Analysis**
   - Integrate existing scripts:
     - `analyze_deletion_monomers.py`
     - `large_scale_indel_analysis.py`
     - `visualize_indel_families_v2.py`
   - Implement metrics:
     - Indel size distribution per family
     - Insertion/deletion ratio per family
     - HOR-associated indels
     - Transition hotspot analysis

### Medium Priority

3. **Read-level HOR Detection**
   - Optional HOR detection in read mode
   - Compare read HORs to reference HORs

4. **Enhanced Reporting**
   - MultiQC-style HTML reports
   - Interactive plots
   - Summary statistics dashboard

5. **Advanced Features**
   - Family co-occurrence analysis
   - Transition pattern statistics
   - Phylogenetic integration

---

## 🧪 Testing

### Running Tests

```bash
cd /mnt/ssd-8tb/atrx_china/CENprofiler
./test_pipeline.sh
```

### Test Data

Uses existing validated datasets:
- **Genome Mode**: `/mnt/ssd-8tb/atrx_china/CEN178profiler/v3_genome_analysis_CORRECT/centromeric_regions.fa`
- **Read Mode**: `/mnt/ssd-8tb/atrx_china/CEN178profiler/results_v2_test/test_sv_reads.fa`

### Expected Results

**Genome Mode:**
- Extracts monomers from centromeric regions
- Classifies to 20 families
- Detects HORs with gap awareness
- Generates chromosome statistics

**Read Mode:**
- Extracts monomers from tandem arrays in reads
- Classifies per monomer
- Identifies family transitions
- (Optional) Analyzes indels

---

## 📦 Dependencies

### Required Tools

- **Nextflow** (≥21.10.0)
- **FasTAN** (`/home/jg2070/bin/FasTAN`)
- **tanbed** (`/home/jg2070/alntools/tanbed`)
- **minimap2**
- **Python 3** with:
  - pandas
  - BioPython
  - numpy
  - matplotlib (for plotting modules)

### Reference Files

- **Representative Monomers**: `Col-CC-V2-CEN178-representative.fasta` (5,597 monomers)
- **Family Assignments**: `itol_manual_phylo_clusters.txt` (20 families)

---

## 🔧 Usage Examples

### Genome Mode

```bash
nextflow run main.nf \
  --mode genome \
  --input genome.fasta \
  --reference_monomers /path/to/representatives.fasta \
  --family_assignments /path/to/families.txt \
  --outdir results/genome \
  --min_copies 3 \
  --min_monomers 3 \
  --max_gap 500
```

### Read Mode

```bash
nextflow run main.nf \
  --mode reads \
  --input reads.fasta \
  --reference_monomers /path/to/representatives.fasta \
  --family_assignments /path/to/families.txt \
  --outdir results/reads \
  --analyze_indels true
```

---

## 📊 Output Structure

```
results/
├── 01_fastan/
│   ├── *.1aln               # FasTAN output
│   ├── *.bed                # Tandem arrays
│   └── *.log
├── 02_monomers/
│   ├── monomers.fa          # Extracted monomers
│   ├── monomer_info.tsv     # Positions/metadata
│   ├── monomer_classifications.tsv  # ⭐ MAIN OUTPUT
│   └── *.log
├── 03_hors/  (genome mode)
│   ├── hors_detected.tsv    # ⭐ HORs
│   ├── large_duplications.tsv
│   └── *.log
├── 04_stats/  (genome mode)
│   ├── chromosome_stats.tsv # ⭐ Chromosome-level
│   ├── family_by_chromosome.tsv
│   └── hor_by_chromosome.tsv
├── 04_indels/  (read mode)
│   ├── indel_stats.tsv
│   └── deletion_monomers.tsv
├── 05_plots/
│   └── plots/               # Visualizations
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    └── execution_trace.txt
```

---

## 🎨 Chromosome-Aware Statistics

The pipeline tracks satellite and HOR prevalence per chromosome:

**Family Distribution:**
```
chromosome  monomer_family  count
Chr1        1               1234
Chr1        3               567
Chr2        1               890
...
```

**HOR Distribution:**
```
chromosome  hor_unit    count
Chr1        3F1         51
Chr1        3F4         47
Chr4        3F3         42
...
```

This enables:
- Chromosome-specific satellite landscapes
- Cross-chromosome comparisons
- Centromere-specific enrichment analysis

---

## 📝 Key Parameters

### HOR Detection (Dual Criteria)

```groovy
params {
    // BOTH requirements must be met:
    min_copies = 3      // Pattern repeats ≥3 times
    min_monomers = 3    // Pattern has ≥3 monomers

    // Gap awareness:
    max_gap = 500       // Gaps >500bp break HORs

    // Pattern search:
    max_pattern_length = 20  // Search up to 20-mer patterns
}
```

**Example:**
- `3F3 × 356 copies` = **PASS** ✓ (3 monomers, 356 copies)
- `1F3 × 1068 copies` = **FAIL** ✗ (1 monomer, 1068 copies)
- `5F1 × 2 copies` = **FAIL** ✗ (5 monomers, 2 copies)

---

## 🚀 Next Development Steps

### Immediate (This Week)

1. **Integrate plotting modules**
   - Copy existing scripts to `bin/`
   - Update `genome_plots.nf` and `read_plots.nf`
   - Test with real data

2. **Integrate indel analysis**
   - Copy existing scripts to `bin/`
   - Update `analyze_indels.nf`
   - Define indel metrics

### Short-term (This Month)

3. **Enhanced visualization**
   - Interactive HTML plots
   - MultiQC integration
   - Summary dashboards

4. **Documentation**
   - Tutorial with example data
   - Parameter tuning guide
   - Troubleshooting section

### Long-term

5. **Advanced features**
   - Read-level HOR detection
   - Comparative analysis (multi-sample)
   - Automated quality control
   - Docker/Singularity containers

---

## 💡 Design Decisions

### Why This Architecture?

1. **Modular processes**: Each step can be modified independently
2. **Dual workflows**: Clear separation of genome vs read logic
3. **Flexible configuration**: Easy to adjust parameters without code changes
4. **Test-driven**: Test script ensures core functionality works

### Why Not Fully Implemented Yet?

Following user requirements:
> "DONT WORRY TOO MUCH NOW ABOUT CONTAINERISATION, JUST MAKE A NICE MODULAR THING THAT SERVES AS BACKBONE BECAUSE THE PIPELINE HAS TO CHANGE STILL A LOT"

- **Focus**: Modular backbone ✅
- **Flexibility**: Easy to modify ✅
- **Integration**: Ready to add existing scripts ✅
- **Future-proof**: Can evolve as requirements change ✅

---

## 🔍 Integration Points for Existing Code

### From CEN178profiler v3 (Genome Mode)

**Location**: `/mnt/ssd-8tb/atrx_china/CEN178profiler/v3_genome_analysis_CORRECT/exact_matching_analysis/`

**Scripts to integrate:**
```python
# Already integrated:
detect_hors_monomer_level.py  # → bin/ ✅

# To integrate:
plot_monomer_level_genome_wide.py          # → genome_plots module
plot_large_duplications_detail.py          # → genome_plots module
plot_large_duplications_overview.py        # → genome_plots module
plot_monomer_level_schematics.py           # → genome_plots module
analyze_monomer_enrichment_monomer_level.py # → genome_plots module
```

### From CEN178profiler v2 (Read Mode)

**Location**: `/mnt/ssd-8tb/atrx_china/CEN178profiler/results_v2_test/scripts/`

**Scripts to integrate:**
```python
analyze_deletion_monomers.py           # → analyze_indels module
large_scale_indel_analysis.py          # → analyze_indels module
visualize_indel_families_v2.py         # → read_plots module
```

---

## ✨ Summary

### What We Have

✅ **Solid foundation**: Modular, flexible Nextflow pipeline
✅ **Core functionality**: FasTAN → Extract → Classify → Detect HORs
✅ **Chromosome awareness**: Track statistics per chromosome
✅ **Gap-aware HORs**: Properly implemented with dual criteria
✅ **Testing framework**: Automated tests with real data
✅ **Documentation**: Comprehensive README and examples

### What's Next

🔨 **Integrate plotting**: Add your existing visualization scripts
🔨 **Integrate indel analysis**: Add your existing analysis scripts
🔨 **Test thoroughly**: Run on full datasets
🔨 **Iterate**: Modify based on results and feedback

---

**Ready to use as backbone for further development! 🚀**

The pipeline is now at a stage where:
- It can be tested with real data
- Existing scripts can be integrated module-by-module
- Parameters can be tuned based on results
- New features can be added without breaking existing functionality
