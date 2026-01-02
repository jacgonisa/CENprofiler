# CENprofiler Workflow Diagrams

## Overview Schematic

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         CENprofiler Pipeline                             │
│                   Centromeric Satellite & HOR Analysis                   │
└─────────────────────────────────────────────────────────────────────────┘

                              INPUT MODE
                                  ↓
                    ┌─────────────┴─────────────┐
                    │                           │
              GENOME MODE                  READ MODE
                    │                           │
                    ↓                           ↓
```

---

## Genome Mode Workflow

```
┌───────────────────────────────────────────────────────────────────────────┐
│                           GENOME MODE                                      │
│            Reference Genome → Satellites → HORs → Statistics              │
└───────────────────────────────────────────────────────────────────────────┘

INPUT: genome.fasta
   │
   │  ┌─────────────────────────────────────────────────────────────┐
   └─→│ [1] FasTAN                                                   │
      │     Tandem Repeat Detection                                 │
      │     • Detects tandem arrays in sequence                     │
      │     • Output: .1aln format                                  │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  *.1aln
      ┌─────────────────────────────────────────────────────────────┐
      │ [2] tanbed                                                   │
      │     Convert to BED Format                                   │
      │     • Parses FasTAN alignment output                        │
      │     • Creates BED with array coordinates                    │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  arrays.bed
      ┌─────────────────────────────────────────────────────────────┐
      │ [3] EXTRACT_MONOMERS                                         │
      │     Extract Individual Monomers                             │
      │     • Filter by period (160-200bp for CEN178)              │
      │     • Extract each repeat unit                              │
      │     • Track positions and metadata                          │
      │     Output: monomers.fa + monomer_info.tsv                  │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  monomers.fa
      ┌─────────────────────────────────────────────────────────────┐
      │ [4] CLASSIFY_MONOMERS                                        │
      │     Family Assignment                                        │
      │     • minimap2: align to reference monomers                 │
      │     • Assign family from phylogenetic clusters              │
      │     • Min identity threshold (default 70%)                  │
      │     Output: monomer_classifications.tsv ⭐                   │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  classifications.tsv
      ┌─────────────────────────────────────────────────────────────┐
      │ [5] DETECT_HORS                                              │
      │     Gap-Aware HOR Detection                                 │
      │     • Find repeating patterns (min 3 monomers)              │
      │     • Require ≥3 copies                                      │
      │     • Check gaps (≤500bp)                                    │
      │     • Classify: homHOR vs hetHOR                            │
      │     Output: hors_detected.tsv ⭐                             │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  hors.tsv
      ┌─────────────────────────────────────────────────────────────┐
      │ [6] CHROMOSOME_STATS                                         │
      │     Chromosome-Level Statistics                             │
      │     • Family distribution per chromosome                    │
      │     • HOR prevalence per chromosome                         │
      │     • Cross-chromosome comparisons                          │
      │     Output: chromosome_stats.tsv                            │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  stats.tsv
      ┌─────────────────────────────────────────────────────────────┐
      │ [7] GENOME_PLOTS                                             │
      │     Visualizations                                          │
      │     • Genome-wide HOR distribution                          │
      │     • HOR schematics (homHOR/hetHOR)                        │
      │     • Large duplication zooms                               │
      │     • Family enrichment analysis                            │
      │     Output: plots/ directory                                │
      └─────────────────────────────────────────────────────────────┘

FINAL OUTPUTS:
  • monomer_classifications.tsv  → All monomers with families
  • hors_detected.tsv            → All HORs with coordinates
  • chromosome_stats.tsv         → Per-chromosome statistics
  • plots/                       → Visualizations
```

---

## Read Mode Workflow

```
┌───────────────────────────────────────────────────────────────────────────┐
│                            READ MODE                                       │
│         Long Reads → Satellites → Families → Indels → Variants            │
└───────────────────────────────────────────────────────────────────────────┘

INPUT: reads.fasta (or BAM)
   │
   │  ┌─────────────────────────────────────────────────────────────┐
   └─→│ [1] FasTAN                                                   │
      │     Tandem Repeat Detection in Reads                        │
      │     • Detects tandem arrays per read                        │
      │     • Identifies CEN178-like repeats (~178bp)               │
      │     Output: .1aln format                                    │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  *.1aln
      ┌─────────────────────────────────────────────────────────────┐
      │ [2] tanbed                                                   │
      │     Convert to BED Format                                   │
      │     • Array coordinates per read                            │
      │     • Period and quality scores                             │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  arrays_per_read.bed
      ┌─────────────────────────────────────────────────────────────┐
      │ [3] EXTRACT_MONOMERS                                         │
      │     Extract Individual Monomers Per Read                    │
      │     • Filter by period (160-200bp)                          │
      │     • Extract each repeat unit with coordinates             │
      │     • Track read ID + array ID + monomer ID                 │
      │     Output: read_monomers.fa + monomer_info.tsv             │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  monomers.fa
      ┌─────────────────────────────────────────────────────────────┐
      │ [4] CLASSIFY_MONOMERS                                        │
      │     Per-Monomer Family Assignment                           │
      │     • minimap2: align to reference monomers                 │
      │     • Assign family based on best match                     │
      │     • Track alignment identity                              │
      │     Output: monomer_classifications.tsv ⭐                   │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  classifications.tsv
      ┌─────────────────────────────────────────────────────────────┐
      │ [5] ANALYZE_INDELS (optional)                                │
      │     Indel & Structural Variant Analysis                     │
      │     • Identify insertion/deletion events                    │
      │     • Family-specific indel patterns                        │
      │     • HOR-associated variants                               │
      │     • Deletion-enriched monomers                            │
      │     Output: indel_stats.tsv                                 │
      └─────────────────┬───────────────────────────────────────────┘
                        │
                        ↓  indel_stats.tsv
      ┌─────────────────────────────────────────────────────────────┐
      │ [6] READ_PLOTS                                               │
      │     Read-Level Visualizations                               │
      │     • Per-read family ribbons                               │
      │     • Family transition heatmaps                            │
      │     • Indel pattern plots                                   │
      │     • Array structure diagrams                              │
      │     Output: plots/ directory                                │
      └─────────────────────────────────────────────────────────────┘

FINAL OUTPUTS:
  • monomer_classifications.tsv  → Per-read monomer families
  • indel_stats.tsv              → Indel/variant patterns
  • plots/                       → Per-read visualizations
```

---

## HOR Detection Algorithm (Genome Mode)

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    GAP-AWARE HOR DETECTION                               │
│                  Monomer-Level Pattern Matching                          │
└─────────────────────────────────────────────────────────────────────────┘

INPUT: Family sequence [3, 3, 3, 4, 5, 3, 3, 3, 4, 5, ...]
                             │
                             ↓
        ┌────────────────────────────────────────────────┐
        │ STEP 1: Scan for Repeating Patterns            │
        │   • Try pattern lengths: 3 to 20 monomers      │
        │   • Count consecutive repetitions              │
        │   • Prefer shorter patterns                    │
        └────────────────┬───────────────────────────────┘
                         │
                         ↓
        ┌────────────────────────────────────────────────┐
        │ STEP 2: Check DUAL CRITERIA                    │
        │   ✓ min_copies ≥ 3                             │
        │   ✓ monomers_per_unit ≥ 3                      │
        │                                                 │
        │   Example PASS:                                │
        │     [F3-F3-F3] × 356 = 3F3 ✓                   │
        │                                                 │
        │   Example FAIL:                                │
        │     [F3] × 1068 = 1F3 ✗ (only 1 monomer)      │
        │     [F1-F2-F3-F4-F5] × 2 ✗ (only 2 copies)    │
        └────────────────┬───────────────────────────────┘
                         │
                         ↓
        ┌────────────────────────────────────────────────┐
        │ STEP 3: GAP CHECKING                           │
        │   For each consecutive monomer pair:           │
        │     gap = monomer[i+1].start - monomer[i].end  │
        │                                                 │
        │   If gap > max_gap (500bp):                    │
        │     → BREAK HOR into separate blocks           │
        │                                                 │
        │   Example:                                      │
        │     mono1: [1000-1178]                         │
        │     mono2: [1180-1358]  gap=2bp  ✓ OK         │
        │     mono3: [1358-1536]  gap=0bp  ✓ OK         │
        │     mono4: [2100-2278]  gap=564bp ✗ BREAK!    │
        └────────────────┬───────────────────────────────┘
                         │
                         ↓
        ┌────────────────────────────────────────────────┐
        │ STEP 4: Classify HOR Type                      │
        │                                                 │
        │   homHOR: All same family                      │
        │     Example: 3F3 (F3-F3-F3)                    │
        │                                                 │
        │   hetHOR: Multiple families                    │
        │     Example: 1F1-1F7-1F4                       │
        └────────────────┬───────────────────────────────┘
                         │
                         ↓
OUTPUT: HOR with metadata
  • Pattern: "3F3"
  • Copies: 356
  • Total monomers: 1068
  • Type: homHOR
  • Length: 190.0 kb
  • Start: 2180401
  • End: 2370426
```

---

## Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         DATA FLOW                                        │
└─────────────────────────────────────────────────────────────────────────┘

INPUT FILES                PIPELINE STAGES              OUTPUT FILES
────────────              ──────────────────            ──────────────

genome.fasta    ──→  [FasTAN]          ──→  *.1aln
                                             │
                                             ↓
                                       [tanbed]         ──→  arrays.bed
                     ┌───────────────────────┘
                     │
        ┌────────────┴────────────┐
        ↓                         ↓
    genome.fasta            arrays.bed
        └──────────┬──────────────┘
                   ↓
            [EXTRACT_MONOMERS]    ──→  monomers.fa
                                       monomer_info.tsv
                   ↓
    ┌──────────────┼──────────────────┐
    ↓              ↓                  ↓
monomers.fa   reference.fa    families.txt
    └──────────────┬──────────────────┘
                   ↓
          [CLASSIFY_MONOMERS]   ──→  monomer_classifications.tsv ⭐
                   ↓
    ┌──────────────┴──────────────┐
    ↓                             ↓
classifications.tsv         monomer_info.tsv
    └──────────────┬──────────────┘
                   ↓
            [DETECT_HORS]     ──→  hors_detected.tsv ⭐
                   │                large_duplications.tsv
                   ↓
    ┌──────────────┴──────────────┐
    ↓                             ↓
classifications.tsv          hors.tsv
    └──────────────┬──────────────┘
                   ↓
        [CHROMOSOME_STATS]    ──→  chromosome_stats.tsv
                   │                family_by_chromosome.tsv
                   ↓                hor_by_chromosome.tsv
                   │
                   ↓
          [GENOME_PLOTS]      ──→  plots/*.png
```

---

## Module Dependency Graph

```
                      ┌──────────┐
                      │  INPUT   │
                      └────┬─────┘
                           │
                           ↓
                    ┌──────────────┐
                    │   FASTAN     │
                    └──────┬───────┘
                           │
                           ↓
                    ┌──────────────┐
                    │   TANBED     │
                    └──────┬───────┘
                           │
                           ↓
                ┌──────────────────────┐
                │  EXTRACT_MONOMERS    │
                └──────────┬───────────┘
                           │
                           ↓
                ┌──────────────────────┐
                │  CLASSIFY_MONOMERS   │
                └──────────┬───────────┘
                           │
                           ↓
        ┌──────────────────┴──────────────────┐
        │                                     │
        ↓                                     ↓
┌──────────────┐                    ┌──────────────┐
│ DETECT_HORS  │                    │ANALYZE_INDELS│
│ (genome only)│                    │ (reads only) │
└──────┬───────┘                    └──────┬───────┘
       │                                   │
       ↓                                   ↓
┌──────────────┐                    ┌──────────────┐
│CHROMOSOME_   │                    │  READ_PLOTS  │
│   STATS      │                    │              │
└──────┬───────┘                    └──────────────┘
       │
       ↓
┌──────────────┐
│ GENOME_PLOTS │
└──────────────┘
```

---

## Example: 3-homHOR Detection

```
┌─────────────────────────────────────────────────────────────────────────┐
│                   EXAMPLE: Detecting 3F3 HOR                             │
└─────────────────────────────────────────────────────────────────────────┘

MONOMER SEQUENCE:
┌────┬────┬────┬────┬────┬────┬────┬────┬────┬────┐
│ F3 │ F3 │ F3 │ F3 │ F3 │ F3 │ F3 │ F3 │ F3 │... │  1068 total
└────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘

PATTERN SEARCH:
┌─────────────┐
│ F3│F3│F3   │  ← Pattern: [F3, F3, F3] (length = 3)
└─────────────┘
     ↓
┌─────────────┐
│ F3│F3│F3   │  ← Copy 1
└─────────────┘
┌─────────────┐
│ F3│F3│F3   │  ← Copy 2
└─────────────┘
┌─────────────┐
│ F3│F3│F3   │  ← Copy 3 ... Copy 356
└─────────────┘

CHECK CRITERIA:
✓ Monomers per unit = 3 (≥3 required)
✓ Copies = 356 (≥3 required)
✓ No gaps >500bp

RESULT:
HOR Pattern: 3F3
Copies: 356
Total Monomers: 1068
Type: homHOR (same family)
Length: 190.0 kb
```

---

## Example: hetHOR Detection

```
┌─────────────────────────────────────────────────────────────────────────┐
│                 EXAMPLE: Detecting 1F1-1F7-1F4 HOR                       │
└─────────────────────────────────────────────────────────────────────────┘

MONOMER SEQUENCE:
┌────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┬────┐
│ F1 │ F7 │ F4 │ F1 │ F7 │ F4 │ F1 │ F7 │ F4 │ F1 │ F7 │ F4 │
└────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┴────┘

PATTERN SEARCH:
┌──────────────────┐
│ F1 │ F7 │ F4    │  ← Pattern: [F1, F7, F4] (length = 3)
└──────────────────┘
      ↓
┌──────────────────┐
│ F1 │ F7 │ F4    │  ← Copy 1
└──────────────────┘
┌──────────────────┐
│ F1 │ F7 │ F4    │  ← Copy 2
└──────────────────┘
┌──────────────────┐
│ F1 │ F7 │ F4    │  ← Copy 3
└──────────────────┘
┌──────────────────┐
│ F1 │ F7 │ F4    │  ← Copy 4
└──────────────────┘

CHECK CRITERIA:
✓ Monomers per unit = 3 (≥3 required)
✓ Copies = 4 (≥3 required)
✓ No gaps >500bp

RESULT:
HOR Pattern: 1F1-1F7-1F4
Copies: 4
Total Monomers: 12
Type: hetHOR (multiple families)
```

---

## Chromosome-Aware Statistics

```
┌─────────────────────────────────────────────────────────────────────────┐
│              CHROMOSOME-LEVEL TRACKING                                   │
└─────────────────────────────────────────────────────────────────────────┘

INPUT: Classified monomers + detected HORs
         │
         ↓
    ┌────────────────────────────────┐
    │  Extract chromosome from ID    │
    │  (e.g., "Chr1_..." → Chr1)     │
    └────────┬───────────────────────┘
             │
             ↓
    ┌────────────────────────────────┐
    │ Group by Chromosome            │
    │                                │
    │ Chr1: 151 HORs, 12K monomers   │
    │ Chr2:  31 HORs,  3K monomers   │
    │ Chr3: 101 HORs,  8K monomers   │
    │ Chr4: 122 HORs, 10K monomers   │
    │ Chr5:  36 HORs,  4K monomers   │
    └────────┬───────────────────────┘
             │
             ↓
    ┌────────────────────────────────┐
    │ Calculate per chromosome:      │
    │  • Total monomers              │
    │  • Classification rate         │
    │  • Family distribution         │
    │  • HOR counts (homHOR/hetHOR)  │
    │  • Dominant patterns           │
    └────────┬───────────────────────┘
             │
             ↓
OUTPUT: chromosome_stats.tsv
        family_by_chromosome.tsv
        hor_by_chromosome.tsv

ENABLES:
  → Compare satellite landscapes
  → Identify chromosome-specific families
  → Track HOR prevalence per chromosome
  → Centromere-specific analysis
```

---

## File Format Specifications

### monomer_classifications.tsv

```
monomer_id              seq_id  array_idx  monomer_idx  monomer_start  monomer_end  monomer_family
Chr1_cen_array0_mon0    Chr1    0          0            14853761       14853939     1
Chr1_cen_array0_mon1    Chr1    0          1            14853939       14854117     1
Chr1_cen_array0_mon2    Chr1    0          2            14854117       14854295     1
...
```

### hors_detected.tsv

```
seq_id  hor_start  hor_end    hor_unit   hor_copies  total_monomers  hor_type  hor_length_kb
Chr1    14853761   14854295   3F1        3           9               homHOR    0.534
Chr1    2314       3916       3F1        3           9               homHOR    1.602
Chr4    2180401    2370426    3F3        356         1068            homHOR    190.025
...
```

---

**For complete documentation, see README.md and PIPELINE_SUMMARY.md**
