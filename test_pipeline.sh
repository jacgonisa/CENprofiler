#!/bin/bash
#
# CENprofiler Test Suite
# Tests both genome and read modes with small datasets
#

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================"
echo " CENprofiler Test Suite"
echo "========================================"
echo ""

# Test data paths
TEST_DIR="test_data"
REFERENCE_MONOMERS="/mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/Col-CC-V2-CEN178-representative.fasta"
FAMILY_ASSIGNMENTS="/mnt/ssd-8tb/atrx_china/kmers_and_other_classification_methods/results_phylo_subsampling/itol_manual_phylo_clusters.txt"

# Check test data exists
if [ ! -d "$TEST_DIR" ]; then
    mkdir -p "$TEST_DIR"
    echo -e "${YELLOW}Creating test data directory${NC}"
fi

# Check reference files exist
if [ ! -f "$REFERENCE_MONOMERS" ]; then
    echo -e "${RED}ERROR: Reference monomers not found at $REFERENCE_MONOMERS${NC}"
    exit 1
fi

if [ ! -f "$FAMILY_ASSIGNMENTS" ]; then
    echo -e "${RED}ERROR: Family assignments not found at $FAMILY_ASSIGNMENTS${NC}"
    exit 1
fi

echo -e "${GREEN}✓${NC} Reference files found"
echo ""

# ===========================================
# TEST 1: Genome Mode (using existing data)
# ===========================================

echo "Test 1: GENOME MODE"
echo "-------------------"

# Use existing test genome data
GENOME_INPUT="/mnt/ssd-8tb/atrx_china/CEN178profiler/v3_genome_analysis_CORRECT/centromeric_regions.fa"

if [ -f "$GENOME_INPUT" ]; then
    echo "Using test genome: $GENOME_INPUT"

    OUTPUT_GENOME="test_data/results_genome"
    rm -rf "$OUTPUT_GENOME"

    echo "Running genome mode pipeline..."
    nextflow run main.nf \
        --mode genome \
        --input "$GENOME_INPUT" \
        --reference_monomers "$REFERENCE_MONOMERS" \
        --family_assignments "$FAMILY_ASSIGNMENTS" \
        --outdir "$OUTPUT_GENOME" \
        --fastan_threads 4 \
        --minimap2_threads 2 \
        -resume

    # Check outputs
    if [ -f "$OUTPUT_GENOME/02_monomers/monomer_classifications.tsv" ]; then
        n_monomers=$(tail -n +2 "$OUTPUT_GENOME/02_monomers/monomer_classifications.tsv" | wc -l)
        echo -e "${GREEN}✓${NC} Genome mode completed: $n_monomers monomers classified"

        if [ -f "$OUTPUT_GENOME/03_hors/hors_detected.tsv" ]; then
            n_hors=$(tail -n +2 "$OUTPUT_GENOME/03_hors/hors_detected.tsv" | wc -l)
            echo -e "${GREEN}✓${NC} HOR detection completed: $n_hors HORs detected"
        fi
    else
        echo -e "${RED}✗${NC} Genome mode failed: Output file not found"
        exit 1
    fi
else
    echo -e "${YELLOW}⊘${NC} Skipping genome test - test genome not found"
fi

echo ""

# ===========================================
# TEST 2: Read Mode (using existing data)
# ===========================================

echo "Test 2: READ MODE"
echo "-----------------"

# Use existing test reads
READ_INPUT="/mnt/ssd-8tb/atrx_china/CEN178profiler/results_v2_test/test_sv_reads.fa"

if [ -f "$READ_INPUT" ]; then
    echo "Using test reads: $READ_INPUT"

    OUTPUT_READS="test_data/results_reads"
    rm -rf "$OUTPUT_READS"

    echo "Running read mode pipeline..."
    nextflow run main.nf \
        --mode reads \
        --input "$READ_INPUT" \
        --reference_monomers "$REFERENCE_MONOMERS" \
        --family_assignments "$FAMILY_ASSIGNMENTS" \
        --outdir "$OUTPUT_READS" \
        --fastan_threads 4 \
        --minimap2_threads 2 \
        --analyze_indels true \
        -resume

    # Check outputs
    if [ -f "$OUTPUT_READS/02_monomers/monomer_classifications.tsv" ]; then
        n_monomers=$(tail -n +2 "$OUTPUT_READS/02_monomers/monomer_classifications.tsv" | wc -l)
        n_classified=$(tail -n +2 "$OUTPUT_READS/02_monomers/monomer_classifications.tsv" | \
                      awk -F'\t' '$NF != "" && $NF != "NA"' | wc -l)

        echo -e "${GREEN}✓${NC} Read mode completed: $n_monomers monomers, $n_classified classified"

        if [ -f "$OUTPUT_READS/04_indels/indel_stats.tsv" ]; then
            echo -e "${GREEN}✓${NC} Indel analysis completed"
        fi
    else
        echo -e "${RED}✗${NC} Read mode failed: Output file not found"
        exit 1
    fi
else
    echo -e "${YELLOW}⊘${NC} Skipping read test - test reads not found"
fi

echo ""

# ===========================================
# Summary
# ===========================================

echo "========================================"
echo " Test Summary"
echo "========================================"
echo ""

if [ -d "$OUTPUT_GENOME" ]; then
    echo "Genome Mode Results:"
    echo "  Location: $OUTPUT_GENOME"

    if [ -f "$OUTPUT_GENOME/02_monomers/monomer_classifications.tsv" ]; then
        total=$(tail -n +2 "$OUTPUT_GENOME/02_monomers/monomer_classifications.tsv" | wc -l)
        classified=$(tail -n +2 "$OUTPUT_GENOME/02_monomers/monomer_classifications.tsv" | \
                    awk -F'\t' '$NF != "" && $NF != "NA"' | wc -l)
        pct=$(awk "BEGIN {printf \"%.1f\", ($classified/$total)*100}")
        echo "  Monomers: $total total, $classified classified ($pct%)"
    fi

    if [ -f "$OUTPUT_GENOME/03_hors/hors_detected.tsv" ]; then
        n_hors=$(tail -n +2 "$OUTPUT_GENOME/03_hors/hors_detected.tsv" | wc -l)
        echo "  HORs: $n_hors detected"
    fi
    echo ""
fi

if [ -d "$OUTPUT_READS" ]; then
    echo "Read Mode Results:"
    echo "  Location: $OUTPUT_READS"

    if [ -f "$OUTPUT_READS/02_monomers/monomer_classifications.tsv" ]; then
        total=$(tail -n +2 "$OUTPUT_READS/02_monomers/monomer_classifications.tsv" | wc -l)
        classified=$(tail -n +2 "$OUTPUT_READS/02_monomers/monomer_classifications.tsv" | \
                    awk -F'\t' '$NF != "" && $NF != "NA"' | wc -l)
        pct=$(awk "BEGIN {printf \"%.1f\", ($classified/$total)*100}")
        echo "  Monomers: $total total, $classified classified ($pct%)"
    fi
    echo ""
fi

echo -e "${GREEN}All tests completed successfully!${NC}"
echo ""
echo "Review results in test_data/ directory"
echo "Pipeline logs: test_data/results_*/pipeline_info/"
echo ""
