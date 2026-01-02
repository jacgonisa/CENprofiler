/*
========================================================================================
    CLASSIFY_MONOMERS: Assign monomers to families
========================================================================================
    Uses minimap2 to map monomers to reference set
    Then assigns families based on phylogenetic classification
========================================================================================
*/

process CLASSIFY_MONOMERS {
    tag "classification"
    publishDir "${params.outdir}/02_monomers", mode: 'copy'

    input:
    path monomers_fasta
    path reference_monomers
    path family_assignments
    path monomer_info

    output:
    path "monomer_classifications.tsv", emit: classifications
    path "classification.log", emit: log
    path "monomers.paf", optional: true, emit: paf

    script:
    """
    # Map monomers to reference using minimap2
    echo "Running minimap2..." > classification.log
    minimap2 \\
        -x map-ont \\
        -t ${task.cpus} \\
        --secondary=no \\
        -c \\
        ${reference_monomers} \\
        ${monomers_fasta} \\
        > monomers.paf \\
        2>> classification.log

    # Parse alignments and assign families
    echo "Assigning families..." >> classification.log
    python3 <<'EOF'
import sys
import pandas as pd
from collections import defaultdict

# Load family assignments
print("Loading family assignments...", file=sys.stderr)
monomer_to_family = {}
with open("${family_assignments}", 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\\t')
        if len(parts) >= 2:
            monomer_id = parts[0]
            family_id = int(parts[1])
            monomer_to_family[monomer_id] = family_id

print(f"Loaded {len(monomer_to_family)} family assignments", file=sys.stderr)

# Load monomer info
monomer_info = pd.read_csv("${monomer_info}", sep='\\t')
print(f"Loaded info for {len(monomer_info)} monomers", file=sys.stderr)

# Parse PAF alignments
print("Parsing alignments...", file=sys.stderr)
alignments = defaultdict(list)

try:
    with open("monomers.paf", 'r') as f:
        for line in f:
            fields = line.strip().split('\\t')
            if len(fields) < 12:
                continue

            query = fields[0]      # monomer_id
            target = fields[5]     # reference monomer
            identity = None

            # Extract identity from alignment
            align_len = int(fields[10])
            block_len = int(fields[11])
            if block_len > 0:
                identity = (align_len / block_len) * 100

            # Get MAPQ
            mapq = int(fields[11]) if len(fields) > 11 else 0

            alignments[query].append({
                'target': target,
                'identity': identity,
                'mapq': mapq
            })
except FileNotFoundError:
    print("WARNING: No PAF file found - all monomers will be unassigned", file=sys.stderr)

print(f"Parsed alignments for {len(alignments)} monomers", file=sys.stderr)

# Assign families based on best alignment
classifications = []

for _, row in monomer_info.iterrows():
    monomer_id = row['monomer_id']

    # Initialize classification
    classification = {
        'monomer_id': monomer_id,
        'best_match': None,
        'alignment_identity': None,
        'mapq': None,
        'monomer_family': None
    }

    # Get best alignment
    if monomer_id in alignments and len(alignments[monomer_id]) > 0:
        # Sort by identity
        best = max(alignments[monomer_id], key=lambda x: x['identity'] if x['identity'] else 0)

        classification['best_match'] = best['target']
        classification['alignment_identity'] = best['identity']
        classification['mapq'] = best['mapq']

        # Assign family if identity passes threshold
        if best['identity'] and best['identity'] >= ${params.min_identity}:
            if best['target'] in monomer_to_family:
                classification['monomer_family'] = monomer_to_family[best['target']]

    classifications.append(classification)

# Merge with monomer info
df_class = pd.DataFrame(classifications)
df_merged = pd.merge(monomer_info, df_class, on='monomer_id', how='left')

# Write output
df_merged.to_csv("monomer_classifications.tsv", sep='\\t', index=False)

# Statistics
total = len(df_merged)
classified = df_merged['monomer_family'].notna().sum()
print(f"Classification complete:", file=sys.stderr)
print(f"  Total monomers: {total}", file=sys.stderr)
print(f"  Classified: {classified} ({classified/total*100:.1f}%)", file=sys.stderr)
print(f"  Unassigned: {total-classified} ({(total-classified)/total*100:.1f}%)", file=sys.stderr)

# Write to log
with open("classification.log", 'a') as f:
    f.write(f"\\nClassification Statistics:\\n")
    f.write(f"Total monomers: {total}\\n")
    f.write(f"Classified: {classified} ({classified/total*100:.1f}%)\\n")
    f.write(f"Unassigned: {total-classified}\\n")

    if classified > 0:
        family_counts = df_merged['monomer_family'].value_counts()
        f.write(f"\\nFamily Distribution:\\n")
        for family, count in family_counts.items():
            f.write(f"  Family {int(family)}: {count}\\n")

EOF

    echo "Classification completed" >> classification.log
    """
}
