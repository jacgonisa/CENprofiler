#!/usr/bin/env python3
"""
Prepare monomer classifications for visualization by adding read_id column
"""
import pandas as pd
import sys

if len(sys.argv) != 3:
    print("Usage: prepare_monomer_data_for_viz.py <input_tsv> <output_tsv>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load monomers
df = pd.read_csv(input_file, sep='\t')

# Add read_id column from seq_id if not present
if 'read_id' not in df.columns and 'seq_id' in df.columns:
    df['read_id'] = df['seq_id']

# Save
df.to_csv(output_file, sep='\t', index=False)
print(f"Prepared {len(df)} monomers with read_id column")
print(f"Saved to: {output_file}")
