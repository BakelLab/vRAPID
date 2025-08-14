#!/usr/bin/env python3

import sys

consensus_file, coverage_file, output_file = sys.argv[1:4]

# Read consensus sequence
with open(consensus_file, 'r') as f:
    consensus_seq = ''.join(line.strip() for line in f if not line.startswith('>'))

# Read coverage
coverage = []
with open(coverage_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            coverage.append(int(parts[2]))

# Convert consensus to list for mutability
seq = list(consensus_seq)

# Replace base with 'N' if coverage < 10
for i, cov in enumerate(coverage):
    if i < len(seq) and cov < 10:
        seq[i] = 'n'

# Convert back to string
seq = ''.join(seq)

# Trim sequence if last 10 bases contain >= 3 Ns
while len(seq) >= 10 and seq[-10:].count('n') >= 3:
    seq = seq[:-1].rstrip('n')

# Write output FASTA
sample = consensus_file.split("/")[0]
with open(output_file, 'w') as o:
    o.write(f'>{sample}\n')
    for i in range(0, len(seq), 80):
        o.write(seq[i:i + 80] + '\n')

