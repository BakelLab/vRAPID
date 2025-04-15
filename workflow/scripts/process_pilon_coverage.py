#!/usr/bin/env python3

import sys

fasta_file, chrs, ins_file, dels_file, coverage_file, output_file = sys.argv[1:7]

# Read the sequence file
with open(fasta_file, 'r') as f:
    seq = ''.join(line.rstrip() for line in f if not line.startswith('>'))

# Read the insertion and deletion files into sets for faster lookups
with open(ins_file, 'r') as f:
    ins = set(int(line.strip()) for line in f)

with open(dels_file, 'r') as f:
    dels = set(int(line.strip()) for line in f)

# Read the coverage file and modify sequence based on coverage
with open(coverage_file, 'r') as f:
    f.readline()  # Skip header lines
    f.readline()
    coverage = [int(line.rstrip()) for line in f if line.strip()]


# Convert seq to a list for mutability
seq = list(seq)

# Modify seq based on coverage and ins/dels lists
qnum = 0
for refnum, coverage_value in enumerate(coverage):
    while qnum in ins:
        qnum += 1
    if refnum in dels:
        continue
    if coverage_value < 10:
        seq[qnum] = 'n'
    qnum += 1

# Convert seq back to a string
seq = ''.join(seq)

# Clean up the sequence by removing trailing 'a' and 'n'
seq = seq.rstrip('a').rstrip('n')
seq = seq.lstrip('a').lstrip('n')

# Trim sequence based on specific conditions
while len(seq) >= 10 and seq[-10:].count('n') >= 3:
    seq = seq[:-1].rstrip('n')

sample = fasta_file.split("/")[0]
with open(output_file, 'w') as o:
            o.write(f'>{sample}|{chrs}\n')
            for i in range(0, len(seq), 80):
                o.write(seq[i:i + 80] + '\n')
