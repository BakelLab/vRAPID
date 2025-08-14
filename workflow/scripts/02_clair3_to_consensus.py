#!/usr/bin/env python3
import sys
import gzip

# Input arguments
vcf_file, fasta_file, chrs, coverage_file, output_file = sys.argv[1:6]

# -----------------------------
# Step 1: Parse VCF for ins/dels
# -----------------------------
ins = set()
dels = set()

with gzip.open(vcf_file, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue  # Skip headers
        chrom, pos_str, vid, ref, alt_str, qual, filt, info = line.strip().split('\t', 7)
        pos = int(pos_str)
        alts = alt_str.split(',')
        for alt in alts:
            if alt == '.':
                continue
            if len(alt) > len(ref):
                ins.add(pos - 1)  # 0-based coordinate
            elif len(ref) > len(alt):
                dels.update(range(pos - 1, pos - 1 + len(ref) - 1))

# -----------------------------
# Step 2: Read coverage file
# -----------------------------
with open(coverage_file, 'r') as f:
    coverage = [int(line.strip().split('\t')[2]) for line in f if line.strip()]

# -----------------------------
# Step 3: Read FASTA sequence
# -----------------------------
with open(fasta_file, 'r') as f:
    seq = ''.join(line.strip() for line in f if not line.startswith('>'))

seq = list(seq)  # make mutable

# -----------------------------
# Step 4: Apply ins/dels and coverage
# -----------------------------
qnum = 0
for refnum, cov in enumerate(coverage):
    while qnum in ins:
        qnum += 1
    if refnum in dels:
        continue
    if cov < 10:
        seq[qnum] = 'n'
    qnum += 1

# -----------------------------
# Step 5: Clean and trim sequence
# -----------------------------
seq = ''.join(seq)
seq = seq.rstrip('a').rstrip('n')
seq = seq.lstrip('a').lstrip('n')

while len(seq) >= 10 and seq[-10:].count('n') >= 3:
    seq = seq[:-1].rstrip('n')

# -----------------------------
# Step 6: Write output FASTA
# -----------------------------
sample = fasta_file.split("/")[0]
with open(output_file, 'w') as o:
    o.write(f'>{sample}\n')
    for i in range(0, len(seq), 80):
        o.write(seq[i:i + 80] + '\n')

