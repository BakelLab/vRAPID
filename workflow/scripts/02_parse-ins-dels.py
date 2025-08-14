#!/usr/bin/env python3

import sys
import gzip

vcf_file, ins_file, dels_file = sys.argv[1:4]

ins = set()
dels = set()

# Open gzipped VCF
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

            # Insertion: ALT is longer than REF
            if len(alt) > len(ref):
                ins.add(pos - 1)  # 0-based coordinate

            # Deletion: REF is longer than ALT
            elif len(ref) > len(alt):
                dels.update(range(pos - 1, pos - 1 + len(ref) - 1))

# Write insertion positions
with open(ins_file, 'w') as f:
    for p in sorted(ins):
        f.write(f"{p}\n")

# Write deletion positions
with open(dels_file, 'w') as f:
    for p in sorted(dels):
        f.write(f"{p}\n")

