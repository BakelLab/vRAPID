#!/usr/bin/env python3

import sys

changes_file, ins_file, dels_file = sys.argv[1:4]

ins = set()
dels = set()

with open(changes_file, 'r') as f:
    lines = list(f)

i = 0
while i < len(lines):
    line = lines[i].strip()
    parts = line.split()

    # Handle insertions (two-line format)
    if len(parts) == 2 and (i + 1) < len(lines) and lines[i + 1].startswith(" "):
        # This is an insertion line
        contig, pos_str = parts[0].split(':')
        pos = int(pos_str)
        ins.add(pos - 1)  # Convert to 0-based
        i += 2  # Skip next line (inserted sequence)
        continue

    # Handle deletions
    if len(parts) >= 4 and parts[3] == '.':
        if '-' in parts[0]:
            start, stop = map(int, parts[0].split(':')[1].split('-'))
        else:
            start = stop = int(parts[0].split(':')[1])
        dels.update(range(start - 1, stop))
    
    # Handle insertions (old '.' base)
    if len(parts) >= 4 and parts[2] == '.':
        if '-' in parts[1]:
            start, stop = map(int, parts[1].split(':')[1].split('-'))
        else:
            start = stop = int(parts[1].split(':')[1])
        ins.update(range(start - 1, stop))

    i += 1

with open(ins_file, 'w') as f:
    for pos in sorted(ins):
        f.write(f"{pos}\n")

with open(dels_file, 'w') as f:
    for pos in sorted(dels):
        f.write(f"{pos}\n")

