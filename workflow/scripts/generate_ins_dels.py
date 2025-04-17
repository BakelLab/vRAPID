#!/usr/bin/env python3

import sys

changes_file, ins_file, dels_file = sys.argv[1:4]

ins = set()
dels = set()

with open(changes_file, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if parts[2] == '.':
            if '-' in parts[1]:
                start, stop = map(int, parts[1].split(':')[1].split('-'))
            else:
                start = stop = int(parts[1].split(':')[1])
            ins.update(range(start - 1, stop))
        if parts[3] == '.':
            if '-' in parts[0]:
                start, stop = map(int, parts[0].split(':')[1].split('-'))
            else:
                start = stop = int(parts[0].split(':')[1])
            dels.update(range(start - 1, stop))

with open(ins_file, 'w') as f:
    for pos in sorted(ins):
        f.write(f"{pos}\n")

with open(dels_file, 'w') as f:
    for pos in sorted(dels):
        f.write(f"{pos}\n")

