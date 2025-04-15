#!/usr/bin/env python3

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pileup", required=True)
parser.add_argument("--changes", required=True)
parser.add_argument("--primers", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--min_ratio", type=float, default=0.8)
args = parser.parse_args()

# Read changes
changes = {}
with open(args.changes) as f:
    for line in f:
        pos, ref, alt = line.strip().split()[:3]
        changes[pos] = (ref, alt)

# Read primer positions
primer_positions = set()
with open(args.primers) as f:
    for line in f:
        primer_positions.add(int(line.strip()))

outfile = os.path.join(args.outdir, f"{args.sample}.{args.chrom}_variable_bases.tsv")

with open(args.pileup) as f, open(outfile, "w") as out:
    out.write("reference\tposition\tflagged\tin_primer\treference_base\tpilon_base\tdepth\treference_base_fraction\tforward_depth"
              "\tforward_fraction\treverse_depth\treverse_fraction\tA\tT\tC\tG\ta\tt\tc\tg\tn\tinsertion\tdeletion\n")
    for line in f:
        ref, pos, refbase, cov, seq, qual = line.strip().split()
        refbase = refbase.lower()
        outseq = ""

        if pos in changes and len(changes[pos][0]) == 1 and len(changes[pos][1]) == 1 and changes[pos][0] != '.':
            if refbase != changes[pos][0].lower():
                sys.exit(f"Pileup doesn't match pilon output at position {pos}")
            new_refbase = changes[pos][1].lower()
        elif pos in changes:
            new_refbase = changes[pos][1].lower()
        else:
            new_refbase = refbase

        counts = {b: 0 for b in ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G', 'n', 'I', 'D']}
        depth = forward_depth = reverse_depth = 0
        seq = list(seq)
        getdel = getins = False

        while seq:
            x = seq.pop(0)
            mod = None
            if x == '.':
                mod = refbase.upper()
                forward_depth += 1
            elif x == ',':
                mod = refbase.lower()
                reverse_depth += 1
            elif x == '+':
                getins = True
                digit = ''
            elif x == '-':
                getdel = True
                digit = ''
            elif x.isdigit() and (getdel or getins):
                digit += x
            elif getdel:
                for _ in range(int(digit) - 1):
                    seq.pop(0)
                mod = 'D'
                getdel = False
            elif getins:
                for _ in range(int(digit) - 1):
                    seq.pop(0)
                mod = 'I'
                getins = False
            elif x == '^':
                seq.pop(0)
            elif x in counts:
                mod = x
                if x.islower():
                    reverse_depth += 1
                else:
                    forward_depth += 1
            elif x in ['$', '*']:
                pass
            else:
                sys.exit(f"Base not recognised: {x}")

            if mod:
                counts[mod] += 1
                depth += 1

        try:
            fraction = (counts[new_refbase.upper()] + counts[new_refbase.lower()]) / depth
        except ZeroDivisionError:
            fraction = None

        try:
            forward_fraction = counts[new_refbase.upper()] / forward_depth
        except ZeroDivisionError:
            forward_fraction = "no_cov"

        try:
            reverse_fraction = counts[new_refbase.lower()] / reverse_depth
        except ZeroDivisionError:
            reverse_fraction = "no_cov"

        if fraction is None:
            flagged = "no_cov"
            outseq += 'n'
        elif int(cov) < 10:
            flagged = "low_cov"
            outseq += 'n'
        elif pos in changes and (len(changes[pos][0]) != 1 or len(changes[pos][1]) != 1):
            flagged = "indel"
        elif (forward_fraction == "no_cov" or forward_fraction < args.min_ratio) and \
             (reverse_fraction == "no_cov" or reverse_fraction < args.min_ratio):
            flagged = "FLAGGED"
            outseq += 'n'
        else:
            flagged = "OK"
            outseq += new_refbase

        inprimer = "Y" if int(pos) in primer_positions else "N"

        outlist = [ref, pos, flagged, inprimer, refbase, new_refbase, cov, fraction, forward_depth,
                   forward_fraction, reverse_depth, reverse_fraction]
        outlist.extend([counts[b] for b in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'n', 'I', 'D']])
        out.write("\t".join(map(str, outlist)) + "\n")

