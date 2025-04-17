import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pileup", required=True)
parser.add_argument("--changes", required=True)
parser.add_argument("--primers", required=True)
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

outfile = os.path.join(args.outdir, f"{args.chrom}_variable_bases.tsv")
outseq = ""  # Moved outside to accumulate full sequence

with open(args.pileup) as f, open(outfile, "w") as out:
    out.write("reference\tposition\tflagged\tin_primer\treference_base\tpilon_base\tdepth\treference_base_fraction\tforward_depth"
              "\tforward_fraction\treverse_depth\treverse_fraction\tA\tT\tC\tG\ta\tt\tc\tg\tn\tinsertion\tdeletion\n")
    
    for line in f:
        ref, pos, refbase, cov, seq, qual = line.split()
        refbase = refbase.lower()
        
        # Modify base according to changes (if present)
        if pos in changes and len(changes[pos][0]) == 1 and len(changes[pos][1]) == 1 and changes[pos][0] != '.':
            if refbase != changes[pos][0].lower():
                sys.exit("pileup doesn't match pilon output")
            new_refbase = changes[pos][1].lower()
        elif pos in changes:
            new_refbase = changes[pos][1].lower()
        else:
            new_refbase = refbase

        counts = {'a': 0, 't': 0, 'c': 0, 'g': 0, 'A': 0, 'T': 0, 'C': 0, 'G': 0, 'n': 0, 'I': 0, 'D': 0}
        depth = 0
        seq = list(seq)
        getdel = False
        getins = False
        forward_depth = 0
        reverse_depth = 0

        # Process sequence
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
            elif x in ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']:
                mod = x
                if x.islower():
                    reverse_depth += 1
                else:
                    forward_depth += 1
            elif x in ['$', '*']:
                pass
            else:
                sys.exit(f"Base not recognised {x}")

            if mod:
                counts[mod] += 1
                depth += 1

        # Calculate fractions
        try:
            fraction = (counts[new_refbase.upper()] + counts[new_refbase.lower()]) / depth
        except KeyError:
            continue
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

        # Determine flagged status
        if fraction is None:
            flagged = "no_cov"
            outseq += 'n'
        elif int(cov) < 10:
            flagged = "low_cov"
            outseq += 'n'
        elif pos in changes and (len(changes[pos][0]) != 1 or len(changes[pos][1]) != 1):
            flagged = 'indel'
        elif (forward_fraction == "no_cov" or forward_fraction < args.min_ratio) and \
             (reverse_fraction == "no_cov" or reverse_fraction < args.min_ratio):
            flagged = "FLAGGED"
            outseq += 'n'
        else:
            flagged = "OK"
            outseq += new_refbase

        # Check if position is in primer regions
        inprimer = "Y" if int(pos) in primer_positions else "N"

        # Write results
        outlist = [ref, pos, flagged, inprimer, refbase, new_refbase, cov, fraction, forward_depth,
                   forward_fraction, reverse_depth, reverse_fraction]
        for i in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'n', 'I', 'D']:
            outlist.append(counts[i])
        out.write('\t'.join(map(str, outlist)) + '\n')

# Write final sequence
with open(f"{args.outdir}/{args.chrom}.final.fna", "w") as out:
    out.write(f">{args.chrom}\n")
    for i in range(0, len(outseq), 80):
        out.write(outseq[i:i+80] + '\n')

