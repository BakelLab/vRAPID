#!/usr/bin/env python3

import os
import sys
import argparse
import logging

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze pileup, apply Pilon changes, and flag variable bases.")
    parser.add_argument("--sample_folder", required=True, help="Path to the sample folder containing Pilon changes")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--chrs", required=True, help="Chromosome or reference ID")
    parser.add_argument("--primer_set", required=True, help="Path to the combined primer set file")
    parser.add_argument("--min_ratio", type=float, default=0.9, help="Minimum allele fraction to avoid flagging")
    return parser.parse_args()

def parse_primers(filepath):
    primer_positions = set()
    with open(filepath) as f:
        f.readline()  # skip header
        for line in f:
            try:
                name, seq, pool, length, tm, gc, start, end = line.split()
                for i in range(min(int(start), int(end)), max(int(start), int(end))):
                    primer_positions.add(i)
            except ValueError:
                continue
    return primer_positions

def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    sample_folder, outdir, chrs = args.sample_folder, args.outdir, args.chrs

    # Read Pilon changes
    changes = {}
    pilon_changes_path = f"{sample_folder}/02_assembly/{chrs}_pilon.changes"
    with open(pilon_changes_path) as f:
        for line in f:
            pos = line.split()[0].split(':')[1]
            ref_base = line.split()[2]
            new_base = line.split()[3]
            changes[pos] = (ref_base, new_base)

    primer_positions = parse_primers(args.primer_set)

    modlist = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'n', 'I', 'D']
    outseq = ""

    pileup_path = f"{outdir}/{chrs}_pileup"
    output_table = f"{outdir}/{chrs}_variable_bases.tsv"

    with open(pileup_path) as f, open(output_table, "w") as out:
        out.write("reference\tposition\tflagged\tin_primer\treference_base\tpilon_base\tdepth\treference_base_fraction\tforward_depth"
                "\tforward_fraction\treverse_detph\treverse_fraction\tA\tT\tC\tG\ta\tt\tc\tg\tn\tinsertion\tdeletion\n") 
        for line in f:
            ref, pos, refbase, cov, seq, qual = line.split()
            refbase = refbase.lower()
            pos_int = int(pos)

            if pos in changes and len(changes[pos][0]) == 1 and len(changes[pos][1]) == 1 and changes[pos][0] != '.':
                if refbase != changes[pos][0].lower():
                    sys.exit("pileup doesn't match Pilon output at position %s" % pos)
                new_refbase = changes[pos][1].lower()
            elif pos in changes:
                new_refbase = changes[pos][1].lower()
            else:
                new_refbase = refbase

            counts = {x: 0 for x in modlist}
            depth = 0
            forward_depth = 0
            reverse_depth = 0
            seq = list(seq)
            getdel = getins = False
            digit = ""

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
                elif x.isdigit() and (getins or getdel):
                    digit += x
                elif getins:
                    for _ in range(int(digit) - 1):
                        seq.pop(0)
                    mod = 'I'
                    getins = False
                elif getdel:
                    for _ in range(int(digit) - 1):
                        seq.pop(0)
                    mod = 'D'
                    getdel = False
                elif x == '^':
                    seq.pop(0)  # skip mapping quality
                elif x in ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']:
                    mod = x
                    if x.islower():
                        reverse_depth += 1
                    else:
                        forward_depth += 1
                elif x in ['$', '*']:
                    continue
                else:
                    sys.exit("Unrecognized base in pileup: " + x)

                if mod:
                    counts[mod] += 1
                    depth += 1

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

            in_primer = "Y" if pos_int in primer_positions else "N"

            outlist = [ref, pos, flagged, in_primer, refbase, new_refbase, cov, fraction,
                       forward_depth, forward_fraction, reverse_depth, reverse_fraction]
            outlist += [counts[i] for i in modlist]
            out.write('\t'.join(map(str, outlist)) + '\n')

    with open(f"{outdir}/{chrs}.final.fna", "w") as out:
        out.write(">%s\n" % chrs)
        for i in range(0, len(outseq), 80):
            out.write(outseq[i:i+80] + '\n')

    logging.info("Finished processing %s", chrs)

if __name__ == "__main__":
    main()

