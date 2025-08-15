#!/usr/bin/env python3
import argparse
import os
import sys
import pysam
from Bio import SeqIO
from itertools import groupby


def get_intervals(positions):
    """Turn a list of positions into a list of continuous intervals."""
    if not positions:
        return []
    positions = sorted(set(positions))
    intervals = []
    for _, group in groupby(enumerate(positions), lambda x: x[1] - x[0]):
        group = list(group)
        intervals.append([group[0][1], group[-1][1]])
    return intervals


def collect_depths(bam_path, ref_name, min_depth, ignore_deletions=False, warn_rg_coverage=False):
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Ensure reference exists in BAM
    if ref_name not in bam.references:
        bam.close()
        raise ValueError(f"Reference '{ref_name}' not found in BAM")

    ref_len = bam.get_reference_length(ref_name)
    total_depths = [0] * ref_len

    # Initialize readgroup-specific depths
    rg_depths = {}
    for rg in bam.header.get("RG", []):
        if rg.get("ID") != "unmatched":
            rg_depths[rg["ID"]] = [0] * ref_len

    low_rg_positions = set()

    for col in bam.pileup(ref_name, max_depth=10000, truncate=False, min_base_quality=0):
        for pileupread in col.pileups:
            if pileupread.is_refskip:
                continue

            rg = pileupread.alignment.get_tag("RG")
            if rg not in rg_depths:
                raise ValueError(f"Unexpected readgroup in BAM: {rg}")

            if pileupread.is_del:
                if not ignore_deletions:
                    total_depths[col.pos] += 1
                    rg_depths[rg][col.pos] += 1
            else:
                total_depths[col.pos] += 1
                rg_depths[rg][col.pos] += 1

        # Mask if below threshold
        if total_depths[col.pos] < min_depth:
            total_depths[col.pos] = 0
        else:
            if all(rg_depths[rg][col.pos] < min_depth for rg in rg_depths):
                total_depths[col.pos] = 0
                low_rg_positions.add(col.pos)

    bam.close()

    if warn_rg_coverage and low_rg_positions:
        sys.stderr.write(
            f"Warning: {bam_path} contains positions with low per-readgroup coverage "
            f"but sufficient combined coverage (min depth={min_depth}).\n"
        )
        for start, end in get_intervals(low_rg_positions):
            sys.stderr.write(f"Low RG coverage region: {start}-{end}\n")

    return total_depths, rg_depths


def main():
    parser = argparse.ArgumentParser(description="Generate a depth mask from a BAM file.")
    parser.add_argument("--depth", type=int, default=20, help="Minimum depth to keep a position unmasked.")
    parser.add_argument("--warn-rg-coverage", action="store_true", help="Warn about low per-readgroup coverage.")
    parser.add_argument("--ignore-deletions", action="store_true", help="Ignore deletions when counting depth.")
    parser.add_argument("--store-rg-depths", action="store_true", help="Write per-readgroup depths to files.")
    parser.add_argument("reference", help="Reference FASTA file.")
    parser.add_argument("bamfile", help="Input BAM file.")
    parser.add_argument("outfile", help="Output mask file.")
    args = parser.parse_args()

    record = next(SeqIO.parse(args.reference, "fasta"))
    ref_name = record.id
    ref_length = len(record.seq)

    depths, rg_depths = collect_depths(
        args.bamfile, ref_name, args.depth, args.ignore_deletions, args.warn_rg_coverage
    )

    if len(depths) != ref_length:
        sys.stderr.write("Warning: depth vector length does not match reference length.\n")

    if args.store_rg_depths:
        for rg, dlist in rg_depths.items():
            with open(f"{args.outfile}.{rg}.depths", "w") as f:
                for pos, depth in enumerate(dlist):
                    f.write(f"{ref_name}\t{rg}\t{pos}\t{depth}\n")

    mask_positions = [pos for pos, d in enumerate(depths) if d == 0]
    with open(args.outfile, "w") as f:
        for start, end in get_intervals(mask_positions):
            f.write(f"{ref_name}\t{start+1}\t{end+1}\n")


if __name__ == "__main__":
    main()

