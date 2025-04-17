#!/usr/bin/env python3

import re

def process_pilon_changes(pilon_files, output_files):
    # Ensure inputs are lists
    if isinstance(pilon_files, str):
        pilon_files = [pilon_files]
    if isinstance(output_files, str):
        output_files = [output_files]

    for pilon_file, output_file in zip(pilon_files, output_files):
        changes = {}

        with open(pilon_file) as f:
            for line in f:
                parts = line.strip().split()
                pos = parts[0].split(':')[1]
                ref_base = parts[2]
                new_base = parts[3]
                changes[pos] = (ref_base, new_base)

        with open(output_file, "w") as out:
            for pos, (ref_base, new_base) in changes.items():
                out.write(f"{pos}\t{ref_base}\t{new_base}\n")


def process_primer_file(primer_file, output_file):
    primer_positions = set()

    with open(primer_file) as f:
        f.readline()  # Skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) < 8:
                continue  # Skip malformed lines
            name, seq, pool, length, tm, gc, start, end = parts

            match = re.search(r'_(\d+\.?\d*)kb_', name)
            if not match:
                raise ValueError(f"Could not parse primer set from name: {name}")

            start, end = int(start), int(end)
            for i in range(min(start, end), max(start, end)):
                primer_positions.add(i)

    with open(output_file, "w") as out:
        for pos in sorted(primer_positions):
            out.write(f"{pos}\n")


def main():
    process_pilon_changes(snakemake.input.pilon_changes, snakemake.output.changes)
    process_primer_file(snakemake.input.primer_file, snakemake.output.primer_positions)


if __name__ == "__main__":
    main()

