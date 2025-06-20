#!/usr/bin/env python3

import subprocess
import os

# Define file paths from Snakemake
fasta_file = snakemake.input.fasta
reads_1 = snakemake.input.reads_1
reads_2 = snakemake.input.reads_2
output_dir = snakemake.output.folder 
directory = os.path.dirname(fasta_file)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Step 1: Extract sequence lengths using seqkit
chrom_lengths_file = os.path.join(directory, "chromosome_lengths.txt")
subprocess.run(
    f"seqkit fx2tab --length --name {fasta_file} > {chrom_lengths_file}",
    shell=True, check=True
)

# Step 2: Loop through each chromosome and run shovill
with open(chrom_lengths_file, "r") as file:
    for line in file:
        chrs_full, length = line.strip().split()
        chrs = chrs_full.split("_")[1]  # Assumes headers like >chr_1, >chr_2, etc.

        print(f"Running shovill for chromosome {chrs} with genome size {length}")

        shovill_command = [
            "shovill",
            "--outdir", output_dir,
            "--R1", reads_1,
            "--R2", reads_2,
            "--gsize", length,
            "--cpus", "1",
            "--force"
        ]

        try:
            subprocess.run(shovill_command, check=True)
        except subprocess.CalledProcessError:
            print(f"Error: shovill failed for chromosome {chrs}")
            raise

# Step 3: Final touch to ensure Snakemake sees the output directory
# Even if shovill does not output anything, ensure directory exists
os.makedirs(output_dir, exist_ok=True)


# Clean up
#subprocess.run(f"rm {directory}/chromosome_lengths.txt", shell=True, check=True)

