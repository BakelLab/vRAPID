#!/usr/bin/env python3

import subprocess
import sys
import os

# Define file paths
fasta_file = snakemake.input.fasta
reads_1 = snakemake.input.reads_1
reads_2 = snakemake.input.reads_2
output_dir = snakemake.output.folder
directory = os.path.dirname(fasta_file)

# Step 1: Extract sequence lengths using seqkit
subprocess.run(f"./seqkit fx2tab --length --name {fasta_file} > {directory}/chromosome_lengths.txt", shell=True, check=True)

# Step 2: Read the chromosome lengths and run shovill for each chromosome
with open(f"{directory}/chromosome_lengths.txt", "r") as file:
    for line in file:
        chrs, length = line.split()  # Split the line into chromosome and length
        chrs = chrs.split("_")[1]
        print(f"Running shovill for chromosome {chrs} with genome size {length}")

        # Run shovill command
 #       shovill_outdir = f"{chr}_shovill/"
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
        except subprocess.CalledProcessError as e:
            print(f"Error: shovill command failed for chromosome {chrs}")
            raise e

# Clean up
#subprocess.run(f"rm {directory}/chromosome_lengths.txt", shell=True, check=True)
