#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import pandas as pd
import pysam


def read_bed_3col(filepath):
    """
    Reads a 3-column BED file (chrom, start, end) into a pandas DataFrame.
    """
    return pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": int, "end": int},
        usecols=[0, 1, 2]
    )


def mask_reference(reference_fasta, bed_file, vcf_gz_file, output_fasta):
    """
    Masks positions in a reference FASTA using:
      1. BED intervals (mask with 'N')
      2. VCF variants from compressed file (mask REF-length bases)
    """
    # Load reference into dict of mutable sequences
    seqs = {rec.id: list(rec.seq) for rec in SeqIO.parse(reference_fasta, "fasta")}

    # Mask from BED
    bed_df = read_bed_3col(bed_file)
    for _, region in bed_df.iterrows():
        chrom, start, end = region["chrom"], region["start"], region["end"]
        if chrom not in seqs:
            continue
        for pos in range(start, end):  # BED is 0-based
            seqs[chrom][pos] = "N"

    # Mask from compressed VCF
    vcf_in = pysam.VariantFile(vcf_gz_file, "r")
    for rec in vcf_in.fetch():
        if rec.chrom not in seqs:
            continue
        for offset in range(len(rec.ref)):  # REF length bases
            seqs[rec.chrom][rec.pos - 1 + offset] = "N"  # VCF is 1-based

    # Write masked FASTA
    with open(output_fasta, "w") as out_fh:
        for chrom, seq_list in seqs.items():
            out_fh.write(f">{chrom}\n{''.join(seq_list)}\n")


def main():
    parser = argparse.ArgumentParser(description="Mask a reference sequence using BED and compressed VCF files.")
    parser.add_argument("reference", help="Reference FASTA file.")
    parser.add_argument("maskfile", help="BED file with mask intervals.")
    parser.add_argument("maskvcf", help="Compressed VCF file (.vcf.gz) with variants to mask.")
    parser.add_argument("output", help="Output masked FASTA file.")
    args = parser.parse_args()

    mask_reference(args.reference, args.maskfile, args.maskvcf, args.output)


if __name__ == "__main__":
    main()

