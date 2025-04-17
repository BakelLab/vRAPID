#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
from PyPDF2 import PdfReader, PdfWriter
import os

def add_variant_plot_to_pdf(pileup_file, input_pdf, output_pdf):
    # Initialize dictionaries and lists
    count_dict = {'a': [], 't': [], 'c': [], 'g': [], 'n': [], 'I': [], 'D': []}
    positions = []

    # Process the pileup file to extract counts
    with open(pileup_file) as f:
        for line in f:
            ref, pos, refbase, cov, seq, qual = line.split()
            seq = seq.lower()
            refbase = refbase.lower()
            counts = {'a': 0, 't': 0, 'c': 0, 'g': 0, 'n': 0, 'I': 0, 'D': 0}
            depth = 0
            seq = list(seq)
            getdel = False
            getins = False
            while seq != []:
                x = seq.pop(0)
                mod = None
                if x == '.' or x == ',':
                    mod = refbase
                elif x == '+':
                    getins = True
                    digit = ''
                elif x == '-':
                    getdel = True
                    digit = ''
                elif x.isdigit() and (getdel or getins):
                    digit += x
                elif getdel:
                    for j in range(int(digit) - 1):
                        seq.pop(0)
                    mod = 'D'
                    getdel = False
                elif getins:
                    for j in range(int(digit) - 1):
                        seq.pop(0)
                    mod = 'I'
                    getins = False
                elif x == '^':
                    seq.pop(0)
                elif x in ['a', 't', 'c', 'g']:
                    mod = x
                elif x in ['$', '*']:
                    pass
                else:
                    sys.exit("Base not recognised " + x)
                if mod is not None:
                    counts[mod] += 1
                    depth += 1
            if refbase == 'n':
                continue
            if depth >= 10 and counts[refbase] / depth <= 0.85:
                for i in counts:
                    count_dict[i].append(round(counts[i] / depth, 2))
                positions.append(f"{pos} ({refbase})\n[{depth}]")

    # Step 1: Create the plot
    fig, ax = plt.subplots()
    prev_key = 'start'
    for i in count_dict:
        if prev_key == 'start':
            bottom_val = [0] * len(positions)
        else:
            bottom_val = [x + y for x, y in zip(bottom_val, count_dict[prev_key])]
        ax.bar(positions, count_dict[i], 0.5, label=i, bottom=bottom_val)
        prev_key = i

    ax.set_ylabel('Proportions')
    ax.set_title(f'Variants from reference')
    ax.set_xlabel('Position (ref. base)\n[read count]')
    plt.xticks(rotation=90, fontsize=6)
    ax.legend()
    fig.subplots_adjust(bottom=0.3)

    # Save the plot to a temporary PDF
    plot_pdf = f"{sample_name}_temp_plot.pdf"
    fig.set_size_inches(8, 6)
    fig.savefig(plot_pdf, dpi=300)
    plt.close(fig)

    # Step 2: Read the existing PDF and the new plot PDF
    pdf_writer = PdfWriter()
    reader = PdfReader(input_pdf)

    # Add all pages from the existing input PDF
    for page in reader.pages:
        pdf_writer.add_page(page)

    # Append the plot to the output PDF
    plot_reader = PdfReader(plot_pdf)
    pdf_writer.add_page(plot_reader.pages[0])

    # Step 3: Write the final output PDF
    with open(output_pdf, "wb") as out_pdf:
        pdf_writer.write(out_pdf)

    # Clean up temporary plot PDF
    os.remove(plot_pdf)

# Step 4: Call the function with parameters from Snakemake input and output
if __name__ == "__main__":
    pileup_file = snakemake.input.pileup   # pileup file passed as snakemake input
    input_pdf = snakemake.input.pdf     # input PDF passed as snakemake input
    output_pdf = snakemake.output.pdf    # output PDF passed as snakemake output
    sample_name = snakemake.params.sample_name
    add_variant_plot_to_pdf(pileup_file, input_pdf, output_pdf)

