#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import statistics


def process_chromosome(chrs, qc_dir, sample, pp):
    coverage_file = snakemake.input.coverage
    cov1, cov2 = [], []
    with open(coverage_file) as f:
        for line in f:
            cov = int(line.split()[2])
            cov1.append(cov)
            cov2.append(min([100, cov]))

    fig, axs = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
    axs[0].plot(cov1)
    axs[1].plot(cov2)

    ns = 0
    length = 0
    fasta_file = snakemake.input.fasta
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith('>'):
                length += len(line.rstrip())
                ns += line.lower().count('n')

    fig.suptitle(f"{chrs}: all reads\nconsensus {length} bp long\n{ns} Ns", fontsize=10)
    fig.set_size_inches(8, 6)
    # Writing the text report
    report_file = snakemake.output.report
    with open(report_file, 'w') as o:
        o.write('== all reads ==\n')
        o.write(f"length: {length}\n")
        o.write(f"Ns: {ns}\n")
        o.write(f"max coverage: {max(cov1)}\n")
        o.write(f"median coverage: {statistics.median(cov1)}\n")

    # Saving the plots to the PDF
    for ax in axs.flat:
        ax.set(xlabel='position', ylabel='coverage')
        ax.label_outer()

    pp.savefig(fig, dpi=300)  # Corrected this line to use pp.savefig
    plt.close(fig)

def main():

    # Define the PDF file where plots will be saved
    pdf_file = snakemake.output.pdf
    with PdfPages(pdf_file) as pp:
        process_chromosome(snakemake.params.chromosome, snakemake.params.qc_dir, snakemake.params.sample_name, pp)

if __name__ == '__main__':
    main()
