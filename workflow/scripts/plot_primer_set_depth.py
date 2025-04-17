import os
import statistics
import matplotlib.pyplot as plt
from PyPDF2 import PdfReader, PdfWriter

def plot_primer_depths(primer_file, sample_folder, input_pdf, output_pdf, qc_dir, chrs_list):
    # Initialize the PDF writer
    pdf_writer = PdfWriter()

    # First, add the original input PDF pages
    if os.path.exists(input_pdf):
        reader = PdfReader(input_pdf)
        for page in reader.pages:
            pdf_writer.add_page(page)

    for chrs in chrs_list:
        # Load coverage data for the current chromosome
        cov1 = []
        coverage_file = os.path.join(qc_dir, f"{chrs}_coverage.txt")
        if not os.path.exists(coverage_file):
            print(f"Warning: Coverage file for {chrs} not found: {coverage_file}")
            continue

        with open(coverage_file) as f:
            for line in f:
                cov = int(line.split()[2])
                cov1.append(cov)

        # Process the primer file
        with open(primer_file) as f:
            labels = []
            vals = []
            xnums = []
            f.readline()  # Skip header
            num = 1
            for line in f:
                parts = line.split()
                if parts[0].endswith("LEFT"):
                    label = parts[0][:-5]
                    l = int(parts[7])
                elif parts[0].endswith("RIGHT"):
                    r = int(parts[7])
                    # Handle edge case if l or r exceed coverage length
                    if l-1 >= len(cov1) or r > len(cov1):
                        median_val = 0
                    else:
                        median_val = statistics.median(cov1[l-1:r])
                    vals.append(median_val)
                    xnums.append(num)
                    labels.append(label)
                    num += 1

        # Plotting
        fig, ax = plt.subplots()
        ax.bar(xnums, vals)
        plt.xticks(xnums, labels, rotation=90, fontsize=6)
        ax.set(xlabel="Primer Set", ylabel="Median Depth")
        ax.tick_params(axis='both', which='major', labelsize=6)
        fig.suptitle(f"Median Depth for Primer Set - {chrs}\n{os.path.basename(primer_file)}", fontsize=10)
        fig.tight_layout()
        fig.subplots_adjust(top=0.85)
        fig.set_size_inches(8, 6)

        # Save plot to temp file
        plot_pdf = f"temp_{chrs}_primer_depth_plot.pdf"
        plt.savefig(plot_pdf, dpi=300)
        plt.close(fig)

        # Add this plot to the PDF writer
        plot_reader = PdfReader(plot_pdf)
        pdf_writer.add_page(plot_reader.pages[0])

        # Clean up temp plot
        os.remove(plot_pdf)

    # Finally, write all pages to the output PDF
    with open(output_pdf, "wb") as out_pdf:
        pdf_writer.write(out_pdf)

    print(f"Combined PDF saved to: {output_pdf}")

if __name__ == "__main__":
    primer_file = snakemake.input.primers
    sample_folder = snakemake.params.sample_folder
    input_pdf = snakemake.input.pdf
    output_pdf = snakemake.output.pdf
    qc_dir = snakemake.params.qc_dir
    chrs_param = snakemake.params.chromosomes
    if isinstance(chrs_param, str):
        chrs = [chrs_param]  # Treat single string as one chromosome
    else:
        chrs = chrs_param

    plot_primer_depths(primer_file, sample_folder, input_pdf, output_pdf, qc_dir, chrs)

