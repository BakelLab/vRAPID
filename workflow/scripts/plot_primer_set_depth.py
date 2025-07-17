import os
import statistics
import matplotlib.pyplot as plt
from PyPDF2 import PdfReader, PdfWriter
from collections import defaultdict
import re

def plot_primer_depths(primer_file, sample_folder, input_pdf, output_pdf, qc_dir, chrs_list):
    # Initialize the PDF writer
    pdf_writer = PdfWriter()

    # Add original input PDF pages if available
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

        # Parse primers and group by base name
        primers_by_base = defaultdict(lambda: {'LEFT': {}, 'RIGHT': {}})

        with open(primer_file) as f:
            f.readline()  # Skip header
            for line in f:
                parts = line.split()
                primer_name = parts[0]
                pos = int(parts[7])

                # Extract base, side, and version
                match = re.match(r"(.+)_(LEFT|RIGHT)(?:_v(\d+))?$", primer_name)
                if not match:
                    continue
                base = match.group(1)  # e.g., covid19_1.5kb_7
                side = match.group(2)  # LEFT or RIGHT
                version = match.group(3) or "0"  # "0" means no version

                primers_by_base[base][side][version] = pos

        # Sort bases by numeric amplicon number
        def sort_key(base_name):
            num_match = re.search(r'(\d+)', base_name)
            return int(num_match.group(1)) if num_match else 0

        labels = []
        vals = []
        xnums = []
        num = 1

        for base in sorted(primers_by_base.keys(), key=sort_key):
            sides = primers_by_base[base]
            left_versions = sides['LEFT']
            right_versions = sides['RIGHT']

            if not left_versions or not right_versions:
                continue  # Skip if incomplete

            # Sort RIGHT primers by version for consistency
            for r_ver in sorted(right_versions.keys(), key=int):
                r_pos = right_versions[r_ver]
                # Match LEFT version if exists, else fallback
                if r_ver in left_versions:
                    l_pos = left_versions[r_ver]
                    label = f"{base}_v{r_ver}" if r_ver != "0" else base
                else:
                    # If only one LEFT exists, pair with all RIGHT versions
                    if len(left_versions) == 1:
                        l_ver, l_pos = next(iter(left_versions.items()))
                        label = f"{base}_v{r_ver}" if r_ver != "0" else base
                    else:
                        continue  # Ambiguous case: multiple LEFT and no match

                if l_pos and r_pos:
                    if l_pos - 1 >= len(cov1) or r_pos > len(cov1):
                        median_val = 0
                    else:
                        median_val = statistics.median(cov1[l_pos-1:r_pos])
                    vals.append(median_val)
                    labels.append(label)
                    xnums.append(num)
                    num += 1

        if not vals:
            print(f"No valid primer pairs for {chrs}, skipping plot.")
            continue

        # Dynamic figure width: 0.4 inch per bar, minimum 12 inches
        fig_width = max(12, len(labels) * 0.4)
        fig_height = 6

        # Plotting
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        ax.bar(xnums, vals, color="steelblue")

        plt.xticks(xnums, labels, rotation=90, fontsize=6)
        ax.set(xlabel="Primer Set", ylabel="Median Depth")
        ax.tick_params(axis='both', which='major', labelsize=6)
        fig.suptitle(f"Median Depth for Primer Set - {chrs}\n{os.path.basename(primer_file)}", fontsize=10)

        # Adjust layout for labels
        fig.tight_layout()
        fig.subplots_adjust(top=0.88, bottom=0.35)

        # Save plot to temp file
        plot_pdf = f"{sample_folder}/temp_{chrs}_primer_depth_plot.pdf"
        plt.savefig(plot_pdf, dpi=300)
        plt.close(fig)

        # Add this plot to the PDF writer
        plot_reader = PdfReader(plot_pdf)
        pdf_writer.add_page(plot_reader.pages[0])

        # Clean up temp plot
        os.remove(plot_pdf)

    # Write combined PDF
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
        chrs = [chrs_param]
    else:
        chrs = chrs_param

    plot_primer_depths(primer_file, sample_folder, input_pdf, output_pdf, qc_dir, chrs)

