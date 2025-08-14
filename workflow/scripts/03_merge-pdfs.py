from pypdf import PdfWriter

# Get list of input PDFs from Snakemake
pdfs = list(snakemake.input)

# Initialize the PDF merger
merger = PdfWriter()

# Append each PDF to the merger
for pdf in pdfs:
    merger.append(pdf)

# Write out the merged PDF to the output path defined in Snakemake
merger.write(snakemake.output.pdf)
merger.close()

