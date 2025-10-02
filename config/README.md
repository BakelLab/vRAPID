# vRAPID Pipeline: Configuration Guide

The vRAPID (Virus Reference-based Assembly Pipeline and IDentification) is a comprehensive bioinformatics workflow designed for the assembly, consensus calling, and annotation of viral pathogen genomes. This pipeline is an expansion of the COVID_pipe pipeline developed by Mitchell Sullivan in the Bakel Lab. It utilizes data structures and naming conventions employed by the Center for Research on Influenza Pathogenesis and Transmission (CRIPT) and the Mount Sinai Pathogen Surveillance Program (MS-PSP) at the Icahn School of Medicine at Mount Sinai.

This configuration directory (`snakefile/config`) contains essential files that define the parameters and settings for the pipeline's execution.

------

## 📁 Directory Structure

The `snakefile/config` directory includes the following key files:

- **`config.yaml`**: The primary configuration file that outlines the pipeline's parameters.
- **`samples.csv`**: A CSV file containing metadata for each sample, such as sample IDs and associated information.

------

## 🛠️ Configuration File: `config.yaml`

The `config.yaml` file is structured to provide flexibility and clarity in defining pipeline parameters. Below is an example structure:

```
run_id: "RUN12345"
reference: "reference_genome.fasta"
forward_primer: "forward_primer.fasta"
reverse_primer: "reverse_primer.fasta"
genbank: "genbank_annotations.gff"
primer_set_2kb: "primer_set_2kb.fasta"
primer_set_1_5kb: "primer_set_1_5kb.fasta"
virus: "Influenza_A"
samples: "samples.csv"
ref_fasta_headers: ["header1", "header2", "header3"]
length: 20000
path: "/path/to/data"
```

**Key Parameters:**

- `run_id`: Unique identifier for the sequencing run.
- `reference`: Path to the reference genome file.
- `forward_primer` & `reverse_primer`: Paths to the forward and reverse primer sequences.
- `genbank`: Path to the GenBank annotation file.
- `primer_set_2kb` & `primer_set_1_5kb`: Paths to primer sets used for variant analysis.
- `virus`: Name of the virus being analyzed.
- `samples`: Path to the sample metadata file (CSV format).
- `ref_fasta_headers`: List of reference FASTA headers to be used.
- `length`: Expected length of the assembled genome.
- `path`: Base path for data storage.

------

## 📊 Sample Metadata File: `samples.csv`

The `samples.csv` file should contain metadata for each sample, structured as follows:

```bash
Sample_ID,Sample_Name
S1,Sample1
S2,Sample2
S3,Sample3
```

**Columns:**

- `Sample_ID`: Unique identifier for each sample.
- `Sample_Name`: Common name or label for the sample.

This metadata is utilized throughout the pipeline to organize and process the samples accordingly.

------

## 🔧 Customizing the Configuration

The configuration files are designed to be easily editable to accommodate different datasets and experimental conditions. Modify the `config.yaml` and `samples.csv` files as needed to reflect your specific project requirements.

------

## 📄 License

The vRAPID pipeline is licensed under the MIT License. See the [LICENSE](https://github.com/BakelLab/vRAPID/blob/snakefile/LICENSE) file for details.

------

For more information and updates, please refer to the [vRAPID GitHub repository](https://github.com/BakelLab/vRAPID?utm_source=chatgpt.com).