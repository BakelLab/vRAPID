# vRAPID Pipeline Documentation

## Overview

vRAPID is a Snakemake-based pipeline designed for **viral genome assembly, variant analysis, quality control, and annotation**. It supports multiple viruses, integrates primer information, and generates both per-sample and run-level reports.

The pipeline automates the following major steps:

1. **Genome Assembly** – Assemble viral genomes using reference-guided methods.
2. **Variant Analysis** – Detect intra-host variants and generate pileup-based analyses.
3. **Quality Control** – Compute coverage, generate QC plots, and run taxonomic classification using Kraken2.
4. **Annotation** – Annotate viral genomes based on GenBank reference files.
5. **Reporting** – Generate MultiQC, per-sample, and run-level reports; push results to PathogenDB.

------

## Directory Structure

```bash
vRAPID/
├── workflow/            # Snakemake rules and scripts
├── config/              # Configuration YAML files for runs
├── db/                  # Reference genomes, primer sets, GenBank annotations
├── envs/                # Conda environments for reproducibility
└── docs/                # Documentation
```

------

## Key Features

- Supports multiple viral references and primer sets
- Variant calling with Clair3
- Detailed per-sample and per-run QC reporting (PDFs, plots, MultiQC)
- Integration with Kraken2 for taxonomic profiling
- Automated upload of results to PathogenDB
- Flexible configuration via `config.yaml`

------

## Getting Started

1. **Configure your run** – Edit `config/config.yaml` with sample info, reference genomes, primers, and paths.
2. **Set up Conda environments** – Use the provided `envs/` environments to ensure reproducibility.
3. **Run Snakemake** – Execute the workflow

------

## Documentation

- [Configuration files](../config) – Define samples, references, and pipeline parameters.
- [Reference databases](../db) – Genome references, primer sets, GenBank annotations, Kraken2 database.
- [Scripts & Rules](../workflow) – Custom Python and R scripts used by Snakemake rules.

------

## License

This project is licensed under the same terms as the main `vRAPID` repository. See the `LICENSE` file for details.