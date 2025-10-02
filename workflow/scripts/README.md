# vRAPID Workflow Scripts

This directory contains all the scripts used in the vRAPID Snakemake workflow for viral genome assembly, variant calling, quality control, and reporting.

## Directory Structure

- **01_***: Initial assembly, preprocessing, and variant calling scripts.
- **02_***: Scripts for processing pileup files and handling variant information.
- **03_***: Quality control and visualization scripts, including coverage, intrahost variants, primer depth, and taxonomic breakdown plots.
- **04_***: Post-analysis scripts for pushing data to PathogenDB and generating run reports.

## Script Usage

Scripts in this directory are called by the Snakemake rules and are generally not intended to be run manually. Each script typically accepts command-line arguments specified by the corresponding Snakemake rule.

### Key Script Categories

1. **01_\* — Assembly and Variant Calling**
   - `01_run_pipeline.py` and related scripts: Prepare input FASTQs and reference genome, perform genome assembly, generate BAM files.
   - `01_variant_analysis.py`: Generate pileup files and extract intra-host variant information.
2. **02_\* — Pileup Processing**
   - `02_process-pileup.py`: Process pileup files to produce tables of variable bases for downstream analysis.
   - Additional scripts for handling primer sets and masked regions.
3. **03_\* — Quality Control**
   - `03_generate-QC-report.py`: Coverage analysis per FASTA header and initial QC reporting.
   - `03_plot-variants.py`: Overlay variant information on QC plots.
   - `03_plot-primer-set-depth.py`: Plot primer depth across the genome.
   - `03_plot-intrahosts.R`: Visualize intra-host variants.
   - `03_plot-taxonomic-breakdown.R`: Generate taxonomic breakdown plots from Kraken2 output.
   - `03_merge-pdfs.py`: Combine QC plots, intrahost plots, and taxonomic plots into a final report.
4. **04_\* — Post-analysis and Reporting**
   - `04_all-virus-assembly-push.py`: Push genome assembly and QC data to PathogenDB.
   - `04_generate-run-report.R`: Generate run-level summary report.
   - `04_generate-workflow-summary.py`: Create workflow-level summary of tools, versions, and results.

## Requirements

- Python (≥3.9 recommended)
- R (for plotting and QC scripts)
- Conda environments as specified in the workflow’s `envs/` directory
- Dependencies are managed per script via the workflow’s Conda environments.

## Notes

- Scripts follow a modular design to enable flexible execution within the Snakemake workflow.
- All file paths and parameters are dynamically set through the Snakemake rules.
- Users should refer to the Snakemake `Snakefile` for the exact order of execution and interdependencies.