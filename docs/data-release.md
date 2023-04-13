# Data Release

## Background

This script was built to prepare the samples for data release back to the submitter after going through the assembly process with the vRAPID pipeline.

### Purpose

This script is the last step for the vRAPID pipeline for the Bakel Lab. It queries PathogenDB's `tCEIRS_assemblies`, `tCEIRS_isolates`, and `tCEIRS_extracts` table to provide a final CSV with the assembly run's status and statistics. Concurrently, it loops through the sample directories that were run, compressing the final consensus, FASTQ files, bam files, TSV file obtained as a result of the variant analysis step, and the quality control PDF for all of the samples into one file.

## Script Usage

### Execution

`data-release.py -p <samples-run-csv-file> -r <Run-ID> -v <virus-type> -f <reference-FASTA-header>`

#### Snakemake

The script can be found in the `Snakefile` found in the `workflow` directory, under the `data_release_prep` rule. All necessary input, parameters, output, and log files are defined under that rule. Note that the snakemake rule also outputs a log file for each sample. Snakemake's log files can always be found in the `log` directory, and identified per sample in the following format:

`13_datarelease.snakemake.log`

### Output

The script will create the output files in the current working directory. 

```bash
<Run-ID>
├── multi_bamqc
├── multiqc_data
├── multiqc_report.html
├── file_movement_message.txt
├── <Run-ID>_cleaned.txt
├── <Run-ID>_report.csv
├── <Run-ID>_<collaborator-ID>_assemblies.csv
├── <Run-ID>_<collaborator-ID>.tar.gz
├── <sample_ID>
│   ├── <01_fastqs>
│   │  ├── <sample_ID>_1.fastq.gz
│   │	│	├── <sample_ID>_2.fastq.gz
│   ├── <02_assembly>
│   │  ├── <sample_ID>.fasta
│   │	│	├── <sample_ID>_ref.bam
│   │	│	├── <sample_ID>_ref.bam.bai
│   │	│	├── *_pilon.*
│   │	│	├── <reference-FASTA-header>_shovill
│   │	│	├── prokka
│   │	│	├── <sample_ID>_annotation
│   │	│	├──  <sample_ID>_genome_annotation.txt
│   ├── <04_variants>
│   │  ├── <sample_ID>.<reference-fasta-header>_variable_bases.tsv
│   │	│	├── <reference-fasta-header>.final.fna
│   │	│	├── <reference-fasta-header>_pileup.gz
│   │	│	├── pileup.gz
│   ├── <03_qualityControl>
│   │  ├──  <reference-FASTA-header>_coverage.txt
│   │  ├──  <reference-FASTA-header>_report.txt
│   │  ├── <sample_ID>_kraken
│   │  ├── <sample_ID>_kraken_input.fastq.gz
│   │  ├── <sample_ID>_kraken_report.out
│   │  ├── <sample_ID>_quality_control.pdf
│   │  ├── <sample_ID>_refbam.flagstat
```

