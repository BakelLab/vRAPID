# Cleanup

## Background

This script was built to remove intermedatery output once all steps have been successfully completed.

### Purpose

This script is the ninth step for the vRAPID pipeline for the Bakel Lab. It **must** be run after the successful completion of all previous steps. The script was designed to remove all intermediatery output from the `02_assmebly` subdirectory of each sample, and to compress the pileup output found in the `04_variants` subdirectory.

## Script Usage

### Execution

`cleanup.py -p <samples.csv>`

#### Snakemake

The script can be found in the `Snakefile` found in the `workflow` directory, under the `clean_up` rule. All necessary input, parameters, output, and log files are defined under that rule. Snakemake's log files can always be found in the `log` directory, and identified per sample in the following format:

`11_cleanup.snakemake.log`

### Output

The script will create the an empty file, `<Run-ID>_cleaned.txt`, upon compeletion.

 ```bash
<Run-ID>
├── multi_bamqc
├── multiqc_data
├── multiqc_report.html
├── file_movement_message.txt
├── <Run-ID>_cleaned.txt
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



