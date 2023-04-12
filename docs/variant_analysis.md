# Variant Analysis

## Background

This script was built to take in the consensus FASTA file and to preform intra-host variant analysis on the virus library.

### Purpose

This script is the second step for the vRAPID pipeline for the Bakel Lab. It uses `samtools mpileup` to create a text pileup output for the `bam` file obtained from the assembly step, more details can be found [here](). 

## Script Usage

### Execution

`variant_analysis.py -rd <path-to-repo-dir> -i <sample_folder> -r <reference-genome-FASTA> -ps_2 <primer-set-2kb> -ps_1_5 <primer-set_1_5kb> `

#### Snakemake

The script can be found in the `Snakefile` found in the `workflow` directory, under the `variant_analysis` rule. All necessary input, parameters, output, and log files are defined under that rule. Note that the snakemake rule also outputs a log file for each sample. Snakemake's log files can always be found in the `log` directory, and identified per sample in the following format:

`02_<sample-name>.aseembly.snakemake.log`

#### Directory structure

For the script to be successfully excuted, the following directory structure is required:

```bash
<Run-ID>
├── <sample_ID>
│   ├── <01_fastqs>
│   │   ├── <sample_ID>_1.fastq.gz
│   │	│	├── <sample_ID>_2.fastq.gz
│   ├── <02_assembly>
│   │   ├── <sample_ID>_ref.bam
│   │	│	├── <sample_ID>_ref.bam.bai
│   │	│	├── <sample_ID>.fasta
```

### Output

The script will create a new subdirectory with the `pileup` file(s) and the TSV containing the final results of the variant analysis.

 ```bash
<Run-ID>
├── <sample_ID>
│   ├── <01_fastqs>
│   │   ├── <sample_ID>_1.fastq.gz
│   │	│	├── <sample_ID>_2.fastq.gz
│   ├── <02_assembly>
│   │   ├── <sample_ID>.fasta
│   │	│	├── <sample_ID>_ref.bam
│   │	│	├── <sample_ID>_ref.bam.bai
│   │	│	├── *_pilon.*
│   │	│	├── <reference-FASTA-header>_shovill
│   │	│	├── prokka
│   ├── <04_variants>
│   │   ├── <sample_ID>.<reference-fasta-header>_variable_bases.tsv
│   │	│	├── <reference-fasta-header>.final.fna
│   │	│	├── <reference-fasta-header>_pileup.gz
│   │	│	├── pileup.gz
 ```

