# Run Genome Annotation

## Background

This script was built to take in the `FASTA` file and to annotate each sample according to virus type.

### Purpose

This script is the fourth step for the vRAPID pipeline for the Bakel Lab. It uses different annotation tools depending on the virus type of the samples being run.

For SARS-CoV-2 samples, the samples are annotated using [`vadr`](https://github.com/ncbi/vadr).

For influenza A & B samples, the samples are annotated using [NCBI's Influenza Annotation Tool](https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/api/help.html).

If an sCoV or an MPX sample is being run through the pipeline, it outputs a file stating that the sample is of either virus type, where no current annotation tool exists. Please note that this is subject to change.

## Script Usage

### Execution

`run_genome_annotation.py -i <sample_ID.fasta> -v <virus-type>`

#### Snakemake

The script can be found in the `Snakefile` found in the `workflow` directory, under the `run_genome_annotation` rule. All necessary input, parameters, output, and log files are defined under that rule. Note that the snakemake rule also outputs a log file for each sample. Snakemake's log files can always be found in the `log` directory, and identified per sample in the following format:

`06_<sample-name>.annotation.snakemake.log`

**SARS-CoV-2 Annotation:**

As stated, `vadr` is used to annotate SARS-CoV-2 samples. To do that, the `vadr` conda environment is activated, the input FASTA file is trimmed, and the specific scripts tied to annotation are called. The scripts can be found within the `workflow/scripts` directory under the following names:

* `vadr_run.py`
* `v-annotate.pl`

**Influenza A/B Annotation:**

The FASTA file is submitted to the virus annotator's API tool. After the job is sucessfully completed, the files are downloaded. After which, the protein IDs are mapped and translated, and the subtype is extracted based on serotype information.

### Output

The script will create the output files within the `02_assembly` directory for each sample.  

**SARS-CoV-2 Annotation:**

 ```bash
<Run-ID>
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

**Influenza A/B Annotation**

```bash
<Run-ID>
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
│   │	│	├──  <sample_ID>.features_cds.fa
│   │	│	├──  <sample_ID>.features_protein.fa
│   │	│	├──  <sample_ID>.features_table.txt
│   │	│	├──  <sample_ID>_genome_annotation.txt
│   │	│	├──  <sample_ID>_subtype.txt
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

