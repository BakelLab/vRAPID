# Run Pipeline

## Background

This script was built to take in Illumina paired end FASTQ files and to assemble them against a reference genome, outputting a final consensus file.

### Purpose

This script is the first step for the vRAPID pipeline for the Bakel Lab. It trims the reads using `cutadapt`, before running the command line version of `FASTQC` on the `cutadapt` output, then aligining the reads against the reference genome using `minimap2`. The script then uses `pilon` to imporve the assembly quality. Lastly, `prokka` is used for a rapid annotation and `shovill` is used for a rapid assembly.

## Script Usage

### Execution

`run_pipeline.py -rd <path-to-repo-dir> -i <sample_folder> -r1 <R1-suffix> -r2 <R2-suffix> -r <reference-genome-FASTA> -pf <forward-primer-file> -pr <reverse-primer-file> -g <reference-genbankfile> -l <reference-genome-length> -headers <reference-FASTA-headers> `

#### Snakemake

The script can be found in the `Snakefile` found in the `workflow` directory, under the `assemble` rule. All necessary input, parameters, output, and log files are defined under that rule. Note that the snakemake rule also outputs a log file for each sample. Snakemake's log files can always be found in the `log` directory, and identified per sample in the following format:

`02_<sample-name>.aseembly.snakemake.log`

#### Directory structure

For the script to be successfully excuted, the following directory structure is required:

```bash
<Run-ID>
├── <sample_ID>
│   ├── <sample_ID>
│   │   ├── <sample_ID>_1.fastq.gz
│   │	│	├── <sample_ID>_2.fastq.gz
```

#### *Note*

If excuted within `snakemake`, the very first rule of the pipeline renames the subdirectory holding the FASTQ files, thus the structure would be slighlty different, see below:

```bash
<Run-ID>
├── <sample_ID>
│   ├── <01_fastqs>
│   │   ├── <sample_ID>_1.fastq.gz
│   │	│	├── <sample_ID>_2.fastq.gz
```

### Output

The script will create a new subdirectory with the final consensus, the `bam` files, `pilon`, `prokka`, and `shovill` output files.

 ```bash
 <Run-ID>
 ├── <sample_ID>
 │   ├── <01_fastqs>
 │   │   ├── <sample_ID>_1.fastq.gz
 │   │	│	├── <sample_ID>_2.fastq.gz
 │   ├── <02_assembly>
 │   │   ├── <sample_ID>.fasta
 │   │	│	├── <sample_ID>_ref.bam.
 │   │	│	├── <sample_ID>_ref.bam.bai
 │   │	│	├── *_pilon.*
 │   │	│	├── <reference-FASTA-header>_shovill
 │   │	│	├── prokka
 ```

