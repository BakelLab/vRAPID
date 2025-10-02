# Readme: Consensus Genome and Variant Calling

This repository contains a **Snakemake workflow** for processing paired-end sequencing data into consensus genomes with quality control, read alignment, variant calling, consensus polishing, and downstream metrics collection.

------

## 📌 Features

This pipeline automates the following steps:

1. **Primer trimming** with Cutadapt
2. **Read quality control** with FastQC
3. **Read alignment** to a reference genome using Minimap2
4. **Sorting and indexing** BAM files with Samtools
5. **Adding read groups** with Samtools/Picard
6. **Variant calling** with Clair3 (dockerized container)
7. **Coverage-based masking** and consensus sequence generation
8. **Consensus polishing** with BCFtools
9. **Concatenating FASTA sequences** across all reference chromosomes
10. **BAM splitting and insert-size metrics** with Bamtools and Picard
11. **Duplicate marking and alignment metrics** using Picard

The final outputs are **consensus FASTA files**, **VCF variant calls**, and **QC metrics** for each sample.

------

## ⚙️ Requirements

### Software

This workflow requires **Snakemake ≥ 9.0.0** and uses Conda/Docker for reproducibility.
Dependencies include:

- Cutadapt
- FastQC
- [Minimap2](https://github.com/lh3/minimap2)
- Samtools
- [Clair3](https://github.com/HKU-BAL/Clair3) (Docker image: `hkubal/clair3:v1.2.0`)
- BCFtools
- [Bamtools](https://github.com/pezmaster31/bamtools)
- Picard
- Python scripts for masking and preconsensus generation (provided in `workflow/scripts/`)

### Environment

- **Conda environment file:** `workflow/envs/env.yml`
- **Docker** (for Clair3 execution)

------

## 📂 Input Files

### 1. Configuration file (`config.yaml`)

Define all references, primers, and thresholds. Example:

```
samples: "samples.csv"
reference_genome: "data/reference.fasta"
ref_fasta_headers: ["Segment1", "Segment2", "Segment3"]
forward_primer: "data/primers/forward.fa"
reverse_primer: "data/primers/reverse.fa"
depth: 10
```

### 2. Sample sheet (`samples.csv`)

A CSV file listing all samples:

```
Sample_ID
SampleA
SampleB
SampleC
```

### 3. Input FASTQ files

For each sample listed in `samples.csv`, place paired FASTQs under:

```
SampleA/00_fastqs/SampleA_1.fastq.gz
SampleA/00_fastqs/SampleA_2.fastq.gz
```

------

## 🧬 Workflow Overview

### Step 1: Trim primers

```
rule run_cutadapt
```

- Removes primer sequences from raw reads.

### Step 2: Quality check

```
rule run_fastqc
```

- Runs FastQC on trimmed reads.

### Step 3: Align reads

```
rule minimap2_clip
```

- Aligns reads to reference genome with Minimap2.
- Removes soft-clipped reads and generates indexed BAM.

### Step 4: Sort and index BAM

```
rule sorted_bam
```

### Step 5: Add read groups

```
rule group_reads_in_bam
```

### Step 6: Variant calling with Clair3

```
rule run_clair3
```

- Uses Dockerized Clair3 to generate a compressed VCF.

### Step 7: Coverage masking

```
rule get_masking_sites
```

- Generates depth-based mask sites via `01_generate-depth-masking.py`.

### Step 8: Pre-consensus sequence

```
rule generate_preconsensus_sequence
```

- Combines VCF and mask regions to create a draft consensus.

### Step 9: Consensus polishing

```
rule polish_consensus
```

- Normalizes variants and generates polished consensus FASTA.

### Step 10: Concatenate FASTAs

```
rule concatenate_FASTAs
```

- Combines all per-chromosome FASTAs into a single consensus per sample.

### Step 11: BAM splitting

```
rule run_bamtools_split
```

- Splits BAMs by reference chromosome.

### Step 12–15: QC metrics with Picard

- Insert-size metrics (`run_picard_insert_metrics`)
- Read group fixing (`run_picard_read_groups`)
- Duplicate marking (`run_picard_duplicates`)
- Alignment summary (`run_picard_alignment_metrics`)

------

## 📤 Outputs

For each sample, the pipeline generates:

- **Consensus FASTA:**

  ```
  {sample}/01_assembly/{sample}.fasta
  ```

- **Variants (VCF):**

  ```
  {sample}/01_assembly/clair3/merge_output.vcf.gz
  ```

- **QC Reports:**

  - FastQC HTML reports
  - Insert size metrics (TXT, PDF)
  - Alignment summary metrics (TXT)
  - Duplicate metrics (TXT)

- **Intermediate files:**
  Sorted BAMs, read-group BAMs, polished FASTAs, masking sites.

------

## 📎 Notes

- Ensure that `config.yaml` and `samples.csv` are properly set before execution.
- Uses strict conda/Docker isolation for reproducibility.