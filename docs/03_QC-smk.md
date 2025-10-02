# Quality Control (QC) Workflow for Sequencing Data

This Snakemake workflow performs **comprehensive quality control** for sequencing data after genome assembly and variant calling. It integrates coverage analysis, variant plotting, primer depth visualization, intrahost variant assessment, and taxonomic classification.

------

## ⚙️ Requirements

### Software

- Snakemake ≥ 9.0.0
- Samtools (depth, flagstat)
- Kraken2 (documentation)
- R with `ggplot2` and other dependencies for plotting scripts
- Python 3 with required packages (as specified in `../envs/env.yml`)
- Scripts in `../scripts/`:
  - `03_generate-QC-report.py`
  - `03_plot-variants.py`
  - `03_plot-primer-set-depth.py`
  - `03_plot-intrahosts.R`
  - `03_plot-taxonomic-breakdown.R`
  - `03_merge-pdfs.py`

### Input Files

- Polished consensus FASTA for each chromosome:

  ```
  {sample}/01_assembly/{chromosomes}_polished.fasta
  ```

- Sorted BAM with read groups:

  ```
  {sample}/01_assembly/{sample}_ref.sorted.rg.bam
  ```

- Variable bases TSV from pileup processing:

  ```
  {sample}/02_variants/{chromosomes}_variable_bases.tsv
  ```

- Pileup files:

  ```
  {sample}/02_variants/{chromosomes}_pileup
  ```

- Original FASTQs for Kraken2:

  ```
  {sample}/00_fastqs/{sample}_1.fastq.gz
  {sample}/00_fastqs/{sample}_2.fastq.gz
  ```

------

## 🧬 Workflow Steps

### 1️⃣ Compute coverage per chromosome

```
rule compute_coverage
```

- Uses `samtools depth` to calculate coverage at each position.
- Output: `{chromosomes}_coverage.txt`

------

### 2️⃣ Generate QC report (coverage plot)

```
rule generate_qc_report
```

- Generates a PDF of coverage per FASTA header.
- Logs QC metrics into a text report.

------

### 3️⃣ Plot variants on coverage

```
rule add_variant_plot_to_pdf
```

- Adds variant information from pileup files to the QC PDF.

------

### 4️⃣ Plot primer depth

```
rule add_primers_plot_to_pdf
```

- Adds primer depth plots to the QC report for quality assessment.

------

### 5️⃣ Intrahost variant visualization

```
rule plot_intrahosts
```

- Uses R script to plot intrahost variants per chromosome.
- Outputs updated PDFs and annotated TSVs.

------

### 6️⃣ Taxonomic classification

```
rule run_kraken
```

- Runs Kraken2 on raw paired-end FASTQs to classify reads.
- Outputs Kraken2 report and detailed classification files.

------

### 7️⃣ Plot taxonomic breakdown

```
rule plot_taxonomic_breakdown
```

- Creates a PDF summarizing the taxonomic composition of reads.

------

### 8️⃣ Merge PDFs

```
rule merge_pdfs
```

- Combines coverage, variant, primer, intrahost, and taxonomic plots into a single PDF per sample.

------

### 9️⃣ BAM flag statistics

```
rule run_flagstat
```

- Uses `samtools flagstat` to summarize BAM statistics.
- Output: `{sample}_refbam.flagstat`

------

## 📤 Outputs

For each sample:

- **Coverage files:** `{chromosomes}_coverage.txt`
- **QC PDFs:** `{sample}_qualityControl.pdf`
- **Variant annotations:** `{chromosomes}_variants_labels.tsv`
- **Taxonomic report:** `taxonomic-breakdown.pdf`
- **BAM statistics:** `{sample}_refbam.flagstat`

------

## 🔧 Notes

- Ensure prior completion of assembly and variant calling workflows.
- Logs are saved per step in `logs/{sample}/03_qualityControl/`.
- PDF reports are incremental; each step appends to the previous PDF.
- Kraken2 requires a valid database specified in `config.yaml` (`kraken_db`).