# Virus Reference-based Assembly Pipeline and IDentification (vRAPID) 

The Virus Reference-based Assembly Pipeline and IDentification (vRAPID) tool runs the assembly, consensus calling and annotation of viral pathogen genomes.

vRAPID relies on data structures and naming conventions used by the [Center for Research on Influenza Pathogenesis and Transmission (CRIPT)](https://www.ceirr-network.org/centers/cript) and the [Mount Sinai Pathogen Surveillance Program (MS-PSP)](https://icahn.mssm.edu/research/genomics/research/microbial-genomics-pathogen-surveillance) at the Icahn School of Medicine at Mount Sinai. vRAPID is an expansion of the [COVID_pipe](https://github.com/mjsull/COVID_pipe) pipeline that was written by Mitchell Sullivan in the Bakel Lab.

#### Prepare the snakemake conda environment

Installation of the required external software packages is largely handled by the pipeline itself, however a conda environment named `snakemake` needs to be present in your environment. We recommend miniconda, which is a free minimal installer for [conda](https://docs.conda.io/en/latest/miniconda.html). Follow the instructions below to start the miniconda installer on Linux. When asked whether the conda environment should automatically be initialized, select 'yes'. Note that Snakemake requires the channel_priority to be set to strict. The post-installation commands to apply this setting are included in the post-installation selection below.

```bash
# Start miniconda installation
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh

# Post-installation commands to enforce strict channel_priority (required for Snakemake)
conda config --set auto_activate_base false
conda config --set channel_priority strict
```

On the Mount Sinai [Minerva](https://labs.icahn.mssm.edu/minervalab/minerva-quick-start/) high-performance computer cluster we additionally recommend setting the conda environment `pkgs` and `envs` directories to your `work` folder, which has a much larger storage allocation (100-200 Gb) than the home folders (5 Gb). The conda configuration commands to set the directories are as follows:

```bash
mkdir -p /sc/arion/work/$USER/conda/pkgs
mkdir -p /sc/arion/work/$USER/conda/envs
conda config --add pkgs_dirs "/sc/arion/work/$USER/conda/pkgs"
conda config --add envs_dirs "/sc/arion/work/$USER/conda/envs"
```

After installation and configuration of conda/miniconda, the following 'conda create' command can be used to set up the required 'snakemake' environment.

```bash
conda create -c conda-forge -c bioconda -n snakemake 'snakemake>=9.0'
conda activate snakemake
pip install snakemake-executor-plugin-lsf
conda deactivate
```

------

## 📦 Overview

The vRAPID pipeline consists of four main modules:

1. **Assembly & Primer Trimming**
2. **Variant Calling & Consensus Polishing**
3. **Per-Sample Quality Control (QC)**
4. **Multi-Sample QC, Reporting, and Database Submission**

**Workflow diagram (conceptual):**

```
Raw FASTQs
   │
   ▼
Trim Primers (cutadapt)
   │
   ▼
QC of Reads (FASTQC)
   │
   ▼
Align Reads to Reference (minimap2)
   │
   ▼
Sorted & Indexed BAM → Read Groups → Deduplication (Picard)
   │
   ▼
Variant Calling (Clair3)
   │
   ▼
Generate Preconsensus → Polish Consensus → Concatenate FASTA
   │
   ▼
Coverage & Pileup Analysis → Variable Bases → Masking Sites
   │
   ▼
Per-Sample QC Plots (coverage, primers, intrahost variants)
   │
   ▼
Kraken2 Taxonomic Classification → Plot Taxonomic Breakdown
   │
   ▼
Merge QC PDFs → Flagstat
   │
   ▼
Multi-Sample QC (QualiMap) → MultiQC
   │
   ▼
Data Submission to PathogenDB
   │
   ▼
Run Summary & Workflow Summary Reports
```

------

## ⚙️ Requirements

### Software

- **Snakemake ≥ 9.0.0**
- **Python 3** (with required packages in `../envs/env.yml`)
- **R 4.x** (with `ggplot2`, `cowplot`, etc.)
- **cutadapt** – primer trimming
- **FASTQC** – read quality check
- **minimap2** – alignment to reference
- **samtools** – depth, flagstat, indexing
- **Picard** – read groups, deduplication, alignment metrics
- **Clair3** – variant calling
- **bcftools** – consensus polishing
- **bamtools** – BAM splitting
- **Kraken2** – taxonomic classification
- **QualiMap** – multi-sample BAM QC
- **MultiQC** – aggregated QC reporting

### Input Files

- Paired-end FASTQs: `{sample}/00_fastqs/{sample}_1.fastq.gz`, `{sample}/00_fastqs/{sample}_2.fastq.gz`
- Reference genome FASTA: `{config[reference_genome]}`
- Primer files: `{config[forward_primer]}`, `{config[reverse_primer]}`
- Sample metadata CSV: `{config[samples]}` (must contain `Sample_ID`)
- Config file: `config.yaml` with run-specific settings, see the [config](config/config.yaml) file here for a detailed parameter set-up.

------

## 🧬 Workflow Modules

### 1️⃣ Assembly & Primer Trimming

- **Trim primers:** `cutadapt`
- **QC reads:** `FASTQC`
- **Align reads:** `minimap2` → sorted BAM → add read groups → mark duplicates
- **Polish consensus sequences:** `Clair3`, `bcftools`
- **Outputs:** polished FASTA, sorted BAM, consensus VCF

### 2️⃣ Variant Calling & Consensus Polishing

- Generate **pileup files** per chromosome
- Process pileup to **variable bases TSV**
- Mask low-coverage regions and generate **preconsensus sequences**
- Polished consensus sequences per chromosome
- Concatenate chromosomes into **sample-level FASTA**

### 3️⃣ Per-Sample Quality Control (QC)

- **Coverage computation** using `samtools depth`
- Generate **QC PDFs**:
  - Coverage plots
  - Primer depth plots
  - Variant visualization
  - Intrahost variant plots
- **Taxonomic classification:** `Kraken2`
- Merge all QC plots into **single PDF per sample**
- Run `samtools flagstat` for BAM statistics

### 4️⃣ Multi-Sample QC, Reporting, and Submission

- **Multi-sample BAM QC:** `QualiMap`
- **MultiQC:** aggregate all QC metrics into a single report
- **Data submission:** scripts push genome assemblies and QC logs to PathogenDB
- **Run summary:** CSV summarizing all samples and processing
- **Workflow summary:** text file with software versions and pipeline overview

------

## 📤 Outputs

For each sample:

- `{sample}/01_assembly/` – polished FASTA, BAMs, and consensus VCF
- `{sample}/02_variants/` – pileup files, variable bases TSV
- `{sample}/03_qualityControl/` – QC PDFs, flagstat, taxonomic reports
- `{sample}/04_status/` – logs for PathogenDB submission

Global outputs:

- `multi_bamqc/multisampleBamQcReport.html` – QualiMap multi-sample report
- `multiqc_report.html` – aggregated QC report
- `{run_id}_run_report.csv` – run summary
- `workflow_summary.txt` – packages and versions used

------

## 🔧 Configuration

- Edit `config.yaml` to set:
  - `samples` → path to sample metadata CSV
  - `reference_genome` → reference FASTA
  - `ref_fasta_headers` → chromosomes/contigs
  - `forward_primer` & `reverse_primer`
  - `depth` → minimum coverage thresholds for masking
  - `kraken_db` → Kraken2 database path
  - `qualimap_input` → TXT file that is tab-delimited of sample and BAM paths for multi-sample QC
  - `run_id` → unique run identifier

------

## 📝 Notes

- **Modular design:** Each module (assembly, variant calling, QC, reporting) can be run independently.
- **Incremental PDFs:** QC plots are appended sequentially to produce a single comprehensive PDF per sample.
- **Logging:** Every rule generates logs under `logs/` for reproducibility.
- **Reproducibility:** Conda environments are specified for each rule via `../envs/env.yml`.