# Multi-Sample QC, Reporting, and Data Submission Workflow

This Snakemake workflow covers **multi-sample QC analysis, report generation, and data submission** to a central database after genome assembly and QC for individual samples.

------

## ⚙️ Requirements

### Software

- Snakemake ≥ 9.0.0
- QualiMap
- MultiQC
- Python 3 and R 4.x with dependencies for custom scripts (`../scripts/`)
- Scripts in `../scripts/`:
  - `04_all-virus-assembly-push.py`
  - `04_generate-run-report.R`
  - `04_generate-workflow-summary.py`

### Input Files

- Individual sample BAMs with read groups:

  ```
  {sample}/01_assembly/{sample}_ref.sorted.rg.bam
  ```

- QC PDFs for each sample:

  ```
  {sample}/03_qualityControl/{sample}_qualityControl.pdf
  ```

- BAM flagstat files:

  ```
  {sample}/03_qualityControl/{sample}_refbam.flagstat
  ```

- Config file: `config.yaml` (contains run_id, qualimap input path, and other metadata)

------

## 🧬 Workflow Steps

### 1️⃣ Run QualiMap Multi-BAM QC

```
rule run_QualiMap_sample
```

- Performs **multi-sample BAM QC** using QualiMap.
- Inputs: BAMs, QC PDFs, flagstat files
- Output: `multi_bamqc/multisampleBamQcReport.html`

**Example command internally:**

```
qualimap multi-bamqc -d {config[qualimap_input]} -r
```

------

### 2️⃣ Generate MultiQC Report

```
rule run_multiqc
```

- Aggregates all QC outputs into a single MultiQC report: `multiqc_report.html`
- Provides an overview of quality metrics across the run.

------

### 3️⃣ Push Data to PathogenDB

```
rule push_data_pathogendb
```

- Submits genome assembly, QC, and BAM statistics to a central **PathogenDB** database.
- Inputs: multi-sample QualiMap report, MultiQC report, BAM flagstat
- Uses Python script `04_all-virus-assembly-push.py` for submission.
- Output: sample-specific log confirming data push: `{sample}.assembly-push.log`

------

### 4️⃣ Generate Run Summary Report

```
rule run_report
```

- Compiles logs from all samples into a **CSV run report**.
- Uses R script `04_generate-run-report.R`
- Output: `{run_id}_run_report.csv`

------

### 5️⃣ Generate Workflow Summary

```
rule summary
```

- Produces a **workflow summary text file** including all packages and versions used.
- Uses Python script `04_generate-workflow-summary.py`
- Output: `workflow_summary.txt`

------

## 📤 Outputs

For the full run:

- **Multi-sample QC reports:**
  - `multi_bamqc/multisampleBamQcReport.html`
  - `multiqc_report.html`
- **Sample submission logs:**
  - `{sample}/04_status/{sample}.assembly-push.log`
- **Run report:**
  - `{run_id}_run_report.csv`
- **Workflow summary:**
  - `workflow_summary.txt`

------

## 🔧 Notes

- Ensure that **assembly**, **variant calling**, and **per-sample QC** workflows are complete before running this stage.
- Logs are saved in `logs/04_*.snakemake.log` for reproducibility.
- MultiQC provides a global overview of the entire sequencing run.
- Data submission scripts require proper configuration in `config.yaml` (database paths, credentials, etc.).