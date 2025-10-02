# Pileup Processing Workflow

This Snakemake workflow performs **variant detection and processing** for assembled consensus genomes using BAM files and reference genomes. It relies on `samtools mpileup` and a custom Python script to extract variable bases for downstream analysis.

------

## ⚙️ Requirements

### Software

- Snakemake ≥ 9.0.0
- Samtools
- Python 3 with required packages for custom scripts (specified in `../envs/env.yml`)
- `02_process-pileup.py` (provided in `workflow/scripts/`)

### Input Files

1. **Reference genome** – as specified in `config.yaml` (`reference_genome`)

2. **Chromosome headers** – as specified in `config.yaml` (`ref_fasta_headers`)

3. **Sample FASTA** – generated from the assembly workflow:

   ```
   {sample}/01_assembly/{sample}.fasta
   ```

4. **Sorted BAM with read groups** – generated from previous workflow steps:

   ```
   {sample}/01_assembly/{sample}_ref.sorted.rg.bam
   ```

5. **Clair3 VCF** – variant calls used in processing:

   ```
   {sample}/01_assembly/clair3/merge_output.vcf.gz
   ```

------

## 🧬 Workflow Steps

### 1️⃣ Generate Pileup

```
rule run_mpileup
```

- Uses `samtools mpileup` to generate pileup files for each chromosome.
- Inputs: sample FASTA and BAM.
- Outputs: `{chromosome}_pileup` file per sample.
- Logs: saved in `logs/{sample}/02_variants/01_{chromosomes}_mpileup.snakemake.log`.

**Example command executed internally:**

```
samtools mpileup -f reference.fasta -r Segment1 SampleA_ref.sorted.rg.bam -o Segment1_pileup
```

------

### 2️⃣ Process Pileup

```
rule process_pileup
```

- Converts pileup files into a **variable bases table (TSV)**.
- Incorporates variant information from the Clair3 VCF.
- Uses the custom script `02_process-pileup.py` with configurable parameters.

**Inputs:**

- Pileup file from previous step
- Clair3 VCF file

**Outputs:**

- `{chromosome}_variable_bases.tsv` in `02_variants/` directory

**Parameters:**

- `min_ratio`: minimum allele fraction to report (default 0.9)
- `primer_set`: primer information from `config.yaml`

**Logs:**

- Saved in `logs/{sample}/02_variants/02_{chromosomes}.process_pileup.snakemake.log`

**Example command executed internally:**

```
python 02_process-pileup.py \
    --sample_folder SampleA \
    --chrs Segment1 \
    --outdir SampleA/02_variants \
    --min_ratio 0.9 \
    --primer_set primers.fa \
    --vcf SampleA/01_assembly/clair3/merge_output.vcf.gz
```

------

## 📤 Outputs

For each sample and chromosome:

- **Pileup files:** `{sample}/02_variants/{chromosomes}_pileup`
- **Variable bases TSV:** `{sample}/02_variants/{chromosomes}_variable_bases.tsv`

These outputs can be used for downstream analysis such as:

- Variant frequency analysis
- Primer validation and trimming assessment
- Generating consensus modifications

------

## 🔧 Notes

- Ensure that the assembly workflow is complete before running this variant workflow.
- The custom Python script assumes standard pileup formatting and Clair3 VCF structure.
- Logs are maintained for troubleshooting and reproducibility.
