### vRAPID

A reference-based pipeline for the assembly, consensus calling and annotation of viral pathogens.

vRAPID is dependent on folder structures and naming schemes used by the genomics sequencing core at the Icahn School of Medicine at Mount Sinai.

vRAPID is an expansion of the [previous pipeline](https://github.com/mjsull/COVID_pipe) used by the Bakel Lab that was written by Mitchell Sullivan. 

#### Installing snakemake on Minerva

Install `snakemake` on minerva before running the pipeline, as follows:

```bash
module purge all
unset PYTHONPATH
unset PERL5LIB
unset R_LIBS
module load anaconda3
conda config --set channel_priority strict
conda create -c conda-forge -c bioconda -n snakemake snakemake=6.8.0 mamba=0.24 tabulate=0.8
```

#### Setting up a vRAPID run on Minerva

The following folder structure should exist in the run directory

```bash
<Run-ID>
└───<sample_ID1>
│   │   <sample_ID1>_1.fastq.gz
│   │   <sample_ID1>_2.fastq.gz
│
└───<sample_ID2>
    │   <sample_ID2>_1.fastq.gz
    │   <sample_ID2>_2.fastq.gz
```

Then run:

`python ~/opt/vRAPID/workflow/scripts/organize_run_samples.py`

This script queries to PathogenDB, separating the samples into subdirectories based on the virus type, expected subtype, and collaborator ID. If no subtype is expected, then the script separates the samples into a `None` subdirectory. This creates the following directory structure within the run directory:

```bash
<Run-ID>
├── <virus1>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
├── <virus2>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz

```

Then within the respective newly created subdirectory, create the following two files:

1. `samples.csv:` with the column name `Sample_ID`, holding a list of the samples within the newly created subdirectory (the new working directory) that are to be run for that virus type and collaborator.
2. `multibamqc_input.txt:` a tab-separated file that holds the list of samples from `samples.csv` in the first column, and the path to the bam file in the second column. Note that the bam files will always be outputted to `<sample_ID>/02_assembly/<sample_ID>_ref.bam`. Please note that the `multibamqc_input.txt` file **does not** require column names.

An example of both file structures can be seen below:

`samples.csv`

```bash
Sample_ID
sample_ID1
sample_ID2
sample_ID3
```

`multibamqc_input.txt`

```bash
sample_ID1        sample_ID1/02_assembly/sample_ID1_ref.bam
sample_ID2        sample_ID2/02_assembly/sample_ID2_ref.bam
sample_ID3        sample_ID3/02_assembly/sample_ID3_ref.bam
```

The final structure should be as follows:

```bash
<Run-ID>
├── <virus1>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
├── <virus2>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│		├── <sample_ID1>
│   │	│	│   │		├── <sample_ID1>_1.fastq.gz
│   │	│	│   │		├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│		├── <sample_ID2>
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
│   │	│	│   │		├── <sample_ID2>_1.fastq.gz
```

Lastly, to successfully run the pipeline, depending on the virus, the `config.yaml` file within the pipeline repository might need to be edited. The following fields are required within the file:

1. **samples**: Name of the file containing a list of the samples being run. Currently set to `samples.csv`.
2. **run_id:** Run ID name, typically identified as `TD######` from the sequencing core's data release e-mail.
3. **reference:** path to the reference within the `db` directory that is found in the repository. Currently set to `sars-cov-2/COVID.fa`.
4. **genbank:** path to the GenBank file within the `db` directory that is found in the repository. Currently set to `sars-cov-2/COVID.gbk`.
5. **reverse_primer:** path to the 3-prime primer file within the `db` directory that is found in the repository. Currently set to `sars-cov-2/SARS-CoV-2_primers_3prime_NI.fa`.
6. **forward_primer:** path to the 5-prime primer file within the `db` directory that is found in the repository. Currently set to `sars-cov-2/SARS-CoV-2_primers_5prime_NI.fa`.
7. **primer_set_2kb:**  path to the 2kb-prime set file within the `db` directory that is found in the repository. Currently set to `sars-cov-2/SARS-CoV-2_primers_2kb_set.tsv`.
8. **primer_set_1_5kb:**  path to the 1.5kb-prime set file within the `db` directory that is found in the repository. Currently set to `sars-cov-2/SARS-CoV-2_primers_1.5kb_set.tsv`.
9. **virus:** the virus name for the samples being run. Currently set to `SARS-CoV-2`. Other options include: `Influenza-A`, `Influenza-B`, `sCoV`, `MPX`.
10. **path:** the full path to where the samples being run are located. See above for the proper structure.
11. **length:** length of the reference genome(s). For multi-segmented viruses like influenza, this can be a list of lengths, ordered with respect to the FASTA headers.
12. **ref_fasta_headers:** The name of FASTA header(s) in the reference genome. For multi-segmented viruses like influenza, this can be a list of headers, ordered with respect to the length.

#### Running vRAPID

Then once the files are generated, the pipeline can be run using the following command:

`submitjob 12 -c 4 -m 64 -q private -P acc_PVI ~/opt/vRAPID/run-vRAPID-analysis -i <run-ID> -p <path>`

Note that both the `<run-ID> ` and the `<path>` arguments in the command above are optional. If they are not supplied, then they are pulled from the `config.yaml` file.

#### vRAPID output

Once the pipeline has completed, the following output files are created:

```bash
<Run-ID>
├── <virus1>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── pipeline.log
│   │	│	├── file_movement_message.txt
│   │	│	├── <Run-ID>_cleaned.txt
│   │	│	├── <Run-ID>_data_release.txt
│   │	│	├── <Run-ID>_software_version.txt
│   │	│	├── <Run-ID>_<collaborator-ID>_assemblies.csv
│   │	│	├── <Run-ID>_<collaborator-ID>.tar.gz
│   │	│	├── logs
│   │	│	├── multi_bamqc
│   │	│	├── multiqc_report.html
│   │	│	├── <sample_ID>
│   │	│	│	├── 01_fastqs
│   │	│	│	│	├── <sample_ID>_1.fastq.gz
│   │	│	│	│	├── <sample_ID>_2.fastq.gz
│   │	│	│	├── 02_assembly
│   │	│	│	├── 03_qualityControl
│   │	│	│	├── 04_variants
│   │	│	│	├── 05_status
```

