### vRAPID

A reference genome-based pipeline for the assembly, consensus calling and annotation of viral pathogen genomes.

vRAPID relies on data structures and naming conventions used by the [Center for Research on Influenza Pathogenesis and Transmission (CRIPT)](https://www.ceirr-network.org/centers/cript) and the [Mount Sinai Pathogen Surveillance Program (MS-PSP)](https://icahn.mssm.edu/research/genomics/research/microbial-genomics-pathogen-surveillance) at the Icahn School of Medicine at Mount Sinai. vRAPID is an expansion of the [COVID_pipe](https://github.com/mjsull/COVID_pipe) pipeline that was written by Mitchell Sullivan in the Bakel Lab. 




#### Prepare the snakemake conda environment

Installation of the required external software packages is largely handled by the pipeline itself, however a conda environment named `snakemake` needs to be present in your environment. We recommend miniconda, which is a free minimal installer for [conda](https://docs.conda.io/en/latest/miniconda.html). Follow the instructions below to start the miniconda installer on Linux. When asked whether the conda environment should automatically be initialized, select 'yes'. Note that Snakemake requires the channel_priority to be set to strict. The post-installation commands to apply this setting are included in the post-installation selection below.

```
# Start miniconda installation
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh

# Post-installation commands to enforce strict channel_priority (required for Snakemake)
conda config --set auto_activate_base false
conda config --set channel_priority strict
```

On the Mount Sinai 'Minerva' cluster we additionally recommend setting the conda environment `pkgs` and `envs` directories to your `work` folder, which has a much larger storage allocation (100-200 Gb) than the home folders (5 Gb). The conda configuration commands to set the directories are as follows:

```
mkdir -p /sc/arion/work/$USER/conda/pkgs
mkdir -p /sc/arion/work/$USER/conda/envs
conda config --add pkgs_dirs "/sc/arion/work/$USER/conda/pkgs"
conda config --add envs_dirs "/sc/arion/work/$USER/conda/envs"
```

After installation and configuration of conda/miniconda, the following 'conda create' command can be used to set up the required 'snakemake' environment.

```
conda create -c conda-forge -c bioconda -n snakemake 'snakemake=6.8.0' 'mamba=0.24' 'tabulate=0.8'
```



#### Running the vRAPID pipeline

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
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── <sample_ID1>
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
├── <virus2>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── <sample_ID1>
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── <sample_ID1>
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz

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
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
├── <virus2>
│   ├── <expect_subtype>
│   │   ├── <collaborator_ID1>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │   ├── <collaborator_ID2>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── <sample_ID1>
│   │	│	│   ├── <sample_ID1>
│   │	│	│   │	├── <sample_ID1>_1.fastq.gz
│   │	│	│   │	├── <sample_ID1>_2.fastq.gz
│   │	│	├── <sample_ID2>
│   │	│	│   ├── <sample_ID2>
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
│   │	│	│   │	├── <sample_ID2>_1.fastq.gz
```

Lastly, to successfully run the pipeline, depending on the virus, the `config.yaml` file within the pipeline repository might need to be edited. Please note that the current default settings are set for SARS-CoV-2. The following fields are required within the file:

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

`submitjob 12 -c 4 -m 64 -q private -P acc_PVI ~/opt/vRAPID/run-vRAPID-analysis -i <run-ID> -p <path> -s <samples-run-csv>`

Note that the `<run-ID> `, `<path>`, and `<samples-run-csv>` arguments in the command above are optional. If they are not supplied, then they are pulled from the `config.yaml` file.

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

#### Accessory vRAPID pipeline scripts

The following is a brief description of the scripts that are called within the `Snakefile`. For more information, please see the  function overview below:

* [run_pipeline.py](run_pipeline.md) : Takes in Illumina paired-end FASTQ files and assembles them against a reference genome, outputting a final `FASTA` consensus genome sequence and aligned read `bam` file, and other intermediary output files.

* [variant_analysis.py](variant_analysis.md): Takes in the consensus FASTA file generated by `run_pipeline.py` and performs an intra-host variant analysis on the virus library. 

* [run_QC.py](run_QC.md): Takes in the `bam` file generated as an output from the `run_pipeline.py` script, and performs a quality control analysis on the virus library.

* [run-genome-annotation.py](run-genome-annotation.md): Takes in the `FASTA` file generated from `run_pipeline.py` and annotates the consensus according to the virus type. For Inluenza samples, this is done using [NCBI's Influenza Annotation Tool](https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/api/help.html). For SARS-CoV-2, this is done using [`vadr`](https://github.com/ncbi/vadr).

* [all-virus_assembly_push.py](all-virus_assembly_push.md): Used to upload viral genome assembly data to the Mount Sinai Pathogen Surveillance Program database (PathogenDB).

* [move_QC_PDB.py](move_QC_PDB.md): Uploads the output that is obtained by running [`multiqc`](https://github.com/ewels/MultiQC), and [`qualimap`](http://qualimap.conesalab.org) to the Mount Sinai Pathogen Surveillance Database (PathogenDB).

* [cleanup.py](cleanup.md): Removes large intermediary output files once the pipeline has completed to conserve space before archival.

* [generate_run_report.R](generate_run_report.md): Creates a final assembly run status and quality report for the samples that were processed within the run, outputting a CSV file.

* [data-release.py](data-release.md): Prepares vRAPID pipeline outputs for data release.
