### vRAPID

The Virus Reference-based Assembly Pipeline and IDentification (vRAPID) tool runs the assembly, consensus calling and annotation of viral pathogen genomes.

vRAPID relies on data structures and naming conventions used by the [Center for Research on Influenza Pathogenesis and Transmission (CRIPT)](https://www.ceirr-network.org/centers/cript) and the [Mount Sinai Pathogen Surveillance Program (MS-PSP)](https://icahn.mssm.edu/research/genomics/research/microbial-genomics-pathogen-surveillance) at the Icahn School of Medicine at Mount Sinai. vRAPID is an expansion of the [COVID_pipe](https://github.com/mjsull/COVID_pipe) pipeline that was written by Mitchell Sullivan in the Bakel Lab. 



## Installation

**Set up the snakemake conda environment**

Installation of the required external software packages is largely handled by the pipeline itself, however a conda environment named `snakemake` needs to be present in your environment. We recommend miniconda, which is a free minimal installer for [conda](https://docs.conda.io/en/latest/miniconda.html). Follow the instructions below to start the miniconda installer on Linux. When asked whether the conda environment should automatically be initialized, select 'yes'. Note that Snakemake requires the channel_priority to be set to strict. The post-installation commands to apply this setting are included in the post-installation selection below.

```bash
# Start miniconda installation
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh

# Post-installation commands to enforce strict channel_priority (required for Snakemake)
conda config --set auto_activate_base false
conda config --set channel_priority strict
```

On the Mount Sinai 'Minerva' cluster we additionally recommend setting the conda environment `pkgs` and `envs`directories to your `work` folder, which has a much larger storage allocation (100-200 Gb) than the home folders (5 Gb). The conda configuration commands to set the directories are as follows:

```bash
mkdir -p /sc/arion/work/$USER/conda/pkgs
mkdir -p /sc/arion/work/$USER/conda/envs
conda config --add pkgs_dirs "/sc/arion/work/$USER/conda/pkgs"
conda config --add envs_dirs "/sc/arion/work/$USER/conda/envs"
```

After installation and configuration of conda/miniconda, the following 'conda create' command can be used to set up the required 'snakemake' environment.

```bash
conda create -c conda-forge -c bioconda -n snakemake 'snakemake>=8.25.0' 'conda>=24.7.1'
```

To make use of snakemake's parallelization of jobs, it is important that you have the LSF plugin installed if you are running the pipeline on Mount Sinai's HPC system, Minerva.

```bash
conda activate snakemake
pip install snakemake-executor-plugin-cluster-generic
conda deactivate snakemake
```

**Clone the vRAPID repository**

Use `git clone` to clone the vRAPID repository, e.g. in your `~/opt/` folder.

```bash
git clone https://github.com/BakelLab/vRAPID.git
```

**Set up the PDB_connect conda repository and .my.cnf configuration file**

Several utility scripts exist that facilitate processing of data for samples sequenced using the Mount Sinai Pathogen Surveillance Program workflow. These scripts connect to the central PathogenDB tracking system to retrieve extract, isolate and assembly information. As such, they require an environment with perl and python database connection libraries. The required libraries can be installed using conda and the provided `PDB_connect.yaml` file in the vRAPID repository.

Assuming that the vRAPID repository was cloned into `~/opt/vRAPID`, the following command will prepare the environment:

```bash
conda env create -f ~/opt/vRAPID/workflow/envs/PDB_connect.yaml
```

In addition to the python and perl libraries, the PathogenDB connection also requires a `~/.my.cnf` configuration file installed in your home directory. This configuration file should contain two entries, one for a read/write account that does not have deletion privileges and another for a root account that does have deletion privileges, as follows:

```
[vanbah01_pathogens_rw]
user=pathogendb_rw
database=vanbah01_pathogens
password=<PASSWORD>
host=data1
port=3306

[vanbah01_pathogens_root]
user=vanbah01_root
database=vanbah01_pathogens
password=<PASSWORD>
host=data1.hpc.mssm.edu
port=3306
```

Make sure to change the `<PASSSWORD>` entries to the actual passwords associated with these accounts.

## Prepare sample folders with sequence data for vRAPID assembly

vRAPID processes data organized into individual run folders for each sample. Each run folder must contain the paired-end sequence read files, and the file names must share the same prefix as the run folder name.

For better organization, it is recommended to group samples by sequencing run and virus type using top-level sub-folders for each virus if the run has multiple viruses. For example, **SARS-CoV-2** (for COVID-19 samples) and **Influenza A** (for IAV). These top-level folders must also include a `config.yaml`file that specifies the configuration parameters for the respective run and virus type.

To create the `config.yaml` file, use the template provided in the GitHub repository at `~/opt/vRAPID/config/config.yaml` and adjust it as needed for your analysis.

Note that if you are looking to run the pipeline in parallel, you will need to adjust the workflow-specific `config.yaml` found in the GitHub repository at `~/opt/vRAPID/config/workflow-config.yaml` to your analysis, ensuring the updated paths are also reflected in the workflow wrapper, `~/opt/vRAPID/run-vRAPID`. This is to be placed in `/hpc/users/USER/work/snakemake/profiles/lsf/vRAPID/config.yaml` (or adjust the file path as needed).

Below is an example of what the run folder structure should look like for a given sequencing run:

```
SEQUENCING-RUN-FOLDER
 |-
 |- SARS-CoV-2/
 |   |- VS/
 |   | 	|- config.yaml
 |   | 	|- VS_12341/
 |   | 	|		|- VS_12341/
 |   |   			|- VS_12341_1.fastq.gz
 |   |   			|- VS_12341_2.fastq.gz
 |   | 	|- VS_12342/
 |   | 	|		|- VS_12342/
 |   |   			|- VS_12342_1.fastq.gz
 |   |   			|- VS_12342_2.fastq.gz
```

**Automated Run Folder Organization for Mount Sinai Pathogen Surveillance Program Samples**

For samples sequenced using the Mount Sinai Pathogen Surveillance Program workflow and stored in the PathogenDB database, you can use the provided utility script `prepare-psp-flugap-run` to automatically generate the required run folder structure.

This workflow includes an additional level of organization by sample submitter, with folders grouped by the submitter's two-letter prefix.

The `prepare-psp-flugap-run` script requires the path to the core run folder and a run ID as inputs. The output folder for vRAPID is set to /sc/arion/projects/PVI/genomes/SARS_CoV_2/PSP/2025_2026/ by default but can be changed with the -o argument. For example:

```bash
# Activate the PathogenDB connection environment and run the preparation script
conda activate PDB_connect
~/opt/vRAPID/prepare-psp-vRAPID-run \
   -i /PATH/TO/CORE/FOLDER \
   -r RUNID \
   -o /sc/arion/projects/PVI/genomes/SARS_CoV_2/PSP/2025_2026/
conda deactivate
```

Note that the run script places a config.yaml file in each toplevel folder that contains the required run-specific parameters (virus type and run ID). Make sure to check these and other parameters before running the pipeline.

Also note that this creates the supplementary files needed to run the pipeline such as the `clusterlogs` directory that holds the log files for each submitted job when the pipeline is run, the `samples.txt` that holds the list of samples to run, and `multibamqc_input.txt` that holds the mapping of each samples and the path to the bam file used when running multibamqc.

For example:

***samples.txt***

```bash
Sample_ID
VS_12341
VS_12342
```

***multibamqc_input.txt***

```bash
VS_12341	VS_12341/02_assembly/VS_12341_ref.bam
VS_12342	VS_12342/02_assembly/VS_12342_ref.bam
```

## Running the vRAPID pipeline

For each top-level run folder containing sample subfolders, execute the following loop to start the snakemake pipeline for each sample:

```bash
# NOTE: make sure that you do not have any conda environments loaded when executing this command
submitjob 12 -c 1 -m 5 -q premium -P acc_PVI ~/opt/vRAPID/run-flugap -i ${RUNID} -p ${PWD} -s samples.txt
```

Note that the config needs to exist in the run directory. See the config directory for the `config.yaml` format for the pipeline. This submits each pipeline step as a seperate job.

#### vRAPID output

Once the pipeline has completed, the following output files are created:

```bash
<Run-ID>
├── <virus1>
│   ├── <expect_subtype>
│   │   ├── <batch_ID1>
│   │	│	├── samples.csv
│   │	│	├── multibamqc_input.txt
│   │	│	├── pipeline.log
│   │	│	├── file_movement_message.txt
│   │	│	├── <Run-ID>_cleaned.txt
│   │	│	├── <Run-ID>_data_release.txt
│   │	│	├── <Run-ID>_software_version.txt
│   │	│	├── <Run-ID>_<batch-ID>_assemblies.csv
│   │	│	├── <Run-ID>_<batch-ID>.tar.gz
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
