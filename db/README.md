## 🧬 vRAPID Pipeline: Reference Data Directory (`db`)

### Overview

This directory houses critical reference files and databases required for various stages of the `vRAPID` pipeline, including genome assembly, variant analysis, and quality control. These resources are referenced within the Snakefile rules to ensure accurate and efficient processing of viral genomic data.

### Key Components

- **Reference Genomes**: Complete viral genome sequences used as references for alignment and assembly.
- **Primer Sets**: Files specifying primer sequences for targeted amplification regions.
- **GenBank Files**: Annotations and metadata associated with reference genomes, facilitating functional annotation and gene prediction.
- **Database Files**: Supporting files for tools like Kraken2, used for taxonomic classification and quality control assessments.

### Usage

The resources in this directory are utilized by various rules within the pipeline, such as:

- **Assembly**: Reference genomes and primer sets guide the assembly process.
- **Variant Analysis**: Reference genomes and primer sets are essential for variant detection and analysis.
- **Quality Control**: GenBank files and database files support quality control measures and taxonomic classification.

### Configuration

Paths to these resources are specified in the `config.yaml` file, ensuring flexibility and ease of updates. For example, entries like `reference`, `forward_primer`, and `reverse_primer` point to the respective files in this directory.

------

## 📁 Directory Structure

The `snakefile/db` directory is organized as follows:

```bash
GITREPO/db/
├── virus/
│   └── viral_reference.fasta
│   └── primer_set_1_5kb.fasta
│   └── primer_set_2kb.fasta
│   └── viral_reference.gbk
```

- **`reference FASTA/`**: Contains the complete viral genome sequences.
- **`Primer Sets`**: Holds the primer sequences for targeted amplification.
- **`GenBank Annotations`**: Includes GenBank files with annotations and metadata.

------

## ⚙️ Integration with the Pipeline

The resources in this directory are integrated into the pipeline through the `config.yaml` file, where paths to these resources are specified. For instance:

```
reference: "db/reference_genomes/viral_reference.fasta"
forward_primer: "db/primer_sets/primer_set_1_5kb.fasta"
reverse_primer: "db/primer_sets/primer_set_2kb.fasta"
genbank: "db/genbank_annotations/viral_reference.gbk"
kraken2_db: "db/kraken2_db/kraken2_standard_db_2024"
```

These configurations ensure that the pipeline can dynamically access the necessary resources during execution.

------

## 📄 License

The contents of this directory are licensed under the same terms as the main `vRAPID` repository. Please refer to the repository's `LICENSE` file for detailed information.

------

For more information on the pipeline's functionality and usage, please refer to the main repository's README and documentation.