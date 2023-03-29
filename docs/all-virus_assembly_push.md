# All Virus Assembly Push

## Background

This script was built to replace the virus specific assembly push scripts previously used for Influenza, SARS-CoV-2/sCoV/MPX. It was also set up for integration into `snakemake` pipelines. The script is based on a standardized dictionary of values (see below) which is parsed through a series of functions. The standardization of the dictionary allows for more flexibility in handling of viruses with different genome structures (such as segmented vs non-segmented). Database insertions/updates are also handled by function and database commits will only occur if no exceptions are raised. Failed/broken samples run in the script should not cause premature exit, instead any issues will be reported in the output log file and the database will rollback rather than update. 

### Purpose

This script is intended to serve as the final individual step in the viral reference-based assembly pipeline for the Bakel Lab. It will parse assembly information, assign scores based on coverage/completeness, upload the information into PathogenDB, and move report files for access on the frontend pathogendb.mssm.edu site. 


## Script Usage
### Execution

`all-virus_assembly_push.py -s <sample> -r <run_id> -c <config>`

#### Snakemake

The script can be found in the `Snakefile` under the `push_data_pathogendb` rule. All necessary input, parameters, output, and log files are defined under that rule. Note that the snakemake rule also outputs a log file that is separate from the log file that is outputted in within the script. Snakemake's log files can always be found in the `log` directory, and identified per sample in the following format:

`05_sample_name.PDBPush.snakemake.log`

#### Config File Structure

````yaml
{
 'samples': 'samples.csv',
 'project': 'PVI',
 'reference': 'sars-cov-2/COVID.fa',
 'genbank': 'sars-cov-2/COVID.gbk',
 'reverse_primer': 'sars-cov-2/SARS-CoV-2_primers_3prime_NI.fa',
 'forward_primer': 'sars-cov-2/SARS-CoV-2_primers_5prime_NI.fa',
 'primer_set_2kb': 'sars-cov-2/SARS-CoV-2_primers_2kb_set.tsv',
 'primer_set_1_5kb': 'sars-cov-2/SARS-CoV-2_primers_1.5kb_set.tsv',
 'virus': 'SARS-CoV-2',
 'path': '/sc/arion/projects/PVI/{pipeline-dir}/{run-id}/{virus}/{collaborator}',
 'length': '29903',
 'ref_fasta_headers': 'MN908947'
 }
````

##### *Note*

*For viruses that are not multisegmented, the `ref_fasta_headers` and `length` must not be in list format.*

#### Directory Structure/File Requirements

Primary/whatever term

`/sc/arion/projects/PVI/{pipeline}/{run-id}/{virus}/{collaborator}/`

##### Sample specific

`{sample}/02_assembly ` containing 

*  `{sample}.fasta`
* `{sample}_subtype.txt`
  * influenza/sCoV only
* `{sample}.features_cds.fa`
  * influenza only
* `{sample}.features_protein.fa`
  * influenza only
* `{sample}.features_table.txt`
  * influenza only

`{sample}/03_qualityControl` containing

* `{sample_refbam.flagstat}`
* `{sample_kraken_report.out}`
* `*_coverage.txt`
* `{sample}_quality_control.pdf`

`{sample}/04_variants` containing 

* `*_variable_bases.tsv`

### Output

Script will insert new records or update the following pathogenDB tables:

* `tCEIRS_assemblies`
* `tCEIRS_assembly_sequences`
* `tCEIRS_Submissions`

In addition a run log, `{sample}.assembly-push.log`, will be generated in the `{sample}/05_status` directory.




## Troubleshooting

### *Note*

*The script outputs a log file per sample which is structured to contain success or failure messages for each vital step. Failure messages will contain full error codes for what went wrong. Please consult the `{sample}.assembly-push.log` as your first step in any troubleshooting.*

### Dependencies 

Snakemake environment set up will cover all required dependencies. 

### Common Issues

* Missing file
  * the script relies on specific file names/paths. If an error reports a missing file check that the file is named/located in the correct directory for the sample
* KeyError
  * the script runs off of a deeply nested dictionary of values. KeyErrors are possible and the first step should be to ensure that the config file and arguments contain all requires information. 



## Structure

Script structure is described below to allow for ease in future implementation of other virus types, parsing functions, or database values.

### Dictionary Logic

#### Data Dictionary

For each sample a dictionary is built to contain all input and parsed scores. This dictionary is the basis for all data handling and will be the source of the final database uploads. At the start if the script the dictionary will be largely empty and will be populated based on the parsed results from the data handling functions. The final dictionary will follow the below structure.

```python
data_dict = {
     'run': 'TestRun2',
     'sample': 'FK_49196',
     'extract_id': '49196',
     'assembly_id': 26888,
     'virus': 'IAV',
     'length': 1000,
     'subtype': 'H15N2',
     
     'file_paths': {
          'base_dir': 'FK_49196',
          'assembly_dir': 'FK_49196/FK_49196/02_assembly',
          'qc_dir': 'FK_49196/FK_49196/03_qualityControl',
          'variant_dir': 'FK_49196/FK_49196/04_variants',
          'fasta_path': 'FK_49196/FK_49196/02_assembly/FK_49196.fasta',
          'flagstat_path': 'FK_49196/FK_49196/03_qualityControl/FK_49196_refbam.flagstat',
          'kraken_path': 'FK_49196/FK_49196/03_qualityControl/FK_49196_kraken_report.out',
          'cds_path': 'FK_49196/FK_49196/02_assembly/FK_49196.features_cds.fa',
          'protein_path': 'FK_49196/FK_49196/02_assembly/FK_49196.features_protein.fa',
          'features_path': 'FK_49196/FK_49196/02_assembly/FK_49196.features_table.txt',
          'subtype_path': 'FK_49196/FK_49196/02_assembly/FK_49196_subtype.txt'},
  
 	'flagstat': {
        'total_reads': 302918, 
        'mapped_reads': 277062},
 	
 	'sequences': {
 		'primary': {
            'HA': {
            'fasta': 'CTATTAACCATGAAGACTATCATTGCTTT...',
            'type': 'Segment',
            'name': 'HA'}},
  		'CDS': {
            'M1': {
                'fasta': 'ATGAGTCTTCTAACCGA...',
                'type': 'CDS',
                'name': 'M1'},
           'M2': {
                'fasta': 'ATGAGTCTTCTAACCGAG...',
                'type': 'CDS',
                'name': 'M2'},
           'NA': {
               'fasta': 'TGAATCCAAATCAAAAG...',
                'type': 'CDS',
                'name': 'NA'}},
  		'Protein': {
            'M1': {
                'fasta': 'MSLLTEVETYVLSIIP...',
                'type': 'Protein',
                'name': 'M1'},
            'M2': {
                'fasta': 'MSLLTEVE...',
                'type': 'Protein',
                'name': 'M2'},
           'NA': {
               'fasta': '*IQIKR**RLA...',
                'type': 'Protein',
                'name': 'NA'}}},
                
 	'percents': {
      	'eukaryote_per': 0.26,
      	'bacteria_per': 2.31,
     	'viral_per': 56.46,
     	'fungi_per': 0.0,
    	'archaea_per': 0.0,
     	'short_unmapped_per': 8.54,
     	'mapped_per': 91.46},
      
 	'scores': {
 		'completeness': 0.99,
      	'quality': 'Check',
      	'status': 'Partial',
      	'variable_bases': 13,
      	'percent': 56.46,
  		'coverage': {
            'coverage_100': 0.99, 
            'coverage_10': 0.99}}}
}
```

#### Scoring Dictionary

The assignment of scored terms such as "Complete", "Partial", "Failed", and "Passed" is based on set virus specific scoring thresholds. Those thresholds are kept in the `scoring_dict` to allow for easier adjusting of the scoring parameters. For Influenza the segments are weighted based on being complete or nearly complete. This allows us to assign a partial score for samples that have complete HA/NA even if the other segments are missing, or to give a partial score to samples that have either HA/NA complete and all other segments other than either HA/NA complete. 

```python
scoring_dict = {
        'IAV': {
            'segment_weights': {
                'HA': {
                    'complete': 10, 
                    'nearly complete': 5},
   				'NA': {
   					'complete': 10, 
   					'nearly complete': 5},
   				'NP': {
   					'complete': 2, 
   					'nearly complete': 1},
   				etc},
  			'variable_bases': 9,
  			'percent': 0.4,
  			'completeness': 20,
  			'coverage_10': 0.9},
 		'IBV': {
 			'segment_weights': {
 				'HA': {
 					'complete': 10, 
 					'nearly complete': 5},
   				'NA': {
   					'complete': 10, 
   					'nearly complete': 5},
   				'NP': {
   					'complete': 2, 
   					'nearly complete': 1},
   				etc},
           'variable_bases': 9,
           'percent': 0.4,
           'completeness': 20,
           'coverage_10': 0.7},
 		'SARS-CoV-2': {
 			'variable_bases': 9,
            'percent': 0.4,
            'completeness': 0.7,
            'coverage_10': 0.7},
 		'sCoV': {
            'variable_bases': 9,
            'percent': 0.4,
            'completeness': 0.7,
            'coverage_10': 0.7},
         'MPX': {
              'variable_bases': 9,
              'percent': 0.4,
              'completeness': 0.7,
              'coverage_10': 0.7}
}
```

### Functions

There are 3 categories of functions; parsing functions, file handling functions, and database functions. The script should remain structured so that functions take values from the `data_dict` and/or `scoring_dict` as input, do work, and then optionally update or add values to the `data_dict`. This will keep the script capable of handling diverse data inputs by forcing standardization onto them.

* Parsing functions are those that process the various flagstat, kraken, fasta files to assign values to the dictionary, including the scoring functions. 

* File handling functions deal with the files structuring or movement for the final reports for links on pathogendb.mssm.edu.

* Database functions handle all database interactions, including checking if records exist, insertions into the database, updates of the database, and deletion of previous records in the case of the `tCEIRS_Submissions` table record handling. 

---
Author: Adriana van de Guchte
Contact: vandeguchtea@gmail.com