#!/usr/bin/env snakemake

##########
# IMPORT #
##########

import pandas as pd
from snakemake.utils import min_version
min_version("9.0.0")

##########
# CONFIG #
##########

configfile: "config.yaml"
mappings  = pd.read_csv(config["samples"])
sampleids = mappings["Sample_ID"].tolist()

########
# DEF #
########

reference   = config['reference_genome'] 
chromosomes = config['ref_fasta_headers']

# Database connection information defaults
PDB_MY_CNF_FILE             = config.get("pdb_my_cnf_file",       "~/.my.cnf")
PDB_MY_CNF_GROUP            = config.get("pdb_my_cnf_group",      "vanbah01_pathogens_root")

#########
# RULES #
#########
      
      
rule push_data_pathogendb:
    message: "Push genome assembly data to pathogenDB"
    input:
        bam            = "{sample}/03_qualityControl/{sample}_refbam.flagstat"
    output:
        up_log         = "{sample}/04_status/{sample}.assembly-push.log",
    params:
        sample_name    = "{sample}",
        config         = "config.yaml",
    log:
        "logs/{sample}/04_PathogenDB-Push.snakemake.log"
    conda:
        "../envs/PDB_connect.yaml"
    script:
        "../scripts/04_all-virus-assembly-push.py"
      
 
rule run_report:
    message: "Generate run summary report"
    input:
        expand("{sample}/04_status/{sample}.assembly-push.log", sample = sampleids)
    output:
        expand("{runid}_run_report.csv", runid = config["run_id"])
    params:
        pdb_db_group = PDB_MY_CNF_GROUP,
        pdb_db_config = PDB_MY_CNF_FILE
    conda:
        "../envs/PDBload.yaml"
    log:
        "logs/04_run-report.snakemake.log"
    script:
        "../scripts/04_generate-run-report.R"


rule run_QualiMap_sample:
    message: "Run Qualimap for each sample" 
    input:
        bam           = expand("{sample}/01_assembly/{sample}_ref.sorted.rg.bam", 
                            sample = sampleids, chromosomes = chromosomes),
        up_log        = expand("{sample}/04_status/{sample}.assembly-push.log",
                            sample = sampleids),
        run_report    = expand("{runid}_run_report.csv",
                            runid = config["run_id"])
    output:
        qualimap_dirs = "multi_bamqc/multisampleBamQcReport.html"
    log:
        "logs/04_qualimap.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
      r"""
      (
          qualimap multi-bamqc -d "{config[qualimap_input]}" -r
      ) &> "{log}"
      """
        
      
rule run_multiqc:
    message: "Generate MultiQC Report for the run"
    input:
        qc_report = "multi_bamqc/multisampleBamQcReport.html"
    output:
        mqc_file  = "multiqc_report.html"
    log:
        "logs/04_multiqc.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
      """
      (
          multiqc .
      ) &> "{log}"
      """


rule summary:
    message: "Get a summary report for the packages used in the run"
    input:
        mqc = "multiqc_report.html",
        csv = expand("{runid}_run_report.csv", runid = config["run_id"])
    output:
        "workflow_summary.txt"
    log:
        "logs/04_workflow-summary.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/04_generate-workflow-summary.py"

    
