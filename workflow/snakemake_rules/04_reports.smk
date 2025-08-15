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

#########
# RULES #
#########


rule run_QualiMap_sample:
    message: "Run Qualimap for each sample" 
    input:
        bam           = expand("{sample}/01_assembly/{sample}_ref.sorted.rg.bam", 
                            sample = sampleids, chromosomes = chromosomes),
        qc_file       = expand("{sample}/03_qualityControl/{sample}_qualityControl.pdf",
                            sample = sampleids, chromosomes = chromosomes),
        flagstat      = expand("{sample}/03_qualityControl/{sample}_refbam.flagstat",
                            sample = sampleids, chromosomes = chromosomes)
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
      
      
rule push_data_pathogendb:
    message: "Push genome assembly data to pathogenDB"
    input:
        qualimap_dirs  = "multi_bamqc/multisampleBamQcReport.html",
        mqc_file       = "multiqc_report.html",
        bam            = "{sample}/03_qualityControl/{sample}_refbam.flagstat"
    output:
        up_log         = "{sample}/04_status/{sample}.assembly-push.log",
    params:
        sample_name    = "{sample}",
        config         = "config.yaml",
    log:
        "logs/{sample}/04_PathogenDB-Push.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/04_all-virus-assembly-push.py"
      
 
rule run_report:
    message: "Generate run summary report"
    input:
        expand("{sample}/04_status/{sample}.assembly-push.log", sample = sampleids)
    output:
        expand("{runid}_run_report.csv", runid = config["run_id"])
    conda:
    	"../envs/env.yml"
    log:
    	"logs/04_run-report.snakemake.log"
    script:
    	"../scripts/04_generate-run-report.R"


rule summary:
    message: "Get a summary report for the packages used in the run"
    input:
        expand("{runid}_run_report.csv", runid = config["run_id"])
    output:
        "workflow_summary.txt"
    log:
        "logs/04_workflow-summary.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/04_generate-workflow-summary.py"

    
