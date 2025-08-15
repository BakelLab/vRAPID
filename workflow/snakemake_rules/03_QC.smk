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


rule compute_coverage:
    message: "Use ref BAM to get coverage per FASTA header"
    input:
        bam      = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
        txt      = "{sample}/01_assembly/{sample}_picard_output.txt"
    output:
        coverage = "{sample}/03_qualityControl/{chromosomes}_coverage.txt"
    params:
        chrs     = "{chromosomes}"
    log:
        "logs/{sample}/03_qualityControl/01_{chromosomes}.coverage.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	samtools depth -aa "{input.bam}" -r "{params.chrs}" > "{output.coverage}"
        ) &> "{log}"
        """


rule generate_qc_report:
    message: "Plot coverage per FASTA header"
    input:
        coverage    = "{sample}/03_qualityControl/{chromosomes}_coverage.txt",
        fasta       = "{sample}/01_assembly/{chromosomes}_polished.fasta",
    output:
        pdf         = temp("{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control2.pdf"),
        report      = "{sample}/03_qualityControl/{chromosomes}_report.txt"
    params:
        qc_dir      = "{sample}/03_qualityControl/",
        sample_name = "{sample}",
        chromosome  = "{chromosomes}"
    log:
        "logs/{sample}/03_qualityControl/02_{chromosomes}.initial-QC-report.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/03_generate-QC-report.py"
    
    
rule add_variant_plot_to_pdf:
    message: "Plot variants from reference using pileup output and append to coverage plot file"
    input:
        pileup      = "{sample}/02_variants/{chromosomes}_pileup",
        pdf         = "{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control2.pdf",
    output:
        pdf         = temp("{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control3.pdf")
    params:
        sample_name = "{sample}",
        chromosomes = "{chromosomes}"
    log:
        "logs/{sample}/03_qualityControl/03_{chromosomes}.plot-variants.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/03_plot-variants.py"  


rule add_primers_plot_to_pdf:
    message: "Plot primer depth and add to the QC report"
    input:
        pdf           = "{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control3.pdf"
    output:
        pdf           = temp("{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control4.pdf"),
    params:
        sample_folder = "{sample}",
        qc_dir        = "{sample}/03_qualityControl",
        chromosomes   = "{chromosomes}"
    log:
        "logs/{sample}/03_qualityControl/04_{chromosomes}.plot-primers.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/03_plot-primer-set-depth.py"
    
    
rule plot_intrahosts:
    message: "Plot intrahosts for each FASTA header"
    input:
        pdf = "{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control4.pdf",
        tsv = "{sample}/02_variants/{chromosomes}_variable_bases.tsv"
    output:
        pdf = temp("{sample}/03_qualityControl/{chromosomes}_var.pdf"),
        tsv = "{sample}/03_qualityControl/{chromosomes}_variants_labels.tsv"
    params:
        intrahost_plot = os.path.join(workflow.basedir, "../workflow/scripts/03_plot-intrahosts.R")
    log:
        "logs/{sample}/03_qualityControl/05_{chromosomes}.plot-intrahost.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
    	r"""
    	(
    		Rscript "{params.intrahost_plot}" -i "{input.tsv}" -o "{output.pdf}" -v "{output.tsv}"
    	) &> "{log}"
    	"""

rule run_kraken:
    message: "Run kraken2 to get the read breakdown"
    input:
        R1     = "{sample}/00_fastqs/{sample}_1.fastq.gz",
        R2     = "{sample}/00_fastqs/{sample}_2.fastq.gz",
    output:
        report = "{sample}/03_qualityControl/{sample}_kraken_report.out",
        out    = "{sample}/03_qualityControl/{sample}_kraken"
    log:
        "logs/{sample}/03_qualityControl/06_kraken.snakemake.log"
    conda:
        "../envs/env.yml"
    threads: 12
    shell:
        r"""
        (
        	kraken2 \
        		--db "{config[kraken_db]}" \
        		--quick \
        		--report "{output.report}" \
        		--threads "{threads}" \
        		--output "{output.out}" \
        		--paired "{input.R1}" "{input.R2}"
        ) &> "{log}"
        """


rule plot_taxonomic_breakdown:
    message: "Plot taxanomic breakdown using kraken2 output"
    input:
        report = "{sample}/03_qualityControl/{sample}_kraken_report.out",
    output:
        pdf   = temp("{sample}/03_qualityControl/taxonomic-breakdown.pdf")
    log:
        "logs/{sample}/03_qualityControl/07_taxonomic-breakdown.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/03_plot-taxonomic-breakdown.R"


rule merge_pdfs:
    message: "Merge QC plots with kraken plot and the intrahost plot(s)"
    input:
        pdf_1 = lambda wildcards: expand("{sample}/03_qualityControl/{sample}-{chromosomes}_quality_control4.pdf",
                                         sample = [wildcards.sample], chromosomes = chromosomes),
        pdf_2 = lambda wildcards: f"{wildcards.sample}/03_qualityControl/taxonomic-breakdown.pdf",
        pdf_3 = lambda wildcards: expand("{sample}/03_qualityControl/{chromosomes}_var.pdf",
                                         sample = [wildcards.sample], chromosomes = chromosomes),
    output:
        pdf   = "{sample}/03_qualityControl/{sample}_qualityControl.pdf"
    log:
        "logs/{sample}/03_qualityControl/08_merge-pdfs.snakemake.log"
    conda:
        "../envs/env.yml"
    script:
        "../scripts/03_merge-pdfs.py"
      
        
rule run_flagstat:
    message: "Run samtool's flagstat"
    input:
        bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
        pdf = "{sample}/03_qualityControl/{sample}_qualityControl.pdf"
    output:
        bam = "{sample}/03_qualityControl/{sample}_refbam.flagstat"
    log:
        "logs/{sample}/03_qualityControl/09_samtools-flagstat.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	samtools flagstat "{input.bam}" > "{output.bam}"
        ) &> "{log}"
        """

