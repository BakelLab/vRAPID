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


rule run_mpileup:
    message: "Run samtool's mpileup"
    input:
        fasta  = "{sample}/01_assembly/{sample}.fasta",
        bam    = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
    output:
        pileup = "{sample}/02_variants/{chromosomes}_pileup"
    log:
        "logs/{sample}/02_variants/01_{chromosomes}_mpileup.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	samtools mpileup -f "{config[reference_genome]}" -r "{wildcards.chromosomes}" "{input.bam}" -o "{output.pileup}"
        ) &> "{log}"
        """    


rule process_pileup:
    message: "Process pileup output to get variable bases TSV"
    input:
        pileup        = "{sample}/02_variants/{chromosomes}_pileup",
        vcf           = "{sample}/01_assembly/clair3/merge_output.vcf.gz"
    output:
        tsv           = "{sample}/02_variants/{chromosomes}_variable_bases.tsv"
    params:
        sample_folder = "{sample}",
        outdir        = "{sample}/02_variants",
        chrom         = "{chromosomes}",
        min_ratio     = 0.9,
        pileup_script = os.path.join(workflow.basedir, "../workflow/scripts/02_process-pileup.py")
    log:
        "logs/{sample}/02_variants/02_{chromosomes}.process_pileup.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	python "{params.pileup_script}" \
            	--sample_folder "{params.sample_folder}" \
            	--chrs "{params.chrom}" \
            	--outdir "{params.outdir}" \
            	--min_ratio "{params.min_ratio}" \
            	--primer_set "{config[primer_file]}" \
            	--vcf "{input.vcf}"
        ) &> "{log}"
        """
