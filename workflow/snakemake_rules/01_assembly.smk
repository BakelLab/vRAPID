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


rule run_cutadapt:
    message: "Trim primers using cutadapt"
    input:
        R1 = "{sample}/00_fastqs/{sample}_1.fastq.gz",
        R2 = "{sample}/00_fastqs/{sample}_2.fastq.gz"
    output:
        R1 = "{sample}/01_assembly/reads.1.fastq.gz",
        R2 = "{sample}/01_assembly/reads.2.fastq.gz"
    log:
        "logs/{sample}/01_assembly/01_run-cutadapt.snakemake.log"
    conda: 
        "../envs/env.yml"
    shell:
        r"""
        (
        	cutadapt \
            	-g file:"{config[forward_primer]}" \
            	-a file:"{config[reverse_primer]}" \
            	-G file:"{config[forward_primer]}" \
            	-A file:"{config[reverse_primer]}" \
            	-o "{output.R1}" -p "{output.R2}" \
            	"{input.R1}" "{input.R2}"
        ) &> "{log}"
        """

    
rule run_fastqc:
    message: "Run FASTQC on FASTQs"
    input:
        R1 = "{sample}/01_assembly/reads.1.fastq.gz",
        R2 = "{sample}/01_assembly/reads.2.fastq.gz"
    output:
        R1 = "{sample}/01_assembly/reads.1_fastqc.html",
        R2 = "{sample}/01_assembly/reads.2_fastqc.html"
    log:
        "logs/{sample}/01_assembly/02_run-fastqc.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	fastqc --nogroup "{input.R1}" "{input.R2}" --outdir $(dirname "{output.R1}")
        ) &> "{log}"
        """


rule minimap2_clip:
    message: "Align reads to reference, clip, and generate final BAM"
    input:
        fR2 = "{sample}/01_assembly/reads.2_fastqc.html",
        cR1 = "{sample}/01_assembly/reads.1.fastq.gz",
        cR2 = "{sample}/01_assembly/reads.2.fastq.gz"
    output:
        bam = "{sample}/01_assembly/{sample}_ref.bam",
        bai = "{sample}/01_assembly/{sample}_ref.bam.bai"
    log:
        "logs/{sample}/01_assembly/03_minimap2-clip.snakemake.log"
    conda: 
        "../envs/env.yml"
    shell:
        r"""
        (
        	# Align and sort
        	minimap2 -t "{threads}" -ax sr "{config[reference_genome]}" "{input.cR1}" "{input.cR2}" \
            	| samtools sort -@ "{threads}" -o "{output.bam}"
        
        	samtools index "{output.bam}"

        	# Create clipped BAM (remove soft-clipped reads)
        	samtools view -h "{output.bam}" \
            	| awk '$6 !~ /S/' \
            	| samtools view -b - > "{output.bam}.tmp"

        	mv "{output.bam}.tmp" "{output.bam}"
        	samtools index "{output.bam}"
        ) &> "{log}"
        """


rule sorted_bam:
    message: "Generate fully sorted BAM"
    input:
        bam = "{sample}/01_assembly/{sample}_ref.bam"
    output:
        bam = "{sample}/01_assembly/{sample}_ref.sorted.bam"
    log:
        "logs/{sample}/01_assembly/04_sort-bam.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	samtools sort -o "{output.bam}" "{input.bam}"
        	samtools index "{output.bam}"
        ) &> "{log}"
        """


rule group_reads_in_bam:
    message: "Run samtools to group reads in BAM file"
    input:
        sorted_bam    = "{sample}/01_assembly/{sample}_ref.sorted.bam"
    output:
        sorted_rg_bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
        sorted_rg_bai = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam.bai"
    log:
        "logs/{sample}/01_assembly/05_read-group-bam.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	samtools addreplacerg -r '@RG\tID:1\tSM:sample\tPL:ILLUMINA' -o "{output.sorted_rg_bam}" "{input.sorted_bam}"
        	samtools index "{output.sorted_rg_bam}"
        ) &> "{log}"
        """


rule run_clair3:
    message: "Run Clair3 on final ref BAM"
    input:
        sorted_rg_bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
        model         = "ilmn"
    output:
        vcf           = "{sample}/01_assembly/clair3/merge_output.vcf.gz",
        folder        = directory("{sample}/01_assembly/clair3")
    params:
        samples       = "{sample}"
    log:
        "logs/{sample}/01_assembly/06_run-clair3.snakemake.log"
    container: "docker://hkubal/clair3:v1.2.0"
    shell:
        r"""
        (
        	run_clair3.sh --bam_fn="{input.sorted_rg_bam}" \
        		--ref_fn="{config[reference_genome]}" \
        		--threads="{threads}" \
        		--platform=ilmn \
        		--model_path=./ilmn \
        		--output="{output.folder}" \
        		--include_all_ctgs
        ) &> "{log}"
        """


rule get_masking_sites:
    message: "Get masking sites for low coverage regions"
    input:
        sorted_rg_bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
    output:
       mask_sites     = "{sample}/01_assembly/{chromosomes}_coverage_mask.txt"
    params:
        depth_masking = os.path.join(workflow.basedir, "../workflow/scripts/01_generate-depth-masking.py")
    log:
        "logs/{sample}/01_assembly/07_{chromosomes}.get-masking-sites.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	python "{params.depth_masking}" --depth "{config[depth]}" "{config[reference_genome]}" "{input.sorted_rg_bam}" "{output.mask_sites}"
        ) &> "{log}"
        """


rule generate_preconsensus_sequence:
    message: "Use coverage information to generate preconsensus sequence"
    input:
        vcf        = "{sample}/01_assembly/clair3/merge_output.vcf.gz",
        mask_sites = "{sample}/01_assembly/{chromosomes}_coverage_mask.txt"
    output:
        fasta      = "{sample}/01_assembly/{chromosomes}_preconsensus.fasta"
    params:
        mask       = os.path.join(workflow.basedir, "../workflow/scripts/01_mask-to-preconsensus.py")
    log:
        "logs/{sample}/01_assembly/08_{chromosomes}.generate-preconsensus-sequence.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
         r"""
         (
         	bcftools index "{input.vcf}"
         	python "{params.mask}" "{config[reference_genome]}" "{input.mask_sites}" "{input.vcf}" "{output.fasta}"
         ) &> "{log}"
         """


rule polish_consensus:
    message: "Polish sequence for each header in reference FASTA"
    input:
        vcf         = "{sample}/01_assembly/clair3/merge_output.vcf.gz",
        fasta       = "{sample}/01_assembly/{chromosomes}_preconsensus.fasta"
    output:
        fasta       = "{sample}/01_assembly/{chromosomes}_polished.fasta",
        normalised  = "{sample}/01_assembly/{chromosomes}_normalised.vcf.gz"
    params:
        chromosomes = "{chromosomes}",
    log:
        "logs/{sample}/01_assembly/09_{chromosomes}.polish-consensus.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	bcftools norm -f "{config[reference_genome]}" -m-any -o "{output.normalised}" "{input.vcf}"
        	bcftools index "{output.normalised}"
        	bcftools consensus -f "{config[reference_genome]}" -o "{output.fasta}" "{output.normalised}"
        ) &> "{log}"
        """ 

        
rule concatenate_FASTAs:
    message: "Concatenate all headers into one FASTA per sample with chromosome name appended"
    input:
        lambda wildcards: expand("{sample}/01_assembly/{chromosomes}_polished.fasta", 
                                 sample = wildcards.sample, 
                                 chromosomes = chromosomes)
    output:
        fasta = "{sample}/01_assembly/{sample}.fasta"
    log:
        "logs/{sample}/01_assembly/10_FASTA-concatenation.snakemake.log"
    shell:
        r"""
        (
        	for fasta in "{input}"; do
        		sample_id=$(dirname "$fasta" | cut -d '/' -f1)
            	chrom=$(basename "$fasta" | cut -d'_' -f1)
            	awk -v s="$sample_id" -v c="$chrom" 'BEGIN{{OFS=""}} /^>/ {{$0 = ">" s "_" c}} {{print}}' "$fasta"
        	done > "{output.fasta}"
        ) &> "{log}"
        """


rule run_bamtools_split:
    message: "Run bamtools split"
    input:
        sorted_rg_bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam"
    output:
        mapped        = temp("{sample}/01_assembly/{sample}_ref.sorted.rg.REF_{chromosomes}.bam")
    log:
        "logs/{sample}/01_assembly/11_{chromosomes}.bamtools-split.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	bamtools split -in "{input.sorted_rg_bam}" -reference
        ) &> "{log}"     
        """
    

rule run_picard_insert_metrics:
    message: "Get insert metrics using Picard"
    input:
        fasta         = "{sample}/01_assembly/{sample}.fasta",
        sorted_rg_bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
        mapped        = expand("{sample}/01_assembly/{sample}_ref.sorted.rg.REF_{chromosomes}.bam", 
                                 sample = sampleids, chromosomes = chromosomes)
    output:
        txt           = temp("{sample}/01_assembly/{sample}_insert_size_metrics.txt"),
        pdf           = temp("{sample}/01_assembly/{sample}_insert_size_histogram.pdf"),
    log:
        "logs/{sample}/01_assembly/12_picard-insert-metrics.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	if [ "$(awk 'END {{print NR}}' "{input.fasta}")" -eq 1 ]; then
            	touch {output.txt} {output.pdf}
        	else
            	picard CollectInsertSizeMetrics \
                	I="{input.sorted_rg_bam}" \
                	O="{output.txt}" \
                	H="{output.pdf}" \
                	M=0.5
        	fi
        ) &> "{log}"
        """    


rule run_picard_read_groups:
    message: "Get read groups using Picard"
    input:
        fasta         = "{sample}/01_assembly/{sample}.fasta",
        sorted_rg_bam = "{sample}/01_assembly/{sample}_ref.sorted.rg.bam",
    output:
        bam           = temp("{sample}/01_assembly/{sample}_ref.sorted.rg.fixed.bam")
    params:
        insample      = "{sample}"
    log:
        "logs/{sample}/01_assembly/13_picard-read-groups.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	if [ "$(awk 'END {{print NR}}' "{input.fasta}")" -eq 1 ]; then
            	touch "{output.bam}"
        	else
            	picard AddOrReplaceReadGroups \
                	I="{input.sorted_rg_bam}" \
                	O="{output.bam}" \
                	RGID=rg1 \
                	RGLB=lib1 \
                	RGPL=ILLUMINA \
                	RGPU=unit1 \
                	RGSM="{params.insample}"
        	fi
        ) &> "{log}"
        """


rule run_picard_duplicates:
    message: "Run Picard duplicates"
    input:
        fasta = "{sample}/01_assembly/{sample}.fasta",
        bam   = "{sample}/01_assembly/{sample}_ref.sorted.rg.fixed.bam",
        pdf   = "{sample}/01_assembly/{sample}_insert_size_histogram.pdf",
    output:
        bam   = temp("{sample}/01_assembly/{sample}_marked_duplicates.bam"),
        txt   = temp("{sample}/01_assembly/{sample}_marked_dup_metrics.txt")
    log:
        "logs/{sample}/01_assembly/14_picard-duplicates.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	if [ "$(awk 'END {{print NR}}' {input.fasta})" -eq 1 ]; then
            	touch "{output.txt}" "{output.bam}"
        	else
            	picard MarkDuplicates \
                	I="{input.bam}" \
                	O="{output.bam}" \
                	M="{output.txt}"
        	fi
        ) &> "{log}"
        """


rule run_picard_alignment_metrics:
    message: "Get alignment metrics using Picard"
    input:
        fasta = "{sample}/01_assembly/{sample}.fasta",
        bam   = "{sample}/01_assembly/{sample}_ref.sorted.rg.fixed.bam",
        txt   = "{sample}/01_assembly/{sample}_marked_dup_metrics.txt"
    output:
        txt   = temp("{sample}/01_assembly/{sample}_picard_output.txt")
    log:
        "logs/{sample}/01_assembly/15_picard-alignment-metrics.snakemake.log"
    conda:
        "../envs/env.yml"
    shell:
        r"""
        (
        	if [ "$(awk 'END {{print NR}}' {input.fasta})" -eq 1 ]; then
            	touch "{output.txt}"
        	else
            	picard CollectAlignmentSummaryMetrics \
                	R="{input.fasta}" \
                	I="{input.bam}" \
                	O="{output.txt}" 
        	fi
        ) &> "{log}"
        """

