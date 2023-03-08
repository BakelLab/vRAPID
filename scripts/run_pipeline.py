#!/usr/bin/env python3
# coding: utf-8
import argparse
import subprocess
import sys, os


def run_illumina(args):
    read1s = []
    read2s = []
    working_dir = os.path.join(args.sample_folder, "02_assembly")
    sample = os.path.basename(args.sample_folder.rstrip('/'))
    ion_reads = None
    if not os.path.exists(working_dir):
    	os.makedirs(working_dir)
    for suffix in os.listdir(args.sample_folder):
    	if suffix == '.DS_Store':
    		pass
    	elif suffix.endswith('Z.fastq'):
    		ion_reads = os.path.join(args.sample_folder, suffix)
    	else:
    		for i in os.listdir(os.path.join(args.sample_folder, suffix)):
    			if i.endswith(args.read1_suffix):
    				read1s.append(os.path.join(args.sample_folder, suffix, i))
    			if i.endswith(args.read2_suffix):
    				read2s.append(os.path.join(args.sample_folder, suffix, i))
    if not ion_reads is None:
    	pass
    elif len(read1s) == 0 or len(read2s) == 0:
    	sys.exit("Couldn't find any reads in sample folder")
    else:
    	subprocess.Popen("cat %s  > %s/combined.1.fastq.gz" % (' '.join(read1s), working_dir), shell=True).wait()
    	subprocess.Popen("cat %s  > %s/combined.2.fastq.gz" % (' '.join(read2s), working_dir), shell=True).wait()
    if not ion_reads is None:
    	print("cutadapt -j %s -g file:%s -a file:%s -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log" % (args.threads, args.forward_primer, args.reverse_primer, working_dir, ion_reads, working_dir))
    	subprocess.Popen("cutadapt -j %s -g file:%s -a file:%s -o %s/reads.1.fq.gz %s > %s/cutadapt.1.log" % (args.threads, args.forward_primer, args.reverse_primer, working_dir, ion_reads, working_dir), shell=True).wait()
    	pilon_read_1 = "%s/reads.1.fq.gz" % working_dir
    	pilon_read_2 = ""
    elif not args.not_amplified:
    	subprocess.Popen("cutadapt -j %s -g file:%s -a file:%s"
    					" -G file:%s -A file:%s "
    					"-o %s/reads.1.fq.gz -p %s/reads.2.fq.gz %s/combined.1.fastq.gz %s/combined.2.fastq.gz > %s/cutadapt.1.log"
    					% (args.threads, args.forward_primer, args.reverse_primer, args.forward_primer, args.reverse_primer, working_dir, working_dir, working_dir, working_dir, working_dir), shell=True).wait()
    	pilon_read_1 = "%s/reads.1.fq.gz" % working_dir
    	pilon_read_2 = "%s/reads.2.fq.gz" % working_dir
    else:
    	pilon_read_1 = "%s/combined.1.fastq.gz" % working_dir
    	pilon_read_2 = "%s/combined.2.fastq.gz" % working_dir
        
    subprocess.Popen("fastqc -t %s --nogroup '%s/%s_reads.1.fq.gz' '%s/%s_reads.2.fq.gz' --outdir '%s'" % (args.threads, working_dir, sample, working_dir, sample, working_dir), shell=True).wait()
    subprocess.Popen("minimap2 -t %s -ax sr %s %s %s | samtools view -b | samtools sort -@ %s -o %s/%s_ref.bam -"
                     " && samtools index %s/%s_ref.bam"
                     % (args.threads, args.reference, pilon_read_1, pilon_read_2, args.threads, working_dir, sample, working_dir, sample), shell=True).wait()

    subprocess.Popen("mv %s/%s_ref.bam %s/%s_ref_clipped.bam"
                     " && mv %s/%s_ref.bam.bai %s/%s_ref_clipped.bam.bai"
                     " && samtools view -H %s/%s_ref_clipped.bam > %s/%s_tmp_header.bam"
                     " && samtools view %s/%s_ref_clipped.bam | awk -F '\t' -v OFS='\t' '$6!~/S/' > %s/%s_tmp_ref.bam"
                     " && cat %s/%s_tmp_header.bam %s/%s_tmp_ref.bam | samtools view -b | samtools sort -@ %s -o %s/%s_ref.bam -"
                     " && samtools index %s/%s_ref.bam"
                     " && rm -f %s/%s_tmp*bam"
                     % (working_dir, sample, working_dir, sample,
                        working_dir, sample, working_dir, sample,
                        working_dir, sample, working_dir, sample,
                        working_dir, sample, working_dir, sample,
                        working_dir, sample, working_dir, sample, working_dir, working_dir, sample,
                        working_dir, sample,
                        working_dir, sample), shell=True).wait()
    
    fastafiles = []
    chromosomes = []
    with open(args.reference, 'r') as f:
        for line in f:
            if line.startswith('>'):
                chromosomes.append(line.strip()[1:])
                fastafiles.append(line.strip()[1:]+".fasta")
    for chrs in chromosomes:
        subprocess.Popen("pilon --fix bases --changes --vcf --threads %s --mindepth 50 --genome %s --frags %s/%s_ref.bam --tracks --output %s/%s_pilon --targets %s"
                % (args.threads, args.reference, working_dir, sample, working_dir, chrs, chrs), shell=True).wait()
        with open('%s/%s_pilon.fasta' % (working_dir, chrs), 'r') as f:
            seq = ''
            for line in f:
                if not line.startswith('>'):
                    seq += line.rstrip()
        seq = list(seq)
        with open('%s/%s_pilon.changes' % (working_dir, chrs), 'r') as f:
            dels = set()
            ins = set()
            for line in f:
                if line.split()[2] == '.':
                    if '-' in line:
                        start, stop = map(int, line.split()[1].split(':')[1].split('-'))
                    else:
                        start = stop = int(line.split()[1].split(':')[1])
                    for num in range(start-1, stop):
                        ins.add(num)
                if line.split()[3] == '.':
                    if '-' in line:
                        start, stop = map(int, line.split()[0].split(':')[1].split('-'))
                    else:
                        start = stop = int(line.split()[0].split(':')[1])
                    for num in range(start-1, stop):
                        dels.add(num)
        with open('%s/%s_pilonCoverage.wig' % (working_dir, chrs), 'r') as f:
            f.readline()
            f.readline()
            qnum = 0
            for refnum, line in enumerate(f):
                while qnum in ins:
                    qnum += 1
                if refnum in dels:
                    continue
                if int(line.rstrip()) < 10:
                    seq[qnum] = 'n'
                qnum += 1
        seq = ''.join(seq)
        if seq.endswith('a'):
            seq = seq.rstrip('a')
        seq = seq.strip('n')
        while True:
            if len(seq) < 10:
                break
            if seq[-10:].count('n') < 3:
                break
            seq = seq[:-1]
            seq = seq.rstrip('n')
        with open('%s/%s.fasta' % (working_dir, chrs), 'w') as o:
            o.write(">%s\n" % chrs)
            for i in range(0, len(seq), 80):
                o.write(seq[i:i + 80] + '\n')
    

    with open('%s/%s.fasta' %(working_dir, sample), 'w') as o:
        if any(x in ["HA.fasta", "NA.fasta", "PA.fasta", "PB1.fasta", "PB2.fasta", "NP.fasta", "NS.fasta", "MP.fasta", "M.fasta"] for x in fastafiles):
            for names in fastafiles:
                name = working_dir+"/"+names
                with open(name) as infile:
                    o.write(infile.read())
                o.write("\n")
                        
        else:
            for names in fastafiles:
                name = working_dir+"/"+names
                with open(name) as infile:
                    o.write(">%s\n" % sample)
                    firstline = infile.readline()
                    sequence = infile.read()
                    o.write(sequence)
                o.write("\n")
               
           
    subprocess.Popen(
        "prokka --force --cpus %s --outdir %s/prokka --prefix %s --kingdom Viruses --proteins %s  %s/%s.fasta "
        % (args.threads, working_dir, sample, args.genbankfile, working_dir, sample), shell=True).wait()
    if not args.not_amplified and ion_reads is None:
        for v in virus_length:
            v = int(v)
            print(v)
            for r in ref_headers:
                subprocess.Popen("shovill --outdir %s/%s_shovill --R1 %s --R2 %s --gsize %d --cpus %s"
                        % (working_dir, r, pilon_read_1, pilon_read_2, v, args.threads), shell=True).wait()
    

__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='COVID pipeline', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='pipeline for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n'
                                            'License: GPLv3\n'
                                            'USAGE: python -i <sample_folder>'
                                            'Sample folder should look like the following'
                                            '<sample_folder>\n'
                                            '└───<reads_15kb_primers>\n'
                                            '│   │   <read_prefix>_1.fastq.gz\n'
                                            '│   │   <read_prefix>_2.fastq.gz\n'
                                            '└───<reads_2kb_primers>\n'
                                            '    │   <read_prefix>_1.fastq.gz\n'
                                            '    │   <read_prefix>_2.fastq.gz\n' % __version__)

parser.add_argument('-rd', '--repo_dir', action='store', help='path to repo dir')
parser.add_argument('-i', '--sample_folder', action='store', help='Sample folder created by process_run.py')
parser.add_argument('-o', '--working_dir', action='store', default="temp", help='working directory (only for CCS reads)')
parser.add_argument('-t', '--threads', action='store', default="4", help='number of threads to use')
parser.add_argument('-s', '--sample', action='store', help='sample name')
parser.add_argument('-v', '--version', action='store_true', help="print version and exit")
parser.add_argument('-a', '--not_amplified', action='store_true', help="Skip cutadapt and assembly")
parser.add_argument('-r1', '--read1_suffix', action='store', default="R1_001.fastq.gz", help='suffix for finding read 1')
parser.add_argument('-r2', '--read2_suffix', action='store', default="R2_001.fastq.gz", help='suffix for finding read 2')
parser.add_argument('-r', '--reference', action='store', default="COVID.fa", help='reference genome for assembly')
parser.add_argument('-pf', '--forward_primer', action='store', default="SARS-CoV-2_primers_5prime_NI.fa", help='five prime primers')
parser.add_argument('-pr', '--reverse_primer', action='store', default="SARS-CoV-2_primers_3prime_NI.fa", help='five prime primers')
parser.add_argument('-g', '--genbankfile', action='store', default="COVID.gbk", help='genbank file of reference genome')
parser.add_argument('-l', '--length', nargs='+', action='store', default="29903", help='length of reference genome')
parser.add_argument('-headers', '--headers', nargs='+', action='store', help='headers of reference genome')

args = parser.parse_args()
virus_length = args.length
ref_headers = args.headers

if args.version:
    sys.stdout.write('Version %s' % __version__)
    sys.exit()

repo_dir = args.repo_dir
#if not args.ccs_reads is None:
#    run_ccs(args)
#elif not args.thermo_fischer is None:
#    run_thermo(args)
#elif not args.sample_folder is None:
if not args.sample_folder is None:     
    run_illumina(args)
else:
    sys.exit("Need to provide with sample folder or ccs reads")
