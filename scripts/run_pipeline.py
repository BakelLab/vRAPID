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
        subprocess.Popen("shovill --outdir %s/shovill --R1 %s --R2 %s --gsize %d --cpus %s"
                     % (working_dir, pilon_read_1, pilon_read_2, virus_length, args.threads), shell=True).wait()

def run_thermo(args):
    working_dir = os.path.join(args.thermo_fischer, "02_assembly")
    illumina_bam_dir = os.path.join(args.thermo_fischer, "bams")
    sample = os.path.basename(args.thermo_fischer.rstrip('/'))
    thermo_bam = []
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    """
    for i in os.listdir(os.path.join(args.thermo_fischer, "bams")):
        if i.endswith(".bam"):
            thermo_bam.append(os.path.join(args.thermo_fischer, "bams", i))
    if len(thermo_bam) == 0:
        sys.exit("Couldn't find any bamss in sample folder")
    elif len(thermo_bam) == 1:
        subprocess.Popen("samtools index %s" % thermo_bam[0], shell=True).wait()
        if os.path.exists(illumina_bam_dir+'/ref.bam'):
            subprocess.Popen("pilon --targets 2019-nCoV --fix bases --changes --vcf --threads %s --mindepth 50 --genome "
                            "%s/db/Ion_AmpliSeq_SARS-CoV-2_reference.fa --frags %s/ref.bam --unpaired %s --tracks --output %s/pilon"
                        % (args.threads, repo_dir, illumina_bam_dir,thermo_bam[0], working_dir), shell=True).wait()
        else:
            subprocess.Popen("pilon --targets 2019-nCoV --fix bases --changes --vcf --threads %s --mindepth 50 --genome "
                            "%s/db/Ion_AmpliSeq_SARS-CoV-2_reference.fa --unpaired %s --tracks --output %s/pilon"
                        % (args.threads, repo_dir, thermo_bam[0], working_dir), shell=True).wait()
    else:
        subprocess.Popen("samtools merge %s > %s/concat.bam" % (" ".join(thermo_bam), working_dir), shell=True).wait()
        subprocess.Popen("samtools index %s/concat.bam" % working_dir, shell=True).wait()
        if os.path.exists(illumina_bam_dir+'/ref.bam'):
            subprocess.Popen("pilon --targets 2019-nCoV --fix bases --changes --vcf --threads %s --mindepth 50 --genome "
                            "%s/db/Ion_AmpliSeq_SARS-CoV-2_reference.fa --frags %s/ref.bam --unpaired %s/concat.bam --tracks --output %s/pilon"
                            % (args.threads, repo_dir, illumina_bam_dir, working_dir, working_dir), shell=True).wait()

        else:
            subprocess.Popen("pilon --targets 2019-nCoV --fix bases --changes --vcf --threads %s --mindepth 50 --genome "
                            "%s/db/Ion_AmpliSeq_SARS-CoV-2_reference.fa --unpaired %s/concat.bam --tracks --output %s/pilon"
                            % (args.threads, repo_dir, working_dir, working_dir), shell=True).wait()
    """
    subprocess.Popen("samtools view -@ %s -h -q 80 %s/*%s.bam | samtools view -@ %s -bS > %s/%s_mapq80.bam" %  (args.threads,illumina_bam_dir,sample,args.threads,illumina_bam_dir,sample), shell=True).wait()
    subprocess.Popen("samtools bam2fq %s/%s_mapq80.bam > %s/%s_mapq80.fastq" %  (illumina_bam_dir,sample,illumina_bam_dir,sample), shell=True).wait()
    subprocess.Popen("samtools bam2fq %s/*%s.bam > %s/%s.fastq" %  (illumina_bam_dir,sample,illumina_bam_dir,sample), shell=True).wait()
    #subprocess.Popen("samtools view -bS %s/%s_mapq80.sam > %s/%s_mapq80.bam" % (illumina_bam_dir,sample,illumina_bam_dir,sample), shell=True).wait()
    #subprocess.Popen("cutadapt -u 20 -u -20 -o %s/%s_2019-nCoV_trimmed20.fastq  %s/%s_2019-nCoV.fastq" %  (illumina_bam_dir,sample,illumina_bam_dir,sample), shell=True).wait()
    #subprocess.Popen("trimmomatic SE -phred33 %s/%s_2019-nCoV.fastq %s/%s_2019-nCoV_trimmed.fastq SLIDINGWINDOW:4:15 MINLEN:36" %  (illumina_bam_dir,sample,illumina_bam_dir,sample), shell=True).wait()
    subprocess.Popen("minimap2 -t %s -ax sr %s %s/%s_mapq80.fastq | samtools view -b | samtools sort -@ %s -o %s/ref_mapq80.bam -"
                     " && samtools index %s/ref_mapq80.bam"
                     % (args.threads, args.reference, illumina_bam_dir,sample, args.threads, working_dir, working_dir), shell=True).wait()
    subprocess.Popen("minimap2 -t %s -ax sr %s %s/%s.fastq | samtools view -b | samtools sort -@ %s -o %s/%s_ref.bam -"
                     " && samtools index %s/%s_ref.bam"
                     % (args.threads, args.reference, illumina_bam_dir,sample, args.threads, working_dir, sample, working_dir, sample), shell=True).wait()
    #subprocess.Popen("samtools view -h -q 50 %s/ref_trimmed20.bam > %s/ref_trimmed20_mapq50.sam" % (working_dir,working_dir), shell=True).wait()
    #subprocess.Popen("samtools view -bS %s/ref_trimmed20_mapq50.sam > %s/ref_trimmed20_mapq50.bam" % (working_dir,working_dir), shell=True).wait()
    #subprocess.Popen("samtools index %s/ref_trimmed20_mapq50.bam" % (working_dir), shell=True).wait()
    #subprocess.Popen("samtools view -ub %s/ref.bam 2019-nCoV | samtools bam2fq - > %s/%s_2019-nCoV.fastq" % (working_dir,working_dir,sample), shell=True).wait()
    #subprocess.Popen("awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' %s/%s_2019-nCoV.fastq > %s/%s_2019-nCoV.meanlength" % working_dir,sample, shell=True).wait()
    subprocess.Popen("pilon --fix bases --changes --vcf --threads %s --mindepth 50 --genome "
                            "%s --unpaired %s/ref_mapq80.bam --tracks --output %s/%s_pilon"
                        % (args.threads, args.reference, working_dir, working_dir,working_dir, sample), shell=True).wait()

    with open('%s/%s_pilon.fasta' % (working_dir, sample)) as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.rstrip()
    seq = list(seq)
    with open('%s/%s_pilon.changes' % (working_dir, sample)) as f:
        dels = set()
        ins = set()
        for line in f:

            if line.split()[2] == '.':
                if '-' in line.split()[1].split(':')[1]:
                    start, stop = map(int, line.split()[1].split(':')[1].split('-'))
                else:
                    start = stop = int(line.split()[1].split(':')[1])
                for num in range(start-1, stop):
                    ins.add(num)

            if line.split()[3] == '.':
                if '-' in line.split()[0].split(':')[1]:
                    start, stop = map(int, line.split()[0].split(':')[1].split('-'))
                else:
                    start = stop = int(line.split()[0].split(':')[1])
                for num in range(start-1, stop):
                    dels.add(num)


    with open('%s/%s_pilonCoverage.wig' % (working_dir, sample)) as f:
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
    with open('/%s.fasta' % (working_dir, sample), 'w') as o:
        o.write(">%s\n" % sample)
        for i in range(0, len(seq), 80):
            o.write(seq[i:i + 80] + '\n')
    subprocess.Popen(
        "prokka --force --cpus %s --outdir %s/prokka --prefix %s --kingdom Viruses --proteins %s  %s/%s.fasta "
        % (args.threads, working_dir, sample, args.genbankfile, working_dir, sample), shell=True).wait()

def run_ccs(args):
    subprocess.Popen("cutadapt -j %s -g file:%s/db/SARS-CoV-2_primers_5prime_anchored.fa -a "
                     "file:%s/db/SARS-CoV-2_primers_3prime_anchored.fa -o %s/%s_reads.1.fq.gz %s > %s/%s_cutadapt.1.log"
                     % (args.threads, repo_dir, repo_dir, working_dir, sample, args.ccs_reads, working_dir, sample), shell=True).wait()
    with open(working_dir + '/%s_cutadapt.1.log' % sample) as f:
        for line in f:
            if line.startswith("Total written (filtered):"):
                bp = int(line.split()[3].replace(',', ''))
                break
    if bp / virus_length > args.coverage_pilon:
        downsample = args.coverage_pilon / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/%s_reads.1.fq.gz %f | gzip > %s/reads.pilon.1.fq.gz"
                         % (working_dir, sample, downsample, working_dir), shell=True).wait()
        pilon_read_1 = "%s/reads.pilon.1.fq.gz" % working_dir
    else:
        pilon_read_1 = "%s/%s_reads.1.fq.gz" % (working_dir, sample)
    subprocess.Popen("minimap2 -t %s -ax map-pb %s %s | samtools view -b | samtools sort -@ %s -o %s/%s_ref.bam -"
                     " && samtools index %s/%s_ref.bam"
                     % (args.threads, args.reference, pilon_read_1, args.threads, working_dir, sample, working_dir, sample), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s --unpaired %s/%s_ref.bam --tracks --output %s/%s_pilon"
                     % (args.threads, args.reference, working_dir, sample, working_dir, sample), shell=True).wait()
    with open('%s/%s_pilon.fasta' % (working_dir, sample)) as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.rstrip()
    seq = list(seq)
    with open('%s/%s_pilonCoverage.wig' % (working_dir, sample)) as f:
        f.readline()
        f.readline()
        for num, line in enumerate(f):
            if int(line.rstrip()) < 20:
                seq[num] = 'n'
    seq = ''.join(seq)
    seq = seq.strip('n')
    with open(working_dir + '/%s.fasta' % args.sample, 'w') as o:
        o.write(">%s\n" % args.sample)
        for i in range(0, len(seq), 80):
            o.write(seq[i:i+80] + '\n')
    subprocess.Popen("prokka --force --cpus %s --outdir %s/prokka --prefix %s --kingdom Viruses --proteins %s  %s/%s.fasta "
                     % (args.threads, working_dir, args.sample, repo_dir, working_dir, args.genbankfile, args.sample), shell=True).wait()
    if bp / virus_length > 60:
        downsample = 60 / (bp / virus_length)
        subprocess.Popen("seqtk sample %s/%s_reads.1.fq.gz %f | gzip > %s/reads.canu.1.fq.gz"
                         % (working_dir, sample, downsample, working_dir), shell=True).wait()
        canu_reads_1 = "%s/reads.canu.1.fq.gz" % working_dir
    else:
        canu_reads_1 = "%s/%s_reads.1.fq.gz" % (working_dir, sample)
    subprocess.Popen("seqkit rmdup -s %s > %s/rmdup.fastq" % (canu_reads_1, working_dir), shell=True).wait()
    subprocess.Popen("canu -d %s/canu -pacbio-corrected %s/rmdup.fastq -p canu genomeSize=%d useGrid=false minOverlapLength=250"
                    % (working_dir, working_dir, virus_length), shell=True).wait()
    subprocess.Popen("minimap2 -t %s -ax map-pb %s/canu/canu.contigs.fasta %s | samtools view -b | samtools sort -@ %s -o %s/assembly.bam -"
                     " && samtools index %s/assembly.bam"
                     % (args.threads, working_dir, pilon_read_1, args.threads, working_dir, working_dir), shell=True).wait()
    subprocess.Popen("pilon --fix bases --threads %s --mindepth 20 --genome %s/canu/canu.contigs.fasta --unpaired %s/assembly.bam --tracks --output %s/assembly_pilon"
                     % (args.threads, working_dir, working_dir, working_dir), shell=True).wait()
    

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
parser.add_argument('-b', '--thermo_fischer', action='store', help='Sample folder with thermofischer bams present')
parser.add_argument('-p', '--ccs_reads', action='store', help='Pacbio CCS reads')
parser.add_argument('-o', '--working_dir', action='store', default="temp", help='working directory (only for CCS reads)')
parser.add_argument('-t', '--threads', action='store', default="4", help='number of threads to use')
parser.add_argument('-s', '--sample', action='store', help='sample name')
parser.add_argument('-c', '--coverage_pilon', default=200, type=int, action='store', help='downsample to this coverage for pilon (only used for CCS reads)')
parser.add_argument('-v', '--version', action='store_true', help="print version and exit")
parser.add_argument('-a', '--not_amplified', action='store_true', help="Skip cutadapt and assembly")
parser.add_argument('-r1', '--read1_suffix', action='store', default="R1_001.fastq.gz", help='suffix for finding read 1')
parser.add_argument('-r2', '--read2_suffix', action='store', default="R2_001.fastq.gz", help='suffix for finding read 2')
parser.add_argument('-r', '--reference', action='store', default="COVID.fa", help='reference genome for assembly')
parser.add_argument('-pf', '--forward_primer', action='store', default="SARS-CoV-2_primers_5prime_NI.fa", help='five prime primers')
parser.add_argument('-pr', '--reverse_primer', action='store', default="SARS-CoV-2_primers_3prime_NI.fa", help='five prime primers')
parser.add_argument('-g', '--genbankfile', action='store', default="COVID.gbk", help='genbank file of reference genome')
parser.add_argument('-l', '--length', action='store', default="29903", help='length of reference genome')

args = parser.parse_args()
virus_length = int(args.length)

if args.version:
    sys.stdout.write('Version %s' % __version__)
    sys.exit()

repo_dir = args.repo_dir
if not args.ccs_reads is None:
    run_ccs(args)
elif not args.thermo_fischer is None:
    run_thermo(args)
elif not args.sample_folder is None:
    run_illumina(args)
else:
    sys.exit("Need to provide with sample folder or ccs reads")
