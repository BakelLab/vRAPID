import os
import sys
import subprocess
import matplotlib
import argparse
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import statistics

def create_plots(sample_folder, amplified, threads, read1_suffix, read2_suffix, krakendb):
    repo_dir = args.repo_dir
    if args.virus == "Influenza-A" or args.virus == "Influenza-B":
        args.virus = args.virus.replace("-", " ")
    qc_dir = sample_folder + '/03_qualityControl'
    working_dir = sample_folder + '/02_assembly'
    sample = os.path.basename(sample_folder.rstrip('/'))
    if not os.path.exists(qc_dir):
        os.makedirs(qc_dir)
    
    subprocess.Popen("bamtools split -in %s/%s_ref.bam -reference" % (working_dir, sample), shell=True).wait()
    subprocess.Popen('''picard CollectInsertSizeMetrics \
            I=%s/%s_ref.bam \
            O=%s/%s_insert_size_metrics.txt \
            H=%s/%s_insert_size_histogram.pdf \
            M=0.5''' % (working_dir, sample, working_dir, sample, working_dir, sample), shell=True).wait()
    subprocess.Popen('''picard MarkDuplicates \
            I=%s/%s_ref.bam \
            O=%s/%s_marked_duplicates.bam \
            M=%s/%s_marked_dup_metrics.txt''' % (working_dir, sample, working_dir, sample, working_dir, sample), shell=True).wait()
    subprocess.Popen('''picard CollectAlignmentSummaryMetrics \
            R=%s/%s.fasta \
            I=%s/%s_ref.bam \
            O=%s/%s_picard_output.txt''' % (working_dir, sample, working_dir, sample, working_dir, sample),shell=True).wait()
    subprocess.Popen('''picard CollectGcBiasMetrics \
            I=%s/%s_ref.bam \
            O=%s/%s_picard_gc_bias_metrics.txt \
            CHART=%s/%s_picard_gc_bias_metrics.pdf \
            S=%s/%s_picard_summary_metrics.txt \
            R=%s/%s.fasta''' % (working_dir, sample, working_dir, sample, working_dir, sample, working_dir, sample, working_dir, sample),shell=True).wait()
    subprocess.Popen('''picard CollectRnaSeqMetrics \
            I=%s/%s_ref.bam \
            O=%s/%s_picard_gc_bias_metrics.txt \
            CHART=%s/%s_picard_rnaseq_metrics.pdf \
            S=%s/%s_picard_summary_metrics.txt \
            R=%s/%s.fasta''' % (working_dir, sample, working_dir, sample, working_dir, sample, working_dir, sample, working_dir, sample),shell=True).wait()
    
    
    chromosomes = []
    with open(args.reference, 'r') as f:
    	for line in f:
    		if line.startswith('>'):
    			chromosomes.append(line.strip()[1:])
    
    pp = PdfPages('%s/%s_quality_control2.pdf' % (qc_dir, sample))
    for chrs in chromosomes:
        outseq = ''
        count = 0
        
        
        # main coverage bit
        subprocess.Popen("samtools depth -aa %s/02_assembly/%s_ref.bam -r %s > %s/%s_coverage.txt" % (sample_folder, sample, chrs, qc_dir, chrs), shell=True).wait()
        with open("%s/%s_coverage.txt" % (qc_dir, chrs)) as f:
            cov1, cov2 = [], []
            for line in f:
                cov = int(line.split()[2])
                cov1.append(cov)
                cov2.append(min([100, cov]))
        fig, axs = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
        axs[0].plot(cov1)
        axs[1].plot(cov2)
        sample = os.path.basename(sample_folder.rstrip('/'))
        ns = 0
        length = 0
        
        with open(sample_folder + '/02_assembly/' + chrs + '.fasta') as f:
            for line in f:
                if not line.startswith('>'):
                    length += len(line.rstrip())
                    ns += line.lower().count('n')
        fig.suptitle("%s: all reads\nconsensus %d bp long\n%d Ns" % (chrs, length, ns), fontsize=10)
        with open('%s/%s_report.txt' % (qc_dir, chrs), 'w') as o:
            o.write('== all reads ==\n')
            o.write("length: %d\n" % length)
            o.write("Ns: %d\n" % ns)
            o.write("max coverage: %d\n" % max(cov1))
            o.write("median coverage: %d\n" % statistics.median(cov1))
        for ax in axs.flat:
            ax.set(xlabel='position', ylabel='coverage')
        for ax in axs.flat:
            ax.label_outer()
        plt.savefig(pp, dpi=300)
        
        with open("%s/04_variants/%s_pileup" % (sample_folder, chrs)) as f:
            count_dict = {'a': [], 't': [], 'c': [], 'g': [],'n': [], 'I': [], 'D': []}
            positions = []
            for line in f:
                ref, pos, refbase, cov, seq, qual = line.split()
                seq = seq.lower()
                refbase = refbase.lower()
                counts = {'a':0, 't':0, 'c':0, 'g':0, 'n':0, 'I':0, 'D':0}
                depth = 0
                seq = list(seq)
                getdel = False
                getins = False
                while seq != []:
                    x = seq.pop(0)
                    mod = None
                    if x == '.' or x == ',':
                        mod = refbase
                    elif x == '+':
                        getins = True
                        digit = ''
                    elif x == '-':
                        getdel = True
                        digit = ''
                    elif x.isdigit() and (getdel or getins):
                        digit += x
                    elif getdel:
                        for j in range(int(digit) -1):
                            seq.pop(0)
                        mod = 'D'
                        getdel = False
                    elif getins:
                        for j in range(int(digit) -1):
                            seq.pop(0)
                        mod = 'I'
                        getins = False
                    elif x == '^':
                        seq.pop(0)
                    elif x in ['a', 't', 'c', 'g']:
                        mod = x
                    elif x in ['$', '*']:
                        pass
                    else:
                        sys.exit("Base not recognised " + x)
                    if not mod is None:
                        counts[mod] += 1
                        depth += 1
                if  refbase == 'n':
                    continue
                if depth >= 10 and counts[refbase] /depth <= 0.85:
                    for i in counts:
                        count_dict[i].append(round(counts[i] /depth,2))
                    positions.append(pos + ' (' + refbase + ')\n['+str(depth)+']')
        fig, ax = plt.subplots()
        prev_key='start'
        for i in count_dict:
            if prev_key=='start':
                bottom_val=[0]*len(positions)
            else:
                bottom_val=[x + y for x, y in zip(bottom_val, count_dict[prev_key])]
            ax.bar(positions, count_dict[i],0.5, label=i,bottom=bottom_val)
            prev_key=i
        
        ax.set_ylabel('Proportions')
        ax.set_title(chrs +': Variants from reference')
        ax.set_xlabel('Position (ref. base)\n[read count]')
        plt.xticks(rotation=90)
        ax.legend()
        fig.subplots_adjust(bottom=0.3)
        plt.savefig(pp, dpi=300)
    	
        # by read coverage bit
        for suffix in os.listdir(sample_folder):
            read1 = None
            read2 = None
            if suffix.endswith('.fastq'):
                continue
            for i in os.listdir(os.path.join(sample_folder, suffix)):
                if i.endswith(read1_suffix):
                    read1 = os.path.join(sample_folder, suffix, i)
                if i.endswith(read2_suffix):
                    read2 = os.path.join(sample_folder, suffix, i)
            if not read1 is None and not read2 is None:
                count += 1
        if count == 1:
            for i in os.listdir(repo_dir + '/db'):
                if i.startswith("%s" % args.virus) and i.endswith(".tsv"):
                    with open(os.path.join(repo_dir, 'db', i)) as f:
                        labels = []
                        vals = []
                        xnums = []
                        f.readline()
                        num = 1
                        for line in f:
                            if line.split()[0].endswith("LEFT"):
                                label = line.split()[0][:-5]
                                l = int(line.split()[7])
                            elif line.split()[0].endswith("RIGHT"):
                                r = int(line.split()[7])
                                vals.append(statistics.median(cov1[l-1:r]))
                                xnums.append(num)
                                labels.append(label)
                                num += 1
                    fig, ax = plt.subplots()
                    ax.barh(xnums, vals)
                    plt.yticks(xnums, labels)
                    ax.set(ylabel="median depth", xlabel="primer set")
                    ax.tick_params(axis='both', which='minor', labelsize=8)
                    fig.suptitle("Median depth for primer set\n%s." % i[:-4], fontsize=10)
                    plt.savefig(pp, dpi=300)
        else:
            for suffix in os.listdir(sample_folder):
                read1 = None
                read2 = None
                for i in os.listdir(os.path.join(sample_folder, suffix)):
                    if i.endswith(read1_suffix):
                        read1 = os.path.join(sample_folder, suffix, i)
                    if i.endswith(read2_suffix):
                        read2 = os.path.join(sample_folder, suffix, i)
                if not read1 is None and not read2 is None:
                    if amplified:
                        subprocess.Popen(
                                "cutadapt -j %s -g file:%s -a file:%s"
                                " -G file:%s -A file:%s "
                                "-o %s/reads.1.fq.gz -p %s/reads.2.fq.gz %s %s > %s/cutadapt.log"
                                % (threads, args.forward_primer, args.reverse_primer, args.forward_primer, args.reverse_primer, qc_dir, qc_dir,
                                    read1, read2, qc_dir), shell=True).wait()
                        pilon_read_1 = "%s/reads.1.fq.gz" % qc_dir
                        pilon_read_2 = "%s/reads.2.fq.gz" % qc_dir
                    else:
                        pilon_read_1 = read1
                        pilon_read_2 = read2
                    subprocess.Popen(
                            "minimap2 -t %s -ax sr %s %s %s | samtools view -b | samtools sort -@ %s -o %s/%s_ref.bam -"
                            " && samtools index %s/%s_ref.bam && samtools depth -aa %s/%s_ref.bam > %s/%s_coverage.txt"
                            % (threads, args.reference, pilon_read_1, pilon_read_2, threads, qc_dir, sample, qc_dir, sample, qc_dir, sample, qc_dir, chrs), shell=True).wait()
                    
                    with open("%s/%s_coverage.txt" % (qc_dir, chrs)) as f:
                        cov1, cov2 = [], []
                        for line in f:
                            cov = int(line.split()[2])
                            cov1.append(cov)
                            cov2.append(min([100, cov]))
                    
                    fig, axs = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0})
                    axs[0].plot(cov1)
                    axs[1].plot(cov2)
                    fig.suptitle(suffix, fontsize=10)
                    for ax in axs.flat:
                        ax.set(xlabel='position', ylabel='coverage')
                    for ax in axs.flat:
                        ax.label_outer()
                    
                    plt.savefig(pp, dpi=300)
                    for i in os.listdir(repo_dir + '/db'):
                        if i.startswith("%s" % args.virus) and i.endswith(".tsv"):
                            with open(os.path.join(repo_dir, 'db', i)) as f:
                                labels = []
                                vals = []
                                xnums = []
                                f.readline()
                                num = 1
                                for line in f:
                                    if line.split()[0].endswith("LEFT"):
                                        label = line.split()[0][:-5]
                                        l = int(line.split()[7])
                                    elif line.split()[0].endswith("RIGHT"):
                                        r = int(line.split()[7])
                                        vals.append(statistics.median(cov1[l-1:r]))
                                        xnums.append(num)
                                        labels.append(label)
                                        num += 1
                            fig, ax = plt.subplots()
                            ax.barh(xnums, vals)
                            plt.yticks(xnums, labels)
                            ax.tick_params(axis='both', which='minor', labelsize=8)
                            ax.set(ylabel="median depth", xlabel="primer set")
                            fig.suptitle("Median depth for primer set\n %s \nand reads %s." % (i[:-4], suffix), fontsize=10)
                            plt.savefig(pp, dpi=300)
                    while True:
                        if cov1[0] >= 10:
                            break
                        cov1.pop(0)
                        if cov1 == []:
                            break
                    gotzero = False
                    while True:
                        count = 0
                        for i in cov1[-10:]:
                            if i < 10:
                                count += 1
                        if count < 3 and gotzero:
                            break
                        x = cov1.pop()
                        if x == 0:
                            gotzero = True
                        if cov1 == []:
                            break
                    ns = 0
                    for i in cov1:
                        if i <10:
                            ns += 1
                    with open('%s_/%s_report.txt' % (qc_dir, chrs), 'a') as o:
                        o.write("== %s ==\n" % suffix)
                        o.write("length: %d\n" % len(cov1))
                        o.write("Ns: %d\n" % ns)
                        if len(cov1) != 0:
                            o.write("max coverage: %d\n" % max(cov1))
                            o.write("median coverage: %d\n" % statistics.median(cov1))
                        else:
                            o.write("max coverage: 0\n")
                            o.write("median coverage: 0\n")
                            
    subprocess.Popen("samtools bam2fq -f 4 -@ %s %s/02_assembly/%s_ref.bam | gzip > %s/%s_kraken_input.fastq.gz" % (threads, sample_folder, sample, qc_dir, sample), shell=True).wait()
    subprocess.Popen("kraken2 --db %s --quick --report %s/%s_kraken_report.out --threads %s --output %s/%s_kraken --paired %s/%s/%s_*_R1_001.fastq.gz %s/%s/%s_*_R2_001.fastq.gz" % (krakendb, qc_dir, sample, threads, qc_dir, sample, sample, "01_fastqs", sample, sample, "01_fastqs", sample), shell=True).wait()
    subprocess.Popen("samtools flagstat  %s/02_assembly/%s_ref.bam > %s/%s_refbam.flagstat" % (sample_folder, sample, qc_dir, sample), shell=True).wait()
    
    pp.close()

    subprocess.Popen("%s -k %s/%s_kraken_report.out -o %s/taxanomic.pdf" % (args.taxanomic_breakdown, qc_dir, sample, qc_dir), shell=True).wait()
    for chrs in chromosomes:
        subprocess.Popen("%s -i %s/04_variants/%s.%s_variable_bases.tsv -o %s/%s_var.pdf" % (args.plot_coverage, sample, sample, chrs, qc_dir, chrs), shell=True).wait()
        subprocess.Popen("pdfunite %s/taxanomic.pdf %s/%s_quality_control2.pdf %s/*_var.pdf %s/%s_quality_control.pdf" % (qc_dir, qc_dir, sample, qc_dir, qc_dir, sample), shell=True).wait()
    subprocess.Popen("rm %s/%s_quality_control2.pdf" % (qc_dir, sample), shell=True).wait()
    subprocess.Popen("rm %s/taxanomic.pdf" % (qc_dir), shell=True).wait()
    for chrs in chromosomes:
        subprocess.Popen("rm %s/%s_var.pdf" % (qc_dir, chrs), shell=True).wait()


__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='vRAPID pipeline QC', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='QC for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n'
                                            'License: GPLv3\n'
                                            'USAGE: python run_QC.py -i sample1' % __version__)

parser.add_argument('-rd', '--repo_dir', action='store', help='path to repo dir')
parser.add_argument('-i', '--illumina_folder', action='store', help='Sample folder created by process_run.py')
parser.add_argument('-a', '--not_amplified', action='store_true', help="Skip cutadapt")
parser.add_argument('-t', '--threads', action='store', default="12", help='number of threads to use')
parser.add_argument('-r1', '--read1_suffix', action='store', default="_1.fastq.gz", help='suffix for finding read 1')
parser.add_argument('-r2', '--read2_suffix', action='store', default="_2.fastq.gz", help='suffix for finding read 2')
parser.add_argument('-kdb', '--kraken_db', action='store', default="/sc/arion/projects/vanbah01b/COVID/db/minikraken2_v2_8GB_201904_UPDATE", help='location of kraken database')
parser.add_argument('-r', '--reference', action='store', default="COVID.fa", help='reference genome for assembly')
parser.add_argument('-pf', '--forward_primer', action='store', default="SARS-CoV-2_primers_5prime_NI.fa", help='forward primer')
parser.add_argument('-pr', '--reverse_primer', action='store', default="SARS-CoV-2_primers_3prime_NI.fa", help='reverse primer')
parser.add_argument('-v', '--virus', action='store', default="SARS-CoV-2", help='virus/lineage: SARS-CoV-2, sCoV, MPOX, Influenza-A, Influenza-B')
parser.add_argument('-pc', '--plot_coverage', action='store', default="SARS-CoV-2", help='path to plotting coverage script')
parser.add_argument('-tb', '--taxanomic_breakdown', action='store', help='path to taxanomic breakdown script')

args = parser.parse_args()

if not args.illumina_folder is None:
    create_plots(args.illumina_folder, not args.not_amplified, args.threads, args.read1_suffix, args.read2_suffix, args.kraken_db)
