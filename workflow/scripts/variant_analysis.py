import sys, os, argparse, subprocess
import pysam

def run_variant_analysis(args):
    sample_folder = args.illumina_folder
    sample = os.path.basename(sample_folder.rstrip('/'))
    outdir = sample_folder + '/04_variants'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    repo_dir = args.repo_dir
    bam = "%s/02_assembly/%s_ref.bam" % (sample_folder, sample)
    ref_fasta = "%s" % args.reference
    subprocess.Popen("samtools mpileup -f %s %s | gzip > %s/pileup" % (ref_fasta, bam, outdir), shell=True).wait()
    modlist = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g', 'n', 'I', 'D']
    
    chromosomes = []
    with open(args.reference, 'r') as f:
    	for line in f:
    		if line.startswith('>'):
    			chromosomes.append(line.strip()[1:])
    print(chromosomes) 
    for chrs in chromosomes:
        outseq = ''
        changes = {}
        subprocess.Popen("samtools mpileup -f %s %s -r %s > %s/%s_pileup" % (ref_fasta, bam, chrs, outdir, chrs), shell=True).wait()
        with open("%s/02_assembly/%s_pilon.changes" % (sample_folder, chrs)) as f:
            for line in f:
                pos = line.split()[0].split(':')[1]
                ref_base = line.split()[2]
                new_base = line.split()[3]
                changes[pos] = (ref_base, new_base)
        primer_positions = set()
        with open(os.path.join("%s" % args.primer_set_2kb)) as f:
            f.readline()
            for line in f:
                name, seq, pool, length, tm, gc, start, end = line.split()
                start, end = int(start), int(end)
                for i in range(min([start, end]), max(start,end)):
                    primer_positions.add(i)
        with open(os.path.join("%s" % args.primer_set_1_5kb)) as f:
            f.readline()
            for line in f:
                name, seq, pool, length, tm, gc, start, end = line.split()
                start, end = int(start), int(end)
                for i in range(min([start, end]), max(start,end)):
                    primer_positions.add(i)
        with open("%s/%s_pileup" % (outdir, chrs)) as f, open("%s/%s.%s_variable_bases.tsv" % (outdir, sample, chrs), "w") as out:
            out.write("reference\tposition\tflagged\tin_primer\treference_base\tpilon_base\tdepth\treference_base_fraction\tforward_depth"
                    "\tforward_fraction\treverse_detph\treverse_fraction\tA\tT\tC\tG\ta\tt\tc\tg\tn\tinsertion\tdeletion\n")
            for line in f:
                ref, pos, refbase, cov, seq, qual = line.split()
                refbase = refbase.lower()
                if pos in changes and len(changes[pos][0]) == 1 and len(changes[pos][1]) == 1 and changes[pos][0] != '.':
                    if refbase != changes[pos][0].lower():
                        sys.exit("pilup doesn't match pilon output")
                    new_refbase = changes[pos][1].lower()
                elif pos in changes:
                    new_refbase = changes[pos][1].lower()
                else:
                    new_refbase = refbase
                counts = {'a':0, 't':0, 'c':0, 'g':0, 'A':0, 'T':0, 'C':0, 'G':0, 'n':0, 'I':0, 'D':0}
                depth = 0
                seq = list(seq)
                getdel = False
                getins = False
                forward_depth = 0
                reverse_depth = 0
                while seq != []:
                    x = seq.pop(0)
                    mod = None
                    if x == '.' :
                        mod = refbase.upper()
                        forward_depth += 1
                    elif x == ',':
                        mod = refbase.lower()
                        reverse_depth += 1
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
                    elif x in ['a', 't', 'c', 'g', 'A', 'T', 'C', 'G']:
                        mod = x
                        if x.islower():
                            reverse_depth += 1
                        else:
                            forward_depth += 1
                    elif x in ['$', '*']:
                        pass
                    else:
                        sys.exit("Base not recognised " + x)
                    if not mod is None:
                        counts[mod] += 1
                        depth += 1
                try:
                    fraction = (counts[new_refbase.upper()]+counts[new_refbase.lower()])/depth
                except KeyError:
                    continue
                except ZeroDivisionError:
                    fraction = None
                try:
                    forward_fraction = counts[new_refbase.upper()] / forward_depth
                except ZeroDivisionError:
                    forward_fraction = "no_cov"
                try:
                    reverse_fraction = counts[new_refbase.lower()] / reverse_depth
                except ZeroDivisionError:
                    reverse_fraction = "no_cov"
                if fraction is None:
                    flagged = "no_cov"
                    outseq += 'n'
                elif int(cov) < 10:
                    flagged = "low_cov"
                    outseq += 'n'
                elif pos in changes and (len(changes[pos][0]) != 1 or len(changes[pos][1]) != 1):
                    flagged = 'indel'
                elif (forward_fraction == "no_cov" or forward_fraction < args.min_ratio) and \
                        (reverse_fraction == "no_cov" or reverse_fraction < args.min_ratio):
                        flagged = "FLAGGED"
                        outseq += 'n'
                else:
                    flagged = "OK"
                    outseq += new_refbase
                if int(pos) in primer_positions:
                    inprimer = "Y"
                else:
                    inprimer = "N"
                outlist = [ref, pos, flagged, inprimer, refbase, new_refbase, cov, fraction, forward_depth, forward_fraction, reverse_depth, reverse_fraction]
                for i in modlist:
                    outlist.append(counts[i])
                out.write('\t'.join(map(str, outlist)) + '\n')
        with open("%s/%s.final.fna" % (outdir, chrs), "w") as out:
            out.write(">%s\n" % chrs)
            for i in range(0, len(outseq), 80):
                out.write(outseq[i:i+80] + '\n')


__version__ = "0.1.1"
parser = argparse.ArgumentParser(prog='vRAPID variant pipeline', formatter_class=argparse.RawDescriptionHelpFormatter,
                                description='QC for the assembly, mapping, base calling of viruses\n' \
                                            'Version: %s\n'
                                            'License: GPLv3\n'
                                            'USAGE: python run_QC.py -i sample1' % __version__)

parser.add_argument('-rd', '--repo_dir', action='store', help='path to repo dir')
parser.add_argument('-i', '--illumina_folder', action='store', help='Sample folder created by process_run.py')
parser.add_argument('-m', '--min_ratio', action='store', default=0.9, help='number of threads to use')
parser.add_argument('-r', '--reference', action='store', default="sars-cov-2/COVID.fa", help='reference genome for assembly')
parser.add_argument('-ps_2', '--primer_set_2kb', action='store', default="sars-cov-2/SARS-CoV-2_primers_2kb_set.tsv", help='2kb primer set')
parser.add_argument('-ps_1_5', '--primer_set_1_5kb', action='store', default="sars-cov-2/SARS-CoV-2_primers_1.5kb_set.tsv", help='1.5kb primer set')

args = parser.parse_args()

if not args.illumina_folder is None:
    run_variant_analysis(args)
