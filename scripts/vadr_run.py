import sys, subprocess
fasta_file = sys.argv[1]
output_dir = sys.argv[2]
out_gff = sys.argv[3]
vadr_db = sys.argv[4]

subprocess.Popen("v-annotate.pl --split --cpu 8 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir %s %s %s --noseqnamemax" % (vadr_db, fasta_file, output_dir), shell=True).wait()

with open(out_gff, 'w') as o:
    o.write("##gff-version 3\n")
    length = 0
    with open(fasta_file) as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.rstrip())
    with open(output_dir + '/VADR'+ '.vadr.pass.tbl') as f:
        for line in f:
            if line.startswith(">Feature"):
                contig = line.split()[1]
                o.write("##sequence-region %s 1 %d" % (contig, length))
            elif not line.startswith("\t"):
                splitline = line.rstrip().split("\t")
                if len(splitline) == 3:
                    feature = splitline[2]
                else:
                    feature = "misc_feature"                     
                start, stop = splitline[:2]
                stop=stop.replace('>','')
                if int(start) < int(stop):
                    strand = "+"
                else:
                    strand = "-"
                o.write("\n" + "\t".join([contig, "vadr", feature, start, stop, strand, str(int(start)%3)]))
                first = True
            else:
                splitline = line.rstrip().split("\t")
                key, value = splitline[3:5]
                if first:
                    o.write("\t")
                    first = False
                else:
                    o.write(";")
                o.write(key + "=" + value)
    o.write("\n##FASTA\n")
    with open(fasta_file) as f:
        for line in f:
            o.write(line)
           
