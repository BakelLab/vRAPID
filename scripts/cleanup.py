#!/usr/bin/env python3

############
# MODULES #
###########

import os
import pandas as pd
import shutil
import glob
from argparse import ArgumentParser

##########
# FILES #
#########

parser = ArgumentParser(description="Parse lineage assignments and versions from lineage tables")
parser.add_argument('-p', '--inpcsv', help='File containing samples in the run in csv format', required=True)
args = parser.parse_args()

inFile = args.inpcsv
mappings = pd.read_csv(inFile)
test = mappings.set_index('Sample_ID')


samples = mappings["Sample_ID"].tolist()

extensions = ('_pilon.changes', '_pilon.fasta', '_pilon.vcf', '_pilonPilon.bed', '.fasta', '.1.log', 'prokka', 'reads.1.fq.gz', 'reads.2.fq.gz', 'shovill', '_insert_size_histogram.pdf', '_insert_size_metrics.txt', '_marked_dup_metrics.txt', '_marked_duplicates.bam', '_picard_output.txt', '_ref_stats', '_ref.bam', '_ref.bam.bai')



try:
    for s in samples:
        s=s+'/02_assembly/'
        for f in list(os.listdir(s)):
            if not f.endswith(extensions):
                os.remove(s+f)

except OSError as e:
	print(e)
