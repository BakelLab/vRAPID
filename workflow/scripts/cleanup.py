#!/usr/bin/env python3

############
# MODULES #
###########

import os
import pandas as pd
import shutil
import glob
from argparse import ArgumentParser
import subprocess
import gzip

##########
# FILES #
#########

parser = ArgumentParser(description="Clean up after pipeline is done")
parser.add_argument('-p', '--inpcsv', help='File containing samples in the run in csv format', required=True)
args = parser.parse_args()

inFile = args.inpcsv
mappings = pd.read_csv(inFile)
test = mappings.set_index('Sample_ID')


samples = mappings["Sample_ID"].tolist()

extensions = ('_pilon.changes', '_pilon.fasta', '_pilon.vcf', '_pilonPilon.bed', '.fasta', '.1.log', 'prokka', 'reads.1.fq.gz', 'reads.2.fq.gz', 'shovill', '_ref_stats', '_ref.bam', '_ref.bam.bai', '_annotation')



try:
    for s in samples:
        assembly=s+'/02_assembly/'
        for f in list(os.listdir(assembly)):
            if not f.endswith(extensions):
                os.remove(assembly+f)
        var=s+"/04_variants"
        for f in list(os.listdir(var)):
            if f.endswith("pileup"):
                v = var+"/"+f
                vgz = var+"/"+f+".gz"
                with open(v, 'rb') as f_in:
                    with gzip.open(vgz, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(v)
            
except OSError as e:
	print(e)
