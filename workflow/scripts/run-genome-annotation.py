#!/usr/bin/env python

import requests
import time
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna
except ImportError:
    generic_dna = None

import argparse
import subprocess
import sys, os

# Arguments parsing
parser = argparse.ArgumentParser(description="Annotate FASTA Sequence if applicable")
parser.add_argument("-i", "--input", required=True, help="Name of fasta file to annotate")
parser.add_argument("-v", "--virus", required=True, default="SARS-CoV-2", help="Virus to run annotator")
parser.add_argument("-u", "--url", default="https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/api/", help="Url for the flu annotation tool. Default: https://www.ncbi.nlm.nih.gov/genomes/FLU/annotation/api/")
parser.add_argument("-p", "--proxy", default="http://nfs01.chimera.hpc.mssm.edu:3128", help="Web proxy url. Default: http://nfs01.chimera.hpc.mssm.edu:3128")
parser.add_argument("-r", "--retries", type=int, default=35, help="Number of times to attempt data retrieval. Default: 35")

args = parser.parse_args()

# function to download results from flu annotator rule
def curl_download(sJobID, sFormat, sProxy, nRetries):
    # get job ID
    sJobID = sJobID[11:-2]
    sCurlResults = "curl -X POST -F 'cmd=download' -F 'job_id=" + sJobID + "' -F 'format=" + sFormat + "' " + sURL
    # tu run when using a proxy
    if sProxy:
        sCurlResults += " --proxy " + sProxy

    sResult = ""
    while nRetries:
        time.sleep(15)
        try:
            result = subprocess.run(sCurlResults, stdout=subprocess.PIPE, shell=True)
            sResult = result.stdout.decode("utf-8")
        except subprocess.CalledProcessError as e:
            print("Error: could not open curl: ", e)
            return ""
        # job status check
        if "Job is not complete. Status is" in sResult:
            nRetries -= 1
            sResult = ""
        else:
            nRetries = 0
    return sResult


def fasta_to_hash(sFasta):
    hFasta = {}
    sFastaHeader = ""
    sFastaSeq = ""
    asFasta = sFasta.split("\n")
    for i in range(len(asFasta)):
        if asFasta[i].startswith(">"):
            if i == len(asFasta) - 1:
                raise Exception("Error: file ends in fasta header without sequence")
            sFastaSeq = sFastaSeq.replace(" ", "")
            if sFastaHeader in hFasta:
                raise Exception(f"Error: duplicate fasta ID '{sFastaHeader}'")
            if sFastaHeader:
                hFasta[sFastaHeader] = sFastaSeq

            # Reset for the next sequence
            sFastaHeader = asFasta[i][1:].strip().split(" ")[0]
            sFastaSeq = ""

        elif i == len(asFasta) - 1:
            sFastaSeq += asFasta[i]
            sFastaSeq = sFastaSeq.replace(" ", "")
            if sFastaHeader in hFasta:
                raise Exception(f"Error: duplicate fasta ID '{sFastaHeader}'")
            hFasta[sFastaHeader] = sFastaSeq
        else:
            if not asFasta[i].strip() or asFasta[i].startswith("#"):
                continue
            if sFastaHeader:
                sFastaSeq += asFasta[i]
    return hFasta


def get_protein_id_mapping(sTable):
    hReturn = {}
    sGene, sID = "", ""
    asLines = sTable.split("\n")
    for sLine in asLines:
        if sLine == "CDS\t\t":
            sGene, sID = "", ""
        if "\t\t\tprotein_id\t" in sLine:
            sID = sLine.split("\t\t\tprotein_id\t")[1]
        if "\t\t\tgene\t" in sLine:
            sGene = sLine.split("\t\t\tgene\t")[1]
        if sGene and sID:
            if sID in hReturn:
                raise Exception("Error: found duplicate protein ID in table file")
            hReturn[sID] = sGene
            sGene, sID = "", ""
    return hReturn


def translate(sSeq):
    seq = Seq(sSeq, generic_dna)
    aaseq = seq.translate()
    aaReturn = ""
    for i in range(0, len(aaseq), 100):
        aaReturn += str(aaseq[i:i+100]) + "\n"
    return aaReturn

if __name__ == "__main__":
    
    # Set up variables from command line input
    fasta = args.input
    sample = args.input.split(".")[0].split("/")[2]
    output = args.input.split(".")[0] + "_trimmed.fa"
    status = args.input.split(".")[0] + "_genome_annotation.txt"
    
    # If the virus is SARS-CoV-2, run genome annotation using VADR
    if args.virus == "SARS-CoV-2":
        vadr_out = args.input.split(".")[0]+"_annotation"
        print(vadr_out)
        
        # Read fasta file, and open output file 
        with open(fasta, 'r') as f, open(status, "w") as o:
            lines = [i for i in f.readlines() if len(i) >= 1]
            # if the file is only a header, it will fail annotation
            if len(lines) <= 1:
                o.write("No sequence detected. Failed Assembly!")
            else:
                subprocess.Popen("$VADRSCRIPTSDIR/miniscripts/fasta-trim-terminal-ambigs.pl --minlen 50 --maxlen 30000 %s > %s" % (args.input, output), shell=True).wait()
                subprocess.Popen("v-annotate.pl -f --split --cpu 8 --glsearch --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn %s %s --noseqnamemax" % (output, vadr_out), shell=True).wait()
                o.write("Genome annotation using VADR for"+sample+" is complete!")

    elif args.virus == "Influenza-A" or args.virus == "Influenza-B":
        # Submit FASTA to IRD and download results
        sJobID = ""
        sProxy = args.proxy
        nRetries = args.retries
        sURL = args.url
        sOutputPrefix = args.input.rsplit(".", 1)[0]
        with open(args.input, "rb") as f:
            files = {"sequence": f}
            data = {"cmd": "submit"}
            submission_retries = nRetries
            while submission_retries:
                r = requests.post(args.url, files=files, data=data, proxies={"http": sProxy})
                for line in r.text.splitlines():
                    if "JSID" in line:
                        sJobID = line.strip()
                        print(sJobID)
                        break
                if sJobID:
                    submission_retries = 0
                else:
                    submission_retries -= 1
            if sJobID:
                print("Job submitted successfully with job ID:", sJobID)
                print(sJobID)
                sResultTable = curl_download(sJobID, 'tbl', sProxy, nRetries)
                sResultCDS = curl_download(sJobID, 'ffn', sProxy, nRetries)
                if sResultTable and sResultCDS:
                    rProteinMappings = get_protein_id_mapping(sResultTable)
                    rCDS = fasta_to_hash(sResultCDS)
                    try:
                        with open("{}.features_table.txt".format(sOutputPrefix), 'w') as outtable:
                            outtable.write(sResultTable)
                    except:
                        raise Exception("Error: can't open '{}.features_table.txt': ".format(sOutputPrefix) + str(sys.exc_info()[1])) 
                    try:
                        with open("{}.features_cds.fa".format(sOutputPrefix), 'w') as outcds:
                            for sID in sorted(rCDS.keys()):
                                if sID in rProteinMappings:
                                    sHeader = ">" + rProteinMappings[sID] + "|" + sID
                                    sSequence = ""
                                    nSeqLength = len(rCDS[sID])
                                    for i in range(0, nSeqLength, 100):
                                        sSequence += rCDS[sID][i:i+100] + "\n"
                                    outcds.write(sHeader + "\n" + sSequence)
                                else:
                                    raise Exception("Error: can't match protein ID '" + sID + "' to a product name")
                    except:
                        raise Exception("Error: can't open '{}.features_cds.fa': ".format(sOutputPrefix) + str(sys.exc_info()[1]))
                    
                    # Print protein sequence
                    try:
                        with open("{}.features_protein.fa".format(sOutputPrefix), 'w') as outprot:
                            for sID in sorted(rCDS.keys()):
                                if sID in rProteinMappings:
                                    sHeader = ">" + rProteinMappings[sID] + "|" + sID
                                    sSequence = translate(rCDS[sID])
                                    outprot.write(sHeader + "\n" + sSequence)
                                else:
                                    raise Exception("Error: can't match protein ID '" + sID + "' to a product name")
                    except:
                        raise Exception("Error: can't open '{}.features_protein.fa': ".format(sOutputPrefix) + str(sys.exc_info()[1]))
        
        # Get Influenza subtype
        h_serotype = {}
        n_serotype_count = 0
        s_serotype = ""
        s_input_file = sOutputPrefix+".features_table.txt"
        try:
            with open(s_input_file, 'r') as f:
                if args.virus == "Influenza-A":
                    print(args.virus)
                    for line in f:
                        if line.strip() == '' or line.startswith('#'):
                            continue
                        if 'INFO: Serotype: ' in line:
                            serotype = line.split('INFO: Serotype: ')[1]
                            h_key = serotype[0]
                            h_value = int(serotype[1:])
                            if h_key in h_serotype:
                                h_serotype[h_key].append(h_value)
                            else:
                                h_serotype[h_key] = [h_value]
                            n_serotype_count += 1
        
                    if n_serotype_count > 2:
                        s_serotype = "Mixed"
                    else:
                        s_h = h_serotype.get("H", ["x"])[0]
                        s_n = h_serotype.get("N", ["x"])[0]
                        s_serotype = f"H{s_h}N{s_n}"
            
                if args.virus == "Influenza-B":
                    for line in f:
                        if ">Feature HA" in line:
                            line = line.rstrip().lower().replace(" ", "|").split("|")
                            s_serotype = line[1]
                            if s_serotype == "ha" or s_serotype == "NA":
                                s_serotype = "NA"
                            else:
                                s_serotype = s_serotype.upper()
        except FileNotFoundError:
            with open(s_input_file, 'w') as f:
                f.write("Error with genome annotation!")
        with open("{}_subtype.txt".format(sOutputPrefix), 'w') as f, open(status, "w") as g:
            f.write(f"{sOutputPrefix}\t{s_serotype}\n")
            g.write("Genome annotation using IRD for "+sample+" is complete!")

    else:
        # If virus doesn't have an annotation tool, write to file
        with open(status, "w") as f:
            f.write("No genome annotation conducted. "+ sample+" is a(n) "+ args.virus+ " sample")

