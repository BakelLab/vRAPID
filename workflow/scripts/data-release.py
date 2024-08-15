import os
import tarfile
import pandas as pd
import mysql.connector
import csv
from argparse import ArgumentParser

# setup to collect all config info for database connection
passpath = os.path.expanduser('~') + '/.my.cnf.pdbrw'
with open(passpath) as cnf_file:
    for line in cnf_file:
        if line.startswith('user='):
            user = line.rstrip()[5:]
        if line.startswith('password='):
            password = line.rstrip()[9:]
        if line.startswith('host='):
            host = line.rstrip()[5:]
        if line.startswith('database='):
            database = line.rstrip()[9:]

# Parse arguments
parser = ArgumentParser(description="Sample names in run")
parser.add_argument('-p', '--inpcsv', help='Sample name in run in CSV format', required=True)
parser.add_argument('-r', '--runid', help='TD number of the assembly run name', required=True)
parser.add_argument('-v', '--virus', nargs='+', default="SARS-CoV-2", help='Virus in run', required=True)
parser.add_argument('-f', '--fastaheaders', nargs='+', default="MN908947", help='Headers in reference FASTA file', required=True)
args = parser.parse_args()

mappings = pd.read_csv(args.inpcsv)
samples = mappings["Sample_ID"].tolist()
run = args.runid
virus = args.virus
fasta_headers = args.fastaheaders

# sets up future list of collaborators
collablist = []

# unique list for collaborators in run
for i in samples:
    collab = i.split('_')[0]
    if collab not in collablist:
        collablist.append(collab)

# Database query setup
query = """
SELECT Sample_Systematic_ID, Sample_Name, Strain_Name, Flu_Type, Expected_Subtype, assembly_status, Assembly_quality, Total_reads, Uniq_mapped_read_percent, Coverage_10, Coverage_100, Viral_percent, Archaea_percent, Fungi_percent, Eukaryota_percent, Bacteria_percent, Variant_pos_sum_15pct, Sequencing_method, Sequencing_region
FROM tCEIRS_assemblies assemblies
JOIN tCEIRS_Extracts extracts ON assemblies.Extract_ID = extracts.Extract_ID
JOIN tCEIRS_Isolates isolates on extracts.Isolate_ID = isolates.Isolate_ID
WHERE assemblies.assembly_run = %s AND isolates.Flu_Type = %s AND extracts.Sample_Systematic_ID LIKE %s
"""

def add_file_to_tar(file_path, tar):
    if os.path.exists(file_path):
        tar.add(file_path)
        print(f"Added {file_path} to tarball")
    else:
        print(f"Skipped {file_path} as it does not exist")

try:
    # build connection to DB, turn off autocommit, set cursor
    db = mysql.connector.connect(host=host, user=user, password=password, database=database)
    db.autocommit = False
    cur = db.cursor()

    # loop through collab list
    for collab in collablist:
        name_of_file = f"{args.runid}_{collab}.tar.gz"
        tar_file = tarfile.open(name_of_file, "w:gz")
        
        with open(f"{run}_{collab}_assemblies.csv", 'w', newline='') as fp:
            myFile = csv.writer(fp)
            headers = ('Sample_Systematic_ID', 'Sample_Name', 'Strain_Name', 'Flu_Type', 'Expected_Subtype', 'assembly_status', 'Assembly_quality', 'Total_reads', 'Uniq_mapped_read_percent', 'Coverage_10', 'Coverage_100', 'Viral_percent', 'Archaea_percent', 'Fungi_percent', 'Eukaryota_percent', 'Bacteria_percent', 'Variant_pos_sum_15pct', 'Sequencing_method', 'Sequencing_region')
            myFile.writerow(headers)
            
            for s in samples:
                for v in virus:
                    # setup query variables
                    condition = (run, v, s + '%')
                    cur.execute(query, condition)
                    result = cur.fetchall()
                    myFile.writerows(result)
                    
                    # Add relevant files to the tarball
                    add_file_to_tar(f"{s}/02_assembly/{s}.fasta", tar_file)
                    add_file_to_tar(f"{s}/02_assembly/reads.1.fq.gz", tar_file)
                    add_file_to_tar(f"{s}/02_assembly/reads.2.fq.gz", tar_file)
                    add_file_to_tar(f"{s}/02_assembly/{s}_ref.bam", tar_file)
                    add_file_to_tar(f"{s}/02_assembly/{s}_ref.bam.bai", tar_file)
                    add_file_to_tar(f"{s}/03_qualityControl/{s}_quality_control.pdf", tar_file)
                    for header in fasta_headers:
                        add_file_to_tar(f"{s}/04_variants/{s}.{header}_variable_bases.tsv", tar_file)
                    fastqs = f"{s}/01_fastqs/"
                    if os.path.exists(fastqs):
                        for fastq in os.listdir(fastqs):
                            if fastq.endswith(".fastq.gz"):
                                add_file_to_tar(os.path.join(fastqs, fastq), tar_file)
                    else:
                        print(f"Skipped {fastqs} as it does not exist")

        tar_file.close()

# print db errors
except mysql.connector.Error as e:
    print("Database Error: {}".format(e))

except Exception as e:
    print("Error: {}".format(e))

# close connection
finally:
    if db.is_connected():
        cur.close()
        db.close()
