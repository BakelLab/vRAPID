#!/usr/bin/env python

"""
Written by: Adriana van de Guchte
Modified by: Zain
Last Update: 01/24/2023
Version: 1.0
Purpose: Push viral assemblies to PDB

"""


# TO DO: figure out common errors/issues in assembly files, which need bypasses and which need error reports

import os, sys
import os.path
import mysql.connector
from Bio import SeqIO
import shutil

# function to check database for previous assembly, will return either assembly_ID or insert note
def db_check(db,extract_id,run):
    cur = db.cursor()
    query = "Select assembly_id from tCEIRS_assemblies where Extract_ID='"+str(extract_id)+"' and assembly_run='"+str(run)+"'"
    cur.execute(query)
    result = cur.fetchall()
    # report out assembly id or flag for insertion
    if len(result)>0:
#         for row in data:
        result = result[0][0]
    else:
        result = 'insert'
    return result

# setup to collect all config info for database connection
path = os.path.expanduser('~') + '/.my.cnf.pdbrw'
with open(path) as cnf_file:
    for line in cnf_file:
        if line.startswith('user='):
            user = line.rstrip()[5:]
        if line.startswith('password='):
            password = line.rstrip()[9:]
        if line.startswith('host='):
            host = line.rstrip()[5:]
        if line.startswith('database='):
            database = line.rstrip()[9:]

if __name__ == "__main__":

    from argparse import ArgumentParser

    def usage(code="0"):
        print("error: " + str(code))
        print("\n\tusage: assembly_push_pdb.py -s <sample> -r <run_id> -d <directory> -p <path> -l <length> -v <virus>")
        sys.exit(0)

    # setup of arguments
    parser = ArgumentParser(description="Parse lineage assignments and versions from lineage tables")
    parser.add_argument('-s', '--ins', help='systematic sample name', required=True)
    parser.add_argument('-r', '--inr', help='run ID', required=True)
    parser.add_argument('-d', '--ind', help='directory', required=True, default='PVI')
    parser.add_argument('-p', '--path', help='path to base directory', required=True)
    parser.add_argument('-l', '--length', help='virus length', required=True, default=float('29870.0'))
    parser.add_argument('-v', '--virus', help='virus name in GenBank File', required=True, default=str('MN908947'))
    args = parser.parse_args()

    # define user inputs
    sample = args.ins
    virus = args.virus
    run = args.inr
    if args.ind=='CHARM':
        base_dir='/sc/arion/projects/CHARM/genomes/assembly'
    else:
        #base_dir='/sc/arion/projects/PVI/'+args.path
        base_dir=(args.path)
    # system exits for pools and negative controls
    extract_id=sample.split('_')[1]
    colloborator_id=sample.split('_')[0]
    if colloborator_id=='Pool':
        sys.exit(0)
    if extract_id.lower()=='neg':
        sys.exit(0)

    # set path to sample QC and variant directories
    sample_path_qc=base_dir+'/'+sample+'/03_qualityControl'
    sample_path_variants=base_dir+'/'+sample+'/04_variants'

    # build base variables
    bacteria_reads=0
    eukaryote_reads=0
    coronavirus_reads=0
    iav_reads=0
    ibv_reads=0
    var15_count=0
    mpx_reads=0

    method = 'Illumina'
    region = 'Whole Genome'

    # set up historic "old run list"
    old_run_list=['TD01619_TD01622_TD01623','TD01627','TD01628_TD01629','TD01630_TD01631','TD01635','TD01637_TD01638',
    'TD01640','TD01640_TD01644_TD01654','TD01640_TD01644_TD01654_TD01753','TD01644','TD01646','TD01654','TD01655','TD01662','TD01686',
    'TD01694','TD01697','TD01704','TD01719','TD01724','TD01735','TD01748','TD01753','TD01778','TD01785','TD01809','TD01833','TD01836','TD01836_rerun','TD01836_TD01836_rerun',
    'TD01866','TD01866_rerun','TD01866_TD01866_rerun','TD01889','TD01916','TD01924']


    try:
        # build connection to DB, turn off autocommit, set cursor
        db = mysql.connector.connect(host=host,user=user,password=password,database=database)
        db.autocommit = False
        cur = db.cursor()

        # checks database for previous entries
        assembly_id = db_check(db,extract_id,run)

        # get total and mapped reads from refbam
        with open("%s/%s_refbam.flagstat" % (sample_path_qc, sample)) as f:
            total_reads = int(f.readline().split()[0])
            f.readline()
            f.readline()
            f.readline()
            mapped_reads = int(f.readline().split()[0])
        # get kraken classification reads
        kraken= sample_path_qc+"/"+str(sample)+"_kraken_report.out"
        print(kraken)
        with open(kraken) as f:
            for line in f:
                read_count_clade = int(line.split()[1])
                classification = line.split()[5]
                if classification=='Bacteria':
                    bacteria_reads=read_count_clade
                if classification=='Eukaryota':
                    eukaryote_reads=read_count_clade
                if classification=='Coronaviridae':
                    coronavirus_reads=read_count_clade
                if classification=='Influenza A virus':
                    iav_reads=read_count_clade
                if classification=='Influenza B virus':
                    ibv_reads=read_count_clade
                if classification=='Monkeypox virus':
                	mpx_reads=read_count_clade

        # get variant calls
        with open("%s/%s_variable_bases.tsv" % (sample_path_variants, virus)) as f:
            f.readline()
            for line in f:
                if line.split('\t')[2]=='FLAGGED':
                    var15_count+=1

        # alert for no coronavirus reads
        if float(args.length) == float('29870.0'):
        	if coronavirus_reads == 0:
        		print(sample+" no coronavirus reads")
        	else:
        		coronavirus_reads=coronavirus_reads+mapped_reads
        # alert for no MPX reads
        if float(args.length) == float('197205.0'):
        	if mpx_reads == 0:
        		print(sample+ " no MPX reads")
        	else:
        		mpx_reads=mpx_reads+mapped_reads
        		
        # set up all percentages
        unique_mapped_per=round((mapped_reads/total_reads)*100,2)
        short_unmapped_per=round((100-unique_mapped_per),2)
        bact_per=round((bacteria_reads/total_reads)*100,2)
        euk_per=round((eukaryote_reads/total_reads)*100,2)
        cor_per=round((coronavirus_reads/total_reads)*100,2)
        iav_per=round((iav_reads/total_reads)*100,2)
        ibv_per=round((ibv_reads/total_reads)*100,2)
        mpx_per=round((mpx_reads/total_reads)*100,2)

        # setup FASTA seq and info
        fasta_path=base_dir+'/'+sample+'/02_assembly/'+sample+".fasta"
        fastaSequence=SeqIO.read(fasta_path, "fasta")
        nCount = fastaSequence.seq.lower().count('n')
        length = len(fastaSequence.seq)
        non_amb_length=length-nCount

        # setup genome completeness and assembly scores
        completeness=round((non_amb_length/float(args.length)),2)

        if completeness>1:
            completeness=1
        seq=fastaSequence.seq.upper()
        if completeness>=0.95:
            assembly_status='Complete'
        else:
            assembly_status='Partial'

        # quality scoring
        if cor_per<4:
            assembly_quality='Failed'

        else:
            if run in old_run_list:
                if var15_count>4:
                    assembly_quality='Check'
                else:
                    assembly_quality='Passed'
            else:
                if var15_count>9:
                    assembly_quality='Check'
                else:
                    assembly_quality='Passed'

        # base queries
        update = "UPDATE `tCEIRS_assemblies` SET `Extract_ID` = %s, `assembly_run` = %s, `assembly_status` = %s,`Total_reads` = %s,`Uniq_mapped_read_percent` = %s,`Short_unmapped_read_percent` = %s,`completeness` = %s,`IAV_percent` = %s,`IBV_percent` = %s,`coronavirus_percent` = %s,`Eukaryota_percent` = %s,`Bacteria_percent` = %s,`Assembly_quality` = %s,`Variant_pos_sum_15pct` = %s, `Sequencing_method` = %s, `Sequencing_region` = %s where `assembly_ID`= %s"
        update_seq = "UPDATE `tCEIRS_assembly_sequences` SET `Sequence_type` = 'Genome', `Sequence_name` = 'SARS-CoV-2', `Sequence` = %s,`Sequence_length` = %s,`Sequence_quality` = %s,`Variant_pos_sum` = %s where `AssemblyID`= %s"
        insert = "INSERT INTO `tCEIRS_assemblies` (`Extract_ID`, `assembly_run`, `assembly_status`, `Total_reads`,`Uniq_mapped_read_percent`,`Short_unmapped_read_percent`,`completeness`,`IAV_percent`,`IBV_percent`,`coronavirus_percent`,`Eukaryota_percent`,`Bacteria_percent`,`Assembly_quality`,`Variant_pos_sum_15pct`, `Sequencing_method`, `Sequencing_region`) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"
        insert_seq ="INSERT INTO `tCEIRS_assembly_sequences`(`AssemblyID`, `Sequence_type`, `Sequence_name`, `Sequence`,`Sequence_length`,`Sequence_quality`,`Variant_pos_sum`) VALUES (%s, 'Genome', 'SARS-CoV-2', %s,%s,%s,%s)"

        # value creation and sql execution based on db check for previous entries
        if assembly_id == 'insert':
            vals = (extract_id,run,assembly_status,total_reads,unique_mapped_per,short_unmapped_per,completeness,iav_per,ibv_per,cor_per,euk_per,bact_per,assembly_quality,var15_count,method,region)
            cur.execute(insert, vals)
            assembly_id=cur.lastrowid
            seq_vals = (assembly_id,str(seq),length,assembly_quality,var15_count)
            cur.execute(insert_seq, seq_vals)
        else:
            vals = (extract_id,run,assembly_status,total_reads,unique_mapped_per,short_unmapped_per,completeness,iav_per,ibv_per,cor_per,euk_per,bact_per,assembly_quality,var15_count,method,region,assembly_id)
            cur.execute(update, vals)
            seq_vals = (str(seq),length,assembly_quality,var15_count,assembly_id)
            cur.execute(update_seq, seq_vals)

        # file path setup for PDB
        output_dir='/sc/arion/projects/CRIPT/crip_surveillance/www/crip-final-assemblies/'+str(assembly_id)

        # clear new directory
        shutil.rmtree(output_dir, ignore_errors=True)
        os.mkdir(output_dir)
        # copy reports into new file path
        shutil.copy(fasta_path, output_dir+'/aid_'+str(assembly_id)+'.fa')
        shutil.copy(sample_path_qc+'/'+sample+'_quality_control.pdf', output_dir+'/'+str(sample)+'_final.report.pdf')
        shutil.copy(sample_path_variants+'/'+virus+'_variable_bases.tsv', output_dir+'/'+str(sample)+'_final.variants.calls.txt')
        
    # exceptions for rollbacks
    except IOError as e:
        print(sample+" failed parsing, missing file: {}".format(e))
        # reverting changes because of exception
        db.rollback()
    except mysql.connector.Error as e:
        print(sample+" failed to update record to database rollback: {}".format(e))
        # reverting changes because of exception
        db.rollback()
    except Exception as e:
        print(sample+" failed: {}".format(e))
        # reverting changes because of exception
        db.rollback()

    # if no exceptions raises, commid to db and report success
    else:
        db.commit()
        print(sample+" committed successfully")

    # closing database connection
    finally:
        if db.is_connected():
            cur.close()
            db.close()

if len(sys.argv)!=8 | len(sys.argv)!=10:
   print(usage())
   sys.exit(0)
