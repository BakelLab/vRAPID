#!/usr/bin/env python

"""
Written by: Adriana van de Guchte
Last Update: 4/10/2023
Version: 1.2
Purpose: Push viral assemblies to PDB

"""

# Dependencies 
import os, sys
import os.path
import shutil
import re
import yaml
import logging
import pandas as pd
import mysql.connector
import glob
from Bio import SeqIO

# ----------- FUNCTIONS --------------
# -------- Data Processing Functions

# function for pulling total and mapped reads from flagstat file
# note: assumes mapped reads occur in line with first mention of word "mapped"
def process_flagstat(file_path):
    return_dict = {}
    with open(file_path, 'r') as f:
        return_dict['total_reads'] = int(f.readline().split()[0])
        for line in f:
            if "mapped" in line:
                return_dict['mapped_reads'] = int(line.split()[0])
                break
    return return_dict

# function for handling the fasta files and assigning to appropirate location in data_dict
def process_fasta_file(file_path, segment_headers, multi,virus):
    return_dict = {'sequences': {}}
    if multi:
        with open(file_path) as f:
            # Initialize variables for the current segment header and sequence
            current_header = None
            current_sequence = ''

            # Loop through each line in the file
            for line in f:
                # Check if the line is a header line
                if line.startswith('>'):
                    # Check if the current sequence is not empty, and if so, add it to the correct segment in the 'data_dict'
                    if current_sequence:
                        for header in segment_headers:
                            if header in current_header:
                                return_dict['sequences'][header] = {'fasta': current_sequence,
                                                                    'type': 'Segment',
                                                                    'name': header}
                                break

                    # Set the current header to the new header line
                    current_header = line.strip()
                    # Reset the current sequence
                    current_sequence = ''
                else:
                    # Add the current line to the current sequence
                    current_sequence += line.strip()

            # After the loop is complete, add the final sequence to the correct segment in the 'return_dict'
            if current_sequence:
                for header in segment_headers:
                    if header in current_header:
                        return_dict['sequences'][header] = {'fasta': current_sequence,
                                                           'type': 'Segment',
                                                           'name': header}
                        break
    else:
        return_dict['sequences']['sequence'] = {'fasta': str(SeqIO.read(file_path, "fasta").seq.upper()),
                                                'type': 'Genome',
                                                'name': virus}
    return return_dict

# function for handling the protein sequences
def parse_protein(input_file):
    with open(input_file) as f:
        lines = f.readlines()
    sequences = {}
    header = ''
    sequence = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if sequence:
                sequences[header] = {'fasta': sequence, 'type': 'Protein', 'name': header}
                sequence = ''
            header = line[1:].split('|')[0]
        else:
            sequence += line
    if sequence:
        sequences[header] = {'fasta': sequence, 'type': 'Protein', 'name': header}
    return {'sequences': sequences}

 # function for handling the cds sequences
def parse_cds(input_file):
    with open(input_file) as f:
        lines = f.readlines()
    sequences = {}
    header = ''
    sequence = ''
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if sequence:
                sequences[header] = {'fasta': sequence, 'type': 'CDS', 'name': header}
                sequence = ''
            header = line[1:].split('|')[0]
        else:
            sequence += line
    if sequence:
        sequences[header] = {'fasta': sequence, 'type': 'CDS', 'name': header}
    return {'sequences': sequences}

# function for handling kraken report
def process_kraken(kraken_path):
    return_dict = {}
    return_dict['kraken_total'] = 0
    return_dict['bacteria_reads'] = 0
    return_dict['eukaryote_reads']=0
    return_dict['virus_reads']=0
    return_dict['fungi_reads']=0
    return_dict['archaea_reads']=0
    with open(kraken_path) as f:
        for line in f:
            read_count_clade = int(line.split()[1])
            classification = line.split()[5]
            if classification=='root':
                return_dict['kraken_total']=read_count_clade
            if classification=='Bacteria':
                return_dict['bacteria_reads']=read_count_clade
            if classification=='Eukaryota':
                return_dict['eukaryote_reads']=read_count_clade
            if classification=='Viruses':
                return_dict['virus_reads']=read_count_clade
            if classification=='Fungi':
                return_dict['fungi_reads']=read_count_clade
            if classification=='Archaea':
                return_dict['archaea_reads']=read_count_clade
    return return_dict

# function for assigning percentages
def process_percentages(kraken_dict,data_dict):
    return_dict = {}
    try:
            return_dict['bacteria_per'] = round((kraken_dict['bacteria_reads']/kraken_dict['kraken_total'])*100,2)
            return_dict['viral_per'] = round((kraken_dict['virus_reads']/kraken_dict['kraken_total'])*100,2)
            return_dict['eukaryote_per'] = round((kraken_dict['eukaryote_reads']/kraken_dict['kraken_total'])*100,2)
            return_dict['mapped_per'] = round((data_dict['flagstat']['mapped_reads']/data_dict['flagstat']['total_reads'])*100,2)
            return_dict['short_unmapped_per'] = round(100-return_dict['mapped_per'],2)
            return_dict['fungi_per'] = round((kraken_dict['fungi_reads']/kraken_dict['kraken_total'])*100,2)
            return_dict['archaea_per'] = round((kraken_dict['archaea_reads']/kraken_dict['kraken_total'])*100,2)
    except ZeroDivisionError:
            return_dict['bacteria_per'] = 0
            return_dict['viral_per'] = 0
            return_dict['eukaryote_per'] = 0
            return_dict['short_unmapped_per'] = 0
            return_dict['mapped_per'] = 0
            return_dict['fungi_per'] = 0
            return_dict['archaea_per'] = 0
    return return_dict
            
            
# function to handle variable bases tsv
def var_count(sample,multi,headers,variant_dir):
    var_count = 0
    if multi:
        for header in headers:
            file_path = os.path.join(variant_dir, sample+'.'+header + "_variable_bases.tsv").replace("\\", "/")
            with open(file_path) as f:
                f.readline()
                for line in f:
                    if line.split('\t')[2]=='FLAGGED':
                        var_count+=1
    else:
        file_pattern = os.path.join(variant_dir, f"*_variable_bases.tsv").replace("\\", "/")
        file_list = glob.glob(file_pattern)
        file_path = file_list[0]
        with open(file_path) as f:
                f.readline()
                for line in f:
                    if line.split('\t')[2]=='FLAGGED':
                        var_count+=1
    return var_count

# function for scoring (completeness, quality, status)
def assembly_scores(data_dict,scoring_dict):
    return_dict = {}
    virus = data_dict['virus']
    if virus == 'IAV' or virus == 'IBV':
        # Define a regular expression pattern to match the completeness information
        completeness_pattern = re.compile(r'nucleotide - (\w+(?: \w+)?)')

        # Define a dictionary to store the completeness information for each segment
        completeness_dict = {}
        
        features = data_dict['file_paths']['features_path']
#         opens the features table and finds the completeness assignment for each segment's nucleotide
        with open(features) as f:
            current_segment = None
            for line in f:
                if line.startswith(">Feature"):
                    # Extract the segment name from the header
                    current_segment = line.strip().split(' ')[1]
                else:
                    # Check if the line contains completeness information
                    completeness_match = completeness_pattern.search(line)
                    if completeness_match:
                        # Extract the completeness assignment from the line
                        completeness = completeness_match.group(1)

                        # Store the completeness information in the dictionary
                        if current_segment not in completeness_dict:
                            completeness_dict[current_segment] = []
                        completeness_dict[current_segment].append(completeness)
                        
#       sets up completeness scoring based on dictionary values
        total_score = 0
#       checks weight score for each segment identified in features table
        for key, value in completeness_dict.items():
            if key in scoring_dict[virus]['segment_weights']:
                if value[0] in scoring_dict[virus]['segment_weights'][key]:
#                   adds weight scores up
                    score = scoring_dict[virus]['segment_weights'][key][value[0]]
                    total_score += score
#         assign completeness score
        return_dict['completeness'] = total_score
        data_dict['scores']['completeness']=return_dict['completeness']
#         placeholder for IAV assembly status
        if total_score == 26:
            return_dict['assembly_status']='Complete'
        elif len(completeness_dict)==8:
            return_dict['assembly_status']='Partial'
        else:
            return_dict['assembly_status']='Failed'
      # check if target virus percent, completeness, and coverage is greater than score requirements
        if data_dict['scores']['percent'] >= scoring_dict[virus]['percent'] and data_dict['scores']['completeness'] >= scoring_dict[virus]['completeness'] and data_dict['scores']['coverage']['coverage_10'] >= scoring_dict[virus]['coverage_10']:
        #     check if variable bases is under score threshold
            if data_dict['scores']['variable_bases'] <= scoring_dict[virus]['variable_bases']:
                quality = 'Passed'
            else:
                quality = 'Check'
        else:
            quality = 'Failed'
        
# handling for non multisegmented viruses        
    else:
#         loops through headers so that subtype name can be ontained
        for header in data_dict['sequences']['primary']:
#         takes count of Ns in fasta to generate non ambiguous length for completeness calculation
            nCount = data_dict['sequences']['primary'][header]['fasta'].lower().count('n')
            seq_length = len(data_dict['sequences']['primary'][header]['fasta'])
            non_amb_length = seq_length - nCount
#             calculates completeness value
            return_dict['completeness']=round((non_amb_length/int(data_dict['length'])),2)
            data_dict['scores']['completeness']=return_dict['completeness']
#         handles completeness values above 1
            if return_dict['completeness']>1:
                return_dict['completeness']=1
#             assigns assembly status of complete or partial based on completeness value
            if return_dict['completeness']>=0.95:
                return_dict['assembly_status']='Complete'
            else:
                return_dict['assembly_status']='Partial'
#        starts loop for scoring of assembly quality
        # check if target virus percent, completeness, and coverage is greater than score requirements
        if data_dict['scores']['percent'] >= scoring_dict[virus]['percent'] and data_dict['scores']['completeness'] >= scoring_dict[virus]['completeness'] and data_dict['scores']['coverage']['coverage_10'] >= scoring_dict[virus]['coverage_10']:
        #     check if variable bases is under score threshold
            if data_dict['scores']['variable_bases'] <= scoring_dict[virus]['variable_bases']:
                quality = 'Passed'
            else:
                quality = 'Check'
        else:
            quality = 'Failed'
    return_dict['assembly_quality'] = quality
    return return_dict


# function for coverage calculation
def calc_coverage(multi,headers,qc_dir):
    # set base variable values 
    return_dict = {'coverage':{}}
    tot_num_lines_10 = 0
    tot_num_lines_100 = 0
    tot_total_lines = 0
    # for multisegmented viruses calc coverage across all segment coverage files
    if multi:
        for header in headers:
            file_path = os.path.join(qc_dir,header+"_coverage.txt").replace("\\", "/")
            with open(file_path,'r') as f:
                lines = f.readlines()
                num_lines_100 = sum(1 for line in lines if int(line.split()[2]) > 100)
                num_lines_10 = sum(1 for line in lines if int(line.split()[2]) > 10)
                total_lines = len(lines)
                tot_num_lines_10 = num_lines_10 + tot_num_lines_10
                tot_num_lines_100 = num_lines_100 + tot_num_lines_100
                tot_total_lines = total_lines + tot_total_lines
        return_dict['coverage']['coverage_100'] = round(tot_num_lines_100/tot_total_lines,2)
        return_dict['coverage']['coverage_10'] = round(tot_num_lines_10/tot_total_lines,2)
    # for non multisegmented viruses, calc coverage within coverage file
    else:
        file_pattern = os.path.join(qc_dir, f"*_coverage.txt")
        file_list = glob.glob(file_pattern)
        file_path = file_list[0]
        with open(file_path,'r') as f:
            lines = f.readlines()
            num_lines_100 = sum(1 for line in lines if int(line.split()[2]) > 100)
            num_lines_10 = sum(1 for line in lines if int(line.split()[2]) > 10)
            total_lines = len(lines)
            return_dict['coverage']['coverage_100'] = round(num_lines_100/total_lines,2)
            return_dict['coverage']['coverage_10'] = round(num_lines_10/total_lines,2)
            
    return return_dict

# -------- File Handling Functions

# merged all influenza segment variable bases tsv files into one file
def merge_tsvs(file_paths):
#     set up final dataframe
    dfs = pd.DataFrame()
    for fpath in file_paths:
        #     checks if segment file exists before merging
        if os.path.exists(fpath):
#             merges segment file into full df
            df = pd.read_csv(fpath, sep='\t')
            dfs = pd.concat([dfs, df], 
                  ignore_index = True)
    return dfs

# moves all files in files_list to directory for pdb access
def copy_files(files_list,output_dir):
    cwd = os.getcwd()
    output_path = os.path.join(cwd,output_dir)
#     clear new directory
    shutil.rmtree(output_dir, ignore_errors=True)
    os.mkdir(output_dir)
#     for each source,destination tuple set up paths
    for file_name, dest_path in files_list:
        dest_file = os.path.join(output_path,dest_path)
        source_file = os.path.join(cwd,file_name)
#         if the source path exists, perform the movement
        if os.path.exists(source_file):
            shutil.copy(source_file, dest_file)


# ------- Database Functions

# function to check database for previous assembly, will return either assembly_ID or insert note
def db_check(db,extract_id,run):
    cur = db.cursor()
    query = "Select assembly_id from tCEIRS_assemblies where Extract_ID='"+str(extract_id)+"' and assembly_run='"+str(run)+"'"
    cur.execute(query)
    result = cur.fetchall()
    
    # report out assembly id or flag for insertion
    if len(result)>0:
        result = result[0][0]
    else:
        result = 'insert'
    return result

# function to check database for previous assembly, will return either seq_ID or insert note
def db_check_seq(db,assembly_id,name,seq_type):
    cur = db.cursor()
    query = "SELECT `SequenceID` FROM `tCEIRS_assembly_sequences` WHERE `AssemblyID` = '"+str(assembly_id)+"' and `Sequence_name` = '"+str(name)+"' and `Sequence_type` = '"+str(seq_type)+"'"
    cur.execute(query)
    result = cur.fetchall()
    
    # report out assembly id or flag for insertion
    if len(result)>0:
        result = result[0][0]
    else:
        result = 'insert'
    return result

# assembly table insertion command
def insert_assembly_table(data_dict,db):
    # build database cursor
    cur = db.cursor()

    # set up base queries
    query = "INSERT INTO `tCEIRS_assemblies` (`Extract_ID`, `assembly_run`, `assembly_status`, `Total_reads`,`Uniq_mapped_read_percent`,`Short_unmapped_read_percent`,`completeness`,`Viral_percent`,`Archaea_percent`,`Fungi_percent`,`Eukaryota_percent`,`Bacteria_percent`,`Assembly_quality`,`Variant_pos_sum_15pct`, `Sequencing_method`, `Sequencing_region`, `assembly_subtype`, `Coverage_100`,`Coverage_10`) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)"

    # build value lists
    vals = (data_dict['extract_id'], data_dict['run'], data_dict['scores']['status'], data_dict['flagstat']['total_reads'],data_dict['percents']['mapped_per'],data_dict['percents']['short_unmapped_per'],data_dict['scores']['completeness'],data_dict['percents']['viral_per'],data_dict['percents']['archaea_per'],data_dict['percents']['fungi_per'],data_dict['percents']['eukaryote_per'],data_dict['percents']['bacteria_per'],data_dict['scores']['quality'],data_dict['scores']['variable_bases'], 'Illumina', 'Whole Genome', data_dict['subtype'],data_dict['scores']['coverage']['coverage_100'],data_dict['scores']['coverage']['coverage_10'])

    # execute queries
    cur.execute(query,vals)
    
    # get assembly id
    assembly_id=cur.lastrowid
    return assembly_id

# assembly table update command
def update_assembly_table(data_dict,db):
    # build database cursor
    cur = db.cursor()

    # set up base queries
    query = "UPDATE `tCEIRS_assemblies` SET `Extract_ID` = %s, `assembly_status` = %s, `Total_reads` = %s,`Uniq_mapped_read_percent` = %s,`Short_unmapped_read_percent` = %s,`completeness` = %s,`Viral_percent` = %s,`Archaea_percent` = %s,`Fungi_percent` = %s,`Eukaryota_percent` = %s,`Bacteria_percent` = %s,`Assembly_quality` = %s,`Variant_pos_sum_15pct` = %s, `Sequencing_method` = %s, `Sequencing_region` = %s, `assembly_subtype` = %s, `Coverage_100` = %s,`Coverage_10` = %s WHERE `assembly_ID` = %s"

    # build value lists
    vals = (data_dict['extract_id'], data_dict['scores']['status'], data_dict['flagstat']['total_reads'],data_dict['percents']['mapped_per'],data_dict['percents']['short_unmapped_per'],data_dict['scores']['completeness'],data_dict['percents']['viral_per'],data_dict['percents']['archaea_per'],data_dict['percents']['fungi_per'],data_dict['percents']['eukaryote_per'],data_dict['percents']['bacteria_per'],data_dict['scores']['quality'],data_dict['scores']['variable_bases'], 'Illumina', 'Whole Genome', data_dict['subtype'],data_dict['scores']['coverage']['coverage_100'],data_dict['scores']['coverage']['coverage_10'],data_dict['assembly_id'])

    # execute queries
    cur.execute(query,vals)

# sequences table insertion command
def insert_seq_table(assembly_id,fasta,name,seqtype,db):
    # build database cursor
    cur = db.cursor()

    # set up base query
    query = "INSERT INTO `tCEIRS_assembly_sequences` (`AssemblyID`, `Sequence_type`, `Sequence_name`, `Sequence`,`Sequence_length`) VALUES (%s,%s,%s,%s,%s)"
    
    # build value lists
    vals = (assembly_id, seqtype, name, fasta, len(fasta))

    # execute queries
    cur.execute(query,vals)
    
    # get seq id
    seq_id=cur.lastrowid
    return seq_id

# sequences table update command
def update_seq_table(seq_id,assembly_id,fasta,name,seqtype,db):
    # build database cursor
    cur = db.cursor()

    # set up base query
    query = "UPDATE `tCEIRS_assembly_sequences` SET `AssemblyID` = %s, `Sequence_type` = %s, `Sequence_name` = %s, `Sequence` = %s,`Sequence_length` = %s WHERE SequenceID = %s"
    
    # build value lists
    vals = (assembly_id, seqtype, name, fasta, len(fasta), seq_id)

    # execute queries
    cur.execute(query,vals)
    
# submissions table insertion command, no update command as that will be the responsibility of manual curation
def insert_sub_table(assembly_id,seq_id,quality,status,db):
    # build database cursor
    cur = db.cursor()
    
#     deletes any existing records in submission table
    try: 
        query = "DELETE FROM `tCEIRS_Submissions` WHERE `Assembly_ID` = '"+str(assembly_id)+"' AND `Seq_ID` = '"+str(seq_id)+"'"
        cur.execute(query)
        message = "Records deleted"
    except:
        message = "No records deleted"
    
    
#  only insert into table if assembly is passed/complete
    if quality == 'Passed' and status == 'Complete':
        # set up base query
        sub_query = "INSERT INTO `tCEIRS_Submissions` (`Assembly_ID`, `Seq_ID`) VALUES (%s,%s)"

        # build value lists
        sub_vals = (assembly_id, seq_id)

        # execute queries
        cur.execute(sub_query,sub_vals)
    return message




# ----------- Body
# database connection
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

# argument useage setup
if __name__ == "__main__":

    from argparse import ArgumentParser

    def usage(code="0"):
        print("error: " + str(code))
        print("\n\tusage: all-virus_assembly_push.py -s <sample> -r <run_id> -c <config>")
        sys.exit(0)

    # setup of arguments
    parser = ArgumentParser(description="Assembly push script for SARS-CoV-2, sCoV, MPX, IAV, and IBV")
    parser.add_argument('-s', '--ins', help='systematic sample name', required=True)
    parser.add_argument('-r', '--inr', help='run id', required=True)
    parser.add_argument('-c', '--inc', help='path to config file', required=True)
    parser.add_argument('-p', '--altpath', help='alternate file path for samples', required=False)

    args = parser.parse_args()

    # define user inputs
    sample = args.ins
    run = args.inr
    config = args.inc
    
    # read config.yaml file for inputs
    with open(config, 'r') as file:
        # Load the YAML file as a Python dictionary
        config_data = yaml.load(file, Loader=yaml.FullLoader)
    
    virus = config_data['virus']
    length = config_data['length']
    headers = config_data['ref_fasta_headers']
    base_dir = config_data['path']

    # overwrites base directory with alternate if one is provided as an argument
    if args.altpath:
        base_dir=args.altpath

    # begin logging
    logpath = os.path.join(sample,'05_status',sample+'.assembly-push.log')
    logging.basicConfig(filename=logpath, filemode='w', level=logging.INFO)

# reformat Influenza naming structure
    if virus == 'Influenza-A':
        virus = 'IAV'
        multi = True
        logging.info('virus set to IAV')
    elif virus == 'Influenza-B':
        virus = 'IBV'
        multi = True
        logging.info('virus set to IBV')
    else:
        multi = False
        logging.info('multi set to False')




    # variable/dictionary setups
    try:
        # build scoring dictionary
        scoring_dict = {
            'IAV':{
                'segment_weights': {
                    'HA': {
                        'complete':10,
                        'nearly complete':5
                    },
                    'NA': {
                        'complete':10,
                        'nearly complete':5
                    },
                    'NP': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'NS': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'PA': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'PB1': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'PB2': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'M': {
                        'complete':2,
                        'nearly complete':1
                    }
                },
                'variable_bases': 9,
                'percent': .4,
                'completeness': 20,
                'coverage_10':.7
            },
            'IBV':{
                'segment_weights': {
                    'HA': {
                        'complete':10,
                        'nearly complete':5
                    },
                    'NA': {
                        'complete':10,
                        'nearly complete':5
                    },
                    'NP': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'NS': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'PA': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'PB1': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'PB2': {
                        'complete':2,
                        'nearly complete':1
                    },
                    'M': {
                        'complete':2,
                        'nearly complete':1
                    }
                },
                'variable_bases': 9,
                'percent': .4,
                'completeness': 20,
                'coverage_10':.7
            },
            'SARS-CoV-2':{
                'variable_bases': 9,
                'percent': .4,
                'completeness':.7,
                'coverage_10':.7
            },
            'sCoV':{
                'variable_bases': 9,
                'percent': .4,
                'completeness':.7,
                'coverage_10':.7
            },
            'MPX':{
                'variable_bases': 9,
                'percent': .4,
                'completeness':.7,
                'coverage_10':.7
            }
        }
        logging.info('Scoring dictionary structure built successfully.') 
    except Exception as e:
        logging.exception(f'Error building scoring dictionary structure: %s', e)


    try:
        # build base dictionary structure
        data_dict = {
            'run': run,
            'sample': sample,
            'extract_id': sample.split('_')[1],
            'assembly_id': None,
            'virus': virus,
            'length': length,
            'subtype': None,
            'file_paths': {
                'base_dir': base_dir,
                'assembly_dir': os.path.join(base_dir, sample, '02_assembly').replace("\\", "/"),
                'qc_dir': os.path.join(base_dir, sample, '03_qualityControl').replace("\\", "/"),
                'variant_dir': os.path.join(base_dir, sample, '04_variants').replace("\\", "/"),
                'fasta_path': os.path.join(sample + '.fasta').replace("\\", "/"),
                'flagstat_path': os.path.join(sample + '_refbam.flagstat').replace("\\", "/"),
                'kraken_path': os.path.join(sample + '_kraken_report.out').replace("\\", "/")
            },
            'flagstat': {
                'total_reads': 0,
                'mapped_reads': 0,
            },
            'sequences': {
                'primary':{}
            },
            'percents': {
                'eukaryote_per': 0,
                'bacteria_per': 0,
                'viral_per': 0,
                'fungi_per': 0,
                'archaea_per': 0,
                'short_unmapped_per':0,
                'mapped_per': 0
            },
            'scores':{
                'completeness':0,
                'quality': None,
                'status': None,
                'variable_bases':0,
                'percent': 0,
                'coverage': 5
            }
        }
        logging.info('Primary dictionary structure built successfully.') 
    except Exception as e:
        logging.exception(f'Error building primary dictionary structure: %s', e)
        

    try:
        # set up structure for primary sequence information depending on header information
        # if virus has segments sets up segment dictionary
        if multi:
            for header in headers:
                data_dict['sequences']['primary'][header] = {
                    'header': header,
                    'type': None,
                    'name': None,
                    'fasta': None}
            logging.info('Segment dictionary and headers built successfully.')
        # if virus does not have segments sets primary dictionary
        else:
            data_dict['sequences']['primary']['sequence'] = {
                'header': None,
                'type': None,
                'name': virus,
                'fasta': None}
            logging.info('Sequence dictionary built assigned successfully.')
    except Exception as e:
        logging.exception(f'Error building primary sequence dictionary: %s', e)

# ------ Primary Data Processing
    # Body

    # assign full file paths
    try:
        data_dict['file_paths']['fasta_path']=os.path.join(data_dict['file_paths']['assembly_dir'],data_dict['file_paths']['fasta_path']).replace("\\", "/")
        data_dict['file_paths']['flagstat_path']=os.path.join(data_dict['file_paths']['qc_dir'],data_dict['file_paths']['flagstat_path']).replace("\\", "/")
        data_dict['file_paths']['kraken_path']=os.path.join(data_dict['file_paths']['qc_dir'],data_dict['file_paths']['kraken_path']).replace("\\", "/")
        logging.info('Fasta, flagstat, and kraken paths assigned successfully.')
    except Exception as e:
        logging.exception(f'Error assigning fasta, flagstat, and kraken paths: %s', e)
        
    # assign Influenza specific file paths
    if virus=='IAV' or virus=='IBV':
        try:
            data_dict['file_paths']['cds_path']=os.path.join(data_dict['file_paths']['assembly_dir'],sample+".features_cds.fa").replace("\\", "/")
            data_dict['file_paths']['protein_path']=os.path.join(data_dict['file_paths']['assembly_dir'],sample+".features_protein.fa").replace("\\", "/")
            data_dict['file_paths']['features_path']=os.path.join(data_dict['file_paths']['assembly_dir'],sample+".features_table.txt").replace("\\", "/")
            data_dict['file_paths']['subtype_path']=os.path.join(data_dict['file_paths']['assembly_dir'],sample+"_subtype.txt").replace("\\", "/")
            logging.info('Influenza specific file paths assigned successfully.')
        except Exception as e:
            logging.exception(f'Error assigning Influenza specific file paths: %s', e)
            

    # assigns total and mapped reads to data_dict
    try:
        data_dict['flagstat']['total_reads'] = process_flagstat(data_dict['file_paths']['flagstat_path'])['total_reads']
        data_dict['flagstat']['mapped_reads'] = process_flagstat(data_dict['file_paths']['flagstat_path'])['mapped_reads']
        logging.info('Flagstat parsing run successfully.')
    except Exception as e:
        logging.exception(f'Error processing flagstat file: %s', e)
        

    # process fasta files
    try:
        fastas = process_fasta_file(data_dict['file_paths']['fasta_path'],headers,multi,virus)
        logging.info('Fasta processing run successfully.')
    except Exception as e:
        logging.exception(f'Error processing fasta file: %s', e)
        
        
    # assign processed fastas to appropriate dictionary key
    try:
        for segment in fastas['sequences']:
            if segment in data_dict['sequences']['primary']:
                data_dict['sequences']['primary'][segment]  = fastas['sequences'][segment]
                logging.info(f'Fasta assignment run successfully.')
    except Exception as e:
        logging.exception(f'Error in fasta assignment: %s', e)

    # assign processed cds and protein segments to appropriate dict location
    if virus == 'IAV' or virus == 'IBV':
        try:
            data_dict['sequences']['CDS'] = {}
            data_dict['sequences']['Protein'] = {}
            cds = parse_cds(data_dict['file_paths']['cds_path'])
            protein = parse_protein(data_dict['file_paths']['protein_path'])
        #     assign cds and proteins to appropriate headers    
            for header in cds['sequences']:
                data_dict['sequences']['CDS'][header]  = cds['sequences'][header]
                logging.info(f'{header} CDS information parsed and assigned successfully.')
            for header in protein['sequences']:
                data_dict['sequences']['Protein'][header]  = protein['sequences'][header]
                logging.info(f'{header} protein infomation parsed and assigned successfully.')
        except Exception as e:
            logging.exception(f'Error in Influenza protein or CDS parsing: %s', e)
                
                
    # assign subtype
    if virus == 'sCoV':
        try:
            data_dict['subtype'] = headers
            logging.info('sCoV subtype assigned successfully.')
        except Exception as e:
            logging.exception(f'Error assining sCoV subtype: %s', e)
            
    elif virus == 'IAV' or virus == 'IBV':
        try:
            with open(data_dict['file_paths']['subtype_path'],'r') as f:
                data_dict['subtype'] = f.readline().split()[1]
                logging.info('Influenza subtype assigned successfully.')
        except Exception as e:
            logging.exception(f'Error assigning Influenza subtype: %s', e)
            
    else:
        data_dict['subtype'] = ''
        
        
    # process kraken report and assign percent scores
    try:
        kraken_report = process_kraken(data_dict['file_paths']['kraken_path'])
        percents = process_percentages(kraken_report,data_dict)
        data_dict['percents']['eukaryote_per'] = percents['eukaryote_per']
        data_dict['percents']['bacteria_per'] = percents['bacteria_per']
        data_dict['percents']['viral_per'] = percents['viral_per']
        data_dict['percents']['short_unmapped_per'] = percents['short_unmapped_per']
        data_dict['percents']['mapped_per'] = percents['mapped_per']
        data_dict['percents']['archaea_per'] = percents['archaea_per']
        data_dict['percents']['fungi_per'] = percents['fungi_per']
        data_dict['scores']['percent'] = data_dict['percents']['viral_per']
        logging.info('Kraken report parsed successfully.')
    except Exception as e:
        logging.exception('Error in kraken parsing: %s', e)


    # assigns variables bases count to data dict
    try:
        data_dict['scores']['variable_bases'] = var_count(sample,multi,headers,data_dict['file_paths']['variant_dir'])
        logging.info('Variable bases count parsed successfully.')
    except Exception as e:
        logging.exception(f'Error in variable base count parsing: %s', e)
        

    # assign genome coverage at or above 100
    try:
        coverage = calc_coverage(multi,headers,data_dict['file_paths']['qc_dir'])
        data_dict['scores']['coverage'] = coverage['coverage']
        logging.info('Coverage function parsed successfully.')
    except Exception as e:
        logging.exception(f'Error in coverage parsing: %s', e)


    # generate final assembly scores (completeness, status, quality)
    try:
        final_scores = assembly_scores(data_dict,scoring_dict)
        logging.info('Final scoring function run successfully.')
    except Exception as e:
        logging.exception(f'Error in scoring parsing: %s', e)

    # assign final scores to data_dict
    try:
        data_dict['scores']['completeness'] = final_scores['completeness']
        data_dict['scores']['status'] = final_scores['assembly_status']
        data_dict['scores']['quality'] = final_scores['assembly_quality']
        # in case of influenza replace completeness weighted score with coverage_10 percent
        if virus=='IAV' or virus=='IBV':
            data_dict['scores']['completeness'] = data_dict['scores']['coverage']['coverage_10']
        logging.info('Assembly scores assigned successfully.')
    except Exception as e:
        logging.exception(f'Error assigning assembly scores: %s', e)
        
        
        
    # removes primary sequence header if no data generated for it (handles missing influenza segments)
    if virus=='IAV' or virus=='IBV':
        for header in headers:
            try:
                # Iterate over the nested keys of each section
                if data_dict['sequences']['primary'][header]['fasta'] is None or data_dict['sequences']['primary'][header]['fasta'] =="":
                    try:
                        del data_dict['sequences']['primary'][header]
                        logging.info(f'{header} fasta empty and removed successfully.')
                    except Exception as e:
                        logging.exception(f'Error removing empty fasta for {header}: %s', e)
                logging.info(f'{header} checked for empty fasta.')
            except Exception as e:
                        logging.exception(f'Error assessing {header} fasta: %s', e)        


    # merge all flu variable base tsvs into one file
    if virus == 'IAV' or virus == 'IBV':
        try:
            file_list = []
        #     builds list of potential variable base files based on segment header names
            for header in headers:
                file_path = os.path.join(data_dict['file_paths']['variant_dir'], sample + '.' + header + "_variable_bases.tsv").replace("\\", "/")
                file_list.append(file_path)
        #     merges all existing segment variable base files
            merged_df = merge_tsvs(file_list)
        #     reads merged dataframe into new tsv
            merged_df.to_csv(os.path.join(data_dict['file_paths']['variant_dir'], sample+'_variable_bases.tsv'), sep="\t",index=False)
            logging.info(f'Influenza merged variable base file built successfully.')
        except Exception as e:
            logging.exception(f'Error building merged influenza variable base file: %s', e)


# ----- Begin Transaction
    try:
        # setup DB connection
        db = mysql.connector.connect(host=host,user=user,password=password,database=database)
        db.autocommit = False
        if db.is_connected():
            logging.info('Database connection opened.')
            
        # checks if sample/run exists in DB already to determine if running insertions or updates 
        if db_check(db,data_dict['extract_id'],data_dict['run']) == 'insert':
        #     if assembly does not exist, inserts all new records
            data_dict['assembly_id'] = insert_assembly_table(data_dict,db)
            # loop through all potential sequences to upload to db for sequence and submission table
            for cat in data_dict['sequences']:
                for header in data_dict['sequences'][cat]:
        #         insert new sequence info and collect newly generated seq id
                    seq_id = insert_seq_table(data_dict['assembly_id'],data_dict['sequences'][cat][header]['fasta'],data_dict['sequences'][cat][header]['name'],data_dict['sequences'][cat][header]['type'],db)
        #             insert new submissions for complete/passed sequences
                    insert_sub_table(data_dict['assembly_id'],seq_id,data_dict['scores']['quality'],data_dict['scores']['status'],db)
        else:
        #     if assembly does exist, collects existing assembly id
            data_dict['assembly_id'] = db_check(db,data_dict['extract_id'],data_dict['run'])
            update_assembly_table(data_dict,db)
            for cat in data_dict['sequences']:
                for header in data_dict['sequences'][cat]:
                    seq_id = db_check_seq(db,data_dict['assembly_id'],data_dict['sequences'][cat][header]['name'],data_dict['sequences'][cat][header]['type'])
                    if seq_id == "insert":
        #                 insert new sequence info and collect newly generated seq id
                        seq_id = insert_seq_table(data_dict['assembly_id'],data_dict['sequences'][cat][header]['fasta'],data_dict['sequences'][cat][header]['name'],data_dict['sequences'][cat][header]['type'],db)
            #             insert new submissions for complete/passed sequences
                        insert_sub_table(data_dict['assembly_id'],seq_id,data_dict['scores']['quality'],data_dict['scores']['status'],db)
                    else:
        #                 update exisiting sequence table records
                        update_seq_table(seq_id,data_dict['assembly_id'],data_dict['sequences'][cat][header]['fasta'],data_dict['sequences'][cat][header]['name'],data_dict['sequences'][cat][header]['type'],db)
        #                 reinsert records into submissions table for complete/passed sequences
                        insert_sub_table(data_dict['assembly_id'],seq_id,data_dict['scores']['quality'],data_dict['scores']['status'],db)
    # Primary rollback if error is encountered
    except Exception as e:
        db.rollback()
        logging.exception(f'Database upload error, database rolled back: %s', e)
        
    # If no errors, commit and move files into destination folder
    else:
        db.commit()
        logging.info('Database upload success.')
        # file movement for final links on PDB
        # setup paths and sources used by all viruses
        try:
            # set output directory location
            output_dir = '/sc/arion/projects/CRIPT/crip_surveillance/www/crip-final-assemblies/'+str(data_dict['assembly_id'])
            #     fasta paths
            fasta_source = data_dict['file_paths']['fasta_path']
            fasta_dest = os.path.join('', 'aid_' + str(data_dict['assembly_id']) + '.fa')
            fasta_dest2 = os.path.join('', str(sample) + '_final.fa')

            #     qc paths
            qc_source = os.path.join(data_dict['file_paths']['qc_dir'], sample + '_quality_control.pdf')
            qc_dest = os.path.join('', str(sample) + '_final.report.pdf')

            #     variant paths
            if not multi:
                # non multisegmented variant source
                variant_source = os.path.join(data_dict['file_paths']['variant_dir'], sample+'.'+config_data['ref_fasta_headers']+'_variable_bases.tsv')
            
            # variant destination for all viruses
            variant_dest = os.path.join('', str(sample) + '_final.variants.calls.txt')

            # set up file movement lists depending on virus type
            if virus == 'IAV' or virus == 'IBV':
                #     define source and destination names
                #     variant path for influenza
                variant_source = os.path.join(data_dict['file_paths']['variant_dir'], sample+'_variable_bases.tsv')

                #     fasta paths
                fasta_dest2 = os.path.join('', str(sample) + '_final.fa')

                #     cds features
                cds_source = os.path.join(data_dict['file_paths']['assembly_dir'], sample + '.features_cds.fa')
                cds_dest = os.path.join('', str(sample) + '_final.features_cds.fa')

                #     protein features
                protein_source = os.path.join(data_dict['file_paths']['assembly_dir'], sample + '.features_protein.fa')
                protein_dest = os.path.join('', str(sample) + '_final.features_protein.fa')

                #     features table
                features_source = os.path.join(data_dict['file_paths']['assembly_dir'], sample + '.features_table.txt')
                features_dest = os.path.join('', str(sample) + '_final.features_table.txt')

                #     subtype
                subtype_source = os.path.join(data_dict['file_paths']['assembly_dir'], sample + '_subtype.txt')
                subtype_dest = os.path.join('', str(sample) + '_final_subtype.txt')

                #     setup list of tuples for function input
                file_moves = [(fasta_source, fasta_dest), (fasta_source, fasta_dest2), (qc_source, qc_dest),
                            (variant_source, variant_dest), (cds_source, cds_dest), (protein_source, protein_dest),
                            (features_source, features_dest), (subtype_source, subtype_dest)]
            else:
                #     setup list of tuples for function input
                file_moves = [(fasta_source, fasta_dest), (qc_source, qc_dest), (variant_source, variant_dest)]
            logging.info('File movements set up successfully.')

        except Exception as e:
            logging.exception(f'Error setting up file movement paths: %s', e)
            
        # run file movement function to finally transfer all files into PDB hosting location
        try:
            copy_files(file_moves,output_dir)
            logging.info('File copy to PDB ran successfully.')
        except Exception as e:
            logging.exception(f'Error in file copy to PDB: %s', e)

    # Close database connection
    finally:
        db.rollback()
        if db.is_connected():
            db.close()
            logging.info('Database connection closed.')

# system exit for improper arguments
if len(sys.argv)!=7|len(sys.argv)!=9:
   print(usage())
   sys.exit(0)
            