#!/usr/bin/env python

"""
Written by: Adriana van de Guchte
Last Update: 4/4/2023
Version: 1.1
Purpose: Move MultiQC and BamQC files from Minerva into PDB holding dir and setup folder system.
Last Update: TD# now an arg.

"""

# required dependencies 
import os, sys
import os.path
import mysql.connector
import shutil

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

# build database connection and cursor
db = mysql.connector.connect(host=host,user=user,password=password,database=database)
cur = db.cursor()


if __name__ == "__main__":

    # # collect TD# from directory name
    directory = os.getcwd()

    # find corresponding xt number
    # split directory by / delimiter and take [-1] item as the TD#
    td = snakemake.config["run_id"]
    # build and execute query to collect XT number from DB based on TD#
    query = "SELECT `Sequence_Plate_ID` FROM `tIlluminaCoreSubmissions` WHERE `Sequence_Run_ID` = %s"
    cur.execute(query, td)
    xt = cur.fetchone()[0]  
    
    # Source path
    BAMQCsource = directory + "/multi_bamqc"
    MultiQCsource = directory + "/multiqc_report.html"

    # Destination path
    # file path setup for PDB
    destination = '/sc/arion/projects/CRIPT/crip_surveillance/www/crip-run-qc'

    # Destination names
    BAMQCdest = destination + "/" + xt
    MultiQCdest = destination + "/" + xt + "/mutliqc_report.html"
    
    # Copy the content of source to destination
    
    try:
        shutil.rmtree(BAMQCdest, ignore_errors=True)
        shutil.copytree(BAMQCsource, BAMQCdest)
        shutil.copy(MultiQCsource, MultiQCdest)
        f = open("file_movement_message.txt", "w")
        f.write("Files copied successfully.")
        f.close()
    
    # If source and destination are same
    except shutil.SameFileError:
        f = open("file_movement_message.txt", "w")
        f.write("Source and destination represents the same file.")
        f.close()
    
    # If there is any permission issue
    except PermissionError:
        f = open("file_movement_message.txt", "w")
        f.write("Permission denied.")
        f.close()
    
    # For other errors
    except Exception as e:
        f = open("file_movement_message.txt", "w")
        f.write("Error occurred while copying file. {}".format(e))
        f.close()

# system exit for improper arguments
if len(sys.argv)!=3:
   print(usage())
   sys.exit(0)
