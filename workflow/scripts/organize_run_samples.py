import glob
import os
from datetime import datetime
import mysql.connector
import re
import shutil

# Get current datestamps
sYearMonth = datetime.now().strftime('%Y-%m')
sCurrentDate = datetime.now().strftime('%Y-%m-%d')
sCurrentTime = datetime.now().strftime('%H-%M-%S')

# Set assembly type
#sAsmType = "asm=spades" if flSpades else ""

# Open database connection and begin transaction
path = os.path.expanduser('~') + '/.my.cnf.pdbrw'
dbPUB = mysql.connector.connect(option_files=path, option_groups="vanbah01_pathogens", database="vanbah01_pathogens")
dbPUB.autocommit = False

# Get a look-up table of SpecimenId to accessionNumber from eRAP tPVI_Surveillance
hIDlookup = {}
oIDlookup = dbPUB.cursor()
oIDlookup.execute("SELECT E.Sample_Systematic_ID,I.Sample_Name,I.Flu_Type,I.Expected_Subtype FROM tCEIRS_Extracts as E JOIN tCEIRS_Isolates as I ON E.Isolate_ID=I.Isolate_ID;")
for rRow in oIDlookup.fetchall():
       sSampleID, sSampleName, sExpectedType, sExpectedSubtype = rRow
       hIDlookup[sSampleID] = {
               "ExpectedType": sExpectedType,
               "ExpectedSubtype": sExpectedSubtype,
               "Name": sSampleName
               }


# Start processing folder list
hCheckPassed = {}
with os.scandir('./') as dir:
    for entry in dir:
        if entry.is_dir() and not entry.name.startswith('.'):
            R1 = glob.glob(f"{entry.name}/{entry.name}/{entry.name}_[A-Z]*_R1_001.fastq.gz")[0]
            R2 = glob.glob(f"{entry.name}/{entry.name}/{entry.name}_[A-Z]*_R2_001.fastq.gz")[0]
            if os.path.exists(R1) and os.path.exists(R2):
                hCheckPassed[entry.name] = 1
            else:
                print(f"Warning: Did not find paired fastq.gz files in folder '{entry.name}', skipping")

# For the folders that pass all checks, move things in place
for sInputFolder in hCheckPassed.keys():
    if re.match('^[A-Z]{2}_\d+$', sInputFolder):
        sPrefix, nExtractID = sInputFolder.split('_')
        ExpectedType = hIDlookup[sInputFolder]['ExpectedType'] if sInputFolder in hIDlookup and 'ExpectedType' in hIDlookup[sInputFolder] else 'Unknown'
        ExpectedSubtype = hIDlookup[sInputFolder]['ExpectedSubtype'] if sInputFolder in hIDlookup and 'ExpectedSubtype' in hIDlookup[sInputFolder] else 'Unknown'
        sSampleName = hIDlookup[sInputFolder]['Name'] if sInputFolder in hIDlookup and 'Name' in hIDlookup[sInputFolder] else 'Unknown'
        os.makedirs(f"{ExpectedType}/{ExpectedSubtype}/{sPrefix}", exist_ok=True)
        shutil.move(sInputFolder, f"{ExpectedType}/{ExpectedSubtype}/{sPrefix}/{sInputFolder}")
    else:
        ExpectedType = hIDlookup[sInputFolder]['ExpectedType'] if sInputFolder in hIDlookup and 'ExpectedType' in hIDlookup[sInputFolder] else 'Unknown'
        os.makedirs(f"{ExpectedType}/Other", exist_ok=True)
        shutil.move(sInputFolder, f"{ExpectedType}/Other/{sInputFolder}")
