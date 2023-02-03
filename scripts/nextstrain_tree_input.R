#create metadatafile and update sequence file for nextstrain

library(odbc)
library(DBI)
library(RMySQL)
library("readxl")
library(lubridate)
library(seqinr)
library('plyr')
library(dplyr)
library(getopt)

args = matrix(c('fasta_file'     , 'f', 2, "character", "Database host",
                'nextstrain_meta', 'm', 2, "character", "Database name",
                'help'           , 'h', 0, "logical",   "Brief help message"
), ncol=5,byrow=T);
opt  = getopt(args);

next_dir<-Sys.getenv("nextstrain_dir")
next_meta<-Sys.getenv("nextstrain_meta")
pvi_fasta<-Sys.getenv("PVI_fasta")
admit_status_report<-Sys.getenv("admit_status_report")
driver_file<-Sys.getenv("MSSQL_driver")


#Read PHI DB credentials
phi_cred<-as.data.frame(t(read.table('~/.my.cnf.phidb',sep='=')))
phi_cred<-phi_cred %>%
  mutate_all(as.character)

con <- dbConnect(odbc(),
                 Driver = driver_file,
                 Server = phi_cred$host,
                 Database = phi_cred$database,
                 UID = phi_cred$user,
                 PWD = phi_cred$password,
                 Port = phi_cred$port)


#Read PDB credentials
pdb_cred<-as.data.frame(t(read.table('~/.my.cnf.pdbrw',sep='=')))
pdb_cred<-pdb_cred %>%
  mutate_all(as.character)

# Connect to PDB
mydb = dbConnect(MySQL(), user=pdb_cred$user, password=pdb_cred$password, dbname=pdb_cred$database, host=pdb_cred$host)

b_rs = dbSendQuery(mydb, "SELECT * FROM tCEIRS_Isolates ti,tCEIRS_Extracts te where ti.Flu_Type LIKE'SARS%' and ti.Sample_Name like 'PV%' and ti.Isolate_ID=te.Isolate_ID")
isolate_table = fetch(b_rs, n=-1);
isolate_table[,1]<-NULL

c_rs = dbSendQuery(mydb, "SELECT * FROM tCEIRS_assemblies")
assemblies_table = fetch(c_rs, n=-1);

isolate_table<-merge(isolate_table,assemblies_table,by='Extract_ID')
isolate_table<-isolate_table[isolate_table$assembly_status=='Complete',] #$isolate_table$coronavirus_percent>=4 & 
assembly_list<-paste0(shQuote(unique(isolate_table$assembly_ID)),collapse= ',')
d_rs<-dbSendQuery(mydb, paste("SELECT * FROM tCEIRS_assembly_sequences where AssemblyID IN (",assembly_list,")",sep=''))
seq_table = fetch(d_rs, n=-1);
isolate_table<-merge(isolate_table,seq_table,by.x = 'assembly_ID',by.y='AssemblyID')

#remove duplicated sample (select samples with low intra host variants)
library(dplyr)
isolate_table<-isolate_table %>%
  group_by(Sample_Name) %>%
  arrange(Variant_pos_sum_15pct) %>%
  slice(1L)
isolate_table$Sample_Name<-trimws(isolate_table$Sample_Name)
isolate_table$Sample_Name<-gsub(' ','_',isolate_table$Sample_Name)

write.fasta(as.list(isolate_table$Sequence), isolate_table$Sample_Name, pvi_fasta, open = "w", nbchar = 60, as.string = FALSE)

b_rs = dbSendQuery(mydb, "select ti.*,th.hospital_abbreviation from tIsolates ti, tHospitals th where ti.hospital_ID=th.hospital_ID")
isolate_table_pdb = fetch(b_rs, n=-1)

b_rs = dbSendQuery(mydb, "select * from tPVI_Surveillance")
pvi_table = fetch(b_rs, n=-1)

isolate_table<-merge(isolate_table,pvi_table[,c('specimen_ID','cml_isolate_id')],by.x='Sample_Name',by.y='specimen_ID',all.x = T)
isolate_table<-merge(isolate_table,isolate_table_pdb[,c('isolate_ID','eRAP_ID','hospital_abbreviation')],by.x='cml_isolate_id',by.y='isolate_ID',all.x=T)
isolate_table<-as.data.frame(isolate_table)
isolate_table$eRAP_ID<-as.character(isolate_table$eRAP_ID)
erap_list<-paste0(shQuote(unique(isolate_table$eRAP_ID)),collapse= ',')

b_rs = dbSendQuery(mydb, paste("select * from tPatientEncounter where eRAP_ID IN (",erap_list,")",sep=''))
encounter_table = fetch(b_rs, n=-1);
encounter_table<-encounter_table[encounter_table$eRAP_ID%in%isolate_table$eRAP_ID,]

library(dplyr)
encounter_table_uniq<-encounter_table %>%
  group_by(eRAP_ID) %>%
  arrange(desc(start_date)) %>%
  slice(1L)

isolate_table<-merge(isolate_table,encounter_table_uniq[,c('eRAP_ID','age','sex','race','ethnic_group','zip_code')],by.x='eRAP_ID',by.y='eRAP_ID',all.x=T)
isolate_table$zip_code<-gsub("-.*","",isolate_table$zip_code)

b_rs = dbSendQuery(mydb, "select * from tZipCodes")
zipcodes_table = fetch(b_rs, n=-1);

isolate_table<-merge(isolate_table,zipcodes_table,by.x='zip_code',by.y='zipcode',all.x=T)

isolates_lookup<-dbGetQuery(con,"Select * from tIsolates")
isolate_table<-merge(isolate_table,isolates_lookup[,c('anon_id','mrn_id')],by.x='cml_isolate_id',by.y='anon_id',all.x = T)

admit_status<-read.csv2(admit_status_report,sep='\t')
admit_status$hosp_status<-'discharged'
admit_status[admit_status$STILL_ADMITTED==1,'hosp_status']<-'still_admitted'
admit_status[admit_status$DECEASED==1,'hosp_status']<-'deceased'

isolate_table<-merge(isolate_table,admit_status[,c('MRN','hosp_status')],by.x='mrn_id',by.y='MRN',all.x = T)

pvi_lookup<-dbGetQuery(con,"Select * from tPVI_Surveillance")
isolate_table<-merge(isolate_table,pvi_lookup[,c('SpecimenId','Gender','Age','CollectionLocation')],by.x='Sample_Name',by.y='SpecimenId',all.x = T)
isolate_table[is.na(isolate_table$age),'age']<-isolate_table[is.na(isolate_table$age),'Age']
isolate_table[is.na(isolate_table$sex),'sex']<-isolate_table[is.na(isolate_table$sex),'Gender']
isolate_table$CollectionLocation<-gsub('Mount Sinai Hospital','MSH',isolate_table$CollectionLocation)
isolate_table$CollectionLocation<-gsub('Mount Sinai Queens','MSQ',isolate_table$CollectionLocation)
isolate_table$CollectionLocation<-gsub('Mount Sinai St. Lukes','STLUKE',isolate_table$CollectionLocation)
isolate_table$CollectionLocation<-gsub('Mount Sinai West','MSW',isolate_table$CollectionLocation)
isolate_table$CollectionLocation<-gsub('Mount Sinai Brooklyn','BRKLYN',isolate_table$CollectionLocation)

isolate_table[is.na(isolate_table$hospital_abbreviation),'hospital_abbreviation']<-isolate_table[is.na(isolate_table$hospital_abbreviation),'CollectionLocation']

nextstrain_metadata<-read.delim2(next_meta,sep='\t',header=T,quote = "")
nextstrain_metadata$sinai_samples<-'?'
nextstrain_metadata$NY_samples<-'?'
nextstrain_metadata$hospital<-'?'
nextstrain_metadata$zipcode<-'?'
nextstrain_metadata$city<-'?'
nextstrain_metadata$borough<-'?'
nextstrain_metadata$neighborhood<-'?'
nextstrain_metadata$race<-'?'
nextstrain_metadata$hosp_status<-'?'
nextstrain_metadata$psp_investigation_id<-'?'

library(dplyr)
nextstrain_metadata<-nextstrain_metadata %>%
  mutate_all(as.character)

isolate_table<-isolate_table %>%
  mutate_all(as.character)

for(row in 1:nrow(isolate_table)){
  
  sample=isolate_table[row,'Sample_Name']
  date=isolate_table[row,'Collection_Date']
  zipcode=isolate_table[row,'zip_code']
  city=isolate_table[row,'city']
  borough=isolate_table[row,'borough']
  neighborhood=isolate_table[row,'neighborhood']
  gender=isolate_table[row,'sex']
  race=isolate_table[row,'race']
  age=isolate_table[row,'age']
  hospital=isolate_table[row,'hospital_abbreviation']
  length=isolate_table[row,'Sequence_length']
  hosp_status=isolate_table[row,'hosp_status']
  psp_investigation_id=isolate_table[row,'PSP_investigation_ID']
  
  
  insert_row=c(sample,'ncov','?','?',date,'North America','USA','New York','?','?','?','?','genome',length,'Human',age,gender,'?','?','?',
               'MSHS Clinical Microbiology Laboratories','MSHS Pathogen Surveillance Program','?','?','?','?','2020-03-16','NY','SINAI',hospital,
               zipcode,city,borough,neighborhood,race,hosp_status,psp_investigation_id)
  
  nextstrain_metadata<-rbind(nextstrain_metadata,insert_row)
}

nextstrain_metadata[nextstrain_metadata$division=='New York' & nextstrain_metadata$sinai_samples=='?','NY_samples']<-'NY'
nextstrain_metadata[grepl('NY1|NY2',nextstrain_metadata$strain),'sinai_samples']<-'NY'
nextstrain_metadata[grepl('NY1|NY2',nextstrain_metadata$strain),'NY_samples']<-'SINAI'

nextstrain_metadata[grepl('clone',nextstrain_metadata$strain) & !is.na(nextstrain_metadata$strain),'sinai_samples']='cultured'
nextstrain_metadata[grepl('clone',nextstrain_metadata$strain) & !is.na(nextstrain_metadata$strain),'NY_samples']='cultured'

system(paste('cp ',next_dir,'defaults/include_orig.txt ',
             next_dir,'/defaults/include.txt',
             sep=''))

system(paste('cp ',next_dir,'/defaults/exclude_orig.txt ',
             next_dir,'/defaults/exclude.txt',
             sep=''))

#Include New York genomes
write(nextstrain_metadata[nextstrain_metadata$division=='New York' & !grepl('USA.*.PV0',nextstrain_metadata$strain),'strain'],
      paste(next_dir,'/defaults/include.txt',sep=''),append = T,sep='\n')

#Include pacbio sequenced genomes
write(nextstrain_metadata[grepl('NY1|NY2',nextstrain_metadata$strain),'strain'],
      paste(next_dir,'/defaults/include.txt',sep=''),append = T,sep='\n')

#Include all the genomes filtered from PDB
write(isolate_table$Sample_Name,
      paste(next_dir,'/defaults/include.txt',sep=''),append = T,sep='\n')

write(nextstrain_metadata[grepl('USA.*.PV0',nextstrain_metadata$strain) & !grepl('NY1|NY2',nextstrain_metadata$strain),'strain'],
      paste(next_dir,'/defaults/exclude.txt',sep=''),append = T,sep='\n')

#Load Chilean genomes

b_rs = dbSendQuery(mydb, "SELECT * FROM tCEIRS_Isolates ti,tCEIRS_Extracts te where ti.Collaborator_ID=7 and ti.Flu_Type='SARS-CoV-2' and ti.Isolate_ID=te.Isolate_ID")
isolate_table = fetch(b_rs, n=-1);
isolate_table[,1]<-NULL

c_rs = dbSendQuery(mydb, "SELECT * FROM tCEIRS_assemblies")
assemblies_table = fetch(c_rs, n=-1);

isolate_table<-merge(isolate_table,assemblies_table,by='Extract_ID')
isolate_table<-isolate_table[isolate_table$coronavirus_percent>=4 & isolate_table$assembly_status=='Complete',]

assembly_list<-paste0(shQuote(unique(isolate_table$assembly_ID)),collapse= ',')
d_rs<-dbSendQuery(mydb, paste("SELECT * FROM tCEIRS_assembly_sequences where AssemblyID IN (",assembly_list,")",sep=''))
seq_table = fetch(d_rs, n=-1);
isolate_table<-merge(isolate_table,seq_table,by.x = 'assembly_ID',by.y='AssemblyID')

write.fasta(as.list(isolate_table$Sequence), isolate_table$Sample_Name, pvi_fasta, open = "a", nbchar = 60, as.string = FALSE)

for(row in 1:nrow(isolate_table)){
  
  sample=isolate_table[row,'Sample_Name']
  date=isolate_table[row,'Collection_Date']
  length=isolate_table[row,'Sequence_length']
  host=isolate_table[row,'Host_Common_Name']
  
  insert_row=c(sample,'ncov','?','?',date,'South America','Chile','Santiago','?','?','?','?','genome',length,host,'?','?','?','?','?',
               'Laboratory of Molecular Virology, Pontificia Universidad Catolica de Chile','MSHS Pathogen Surveillance Program',
               '?','?','?','?','2020-03-16','CHILE','?','?','?','?','?','?','?','?','?')
  
  nextstrain_metadata<-rbind(nextstrain_metadata,insert_row)
}

write(isolate_table$Sample_Name,
      paste(next_dir,'/defaults/include.txt',sep=''),append = T,sep='\n')


b_rs = dbSendQuery(mydb, "SELECT * FROM tCEIRS_Isolates ti,tCEIRS_Extracts te where ti.Collaborator_ID=20 and ti.Flu_Type='SARS-CoV-2' and ti.Isolate_ID=te.Isolate_ID")
isolate_table = fetch(b_rs, n=-1);
isolate_table[,1]<-NULL

c_rs = dbSendQuery(mydb, "SELECT * FROM tCEIRS_assemblies")
assemblies_table = fetch(c_rs, n=-1);

isolate_table<-merge(isolate_table,assemblies_table,by='Extract_ID')
isolate_table<-isolate_table[(isolate_table$Assembly_quality=='Passed' & isolate_table$assembly_status=='Complete') | isolate_table$assembly_ID=='6906',]

assembly_list<-paste0(shQuote(unique(isolate_table$assembly_ID)),collapse= ',')
d_rs<-dbSendQuery(mydb, paste("SELECT * FROM tCEIRS_assembly_sequences where AssemblyID IN (",assembly_list,")",sep=''))
seq_table = fetch(d_rs, n=-1);
isolate_table<-merge(isolate_table,seq_table,by.x = 'assembly_ID',by.y='AssemblyID')

write.fasta(as.list(isolate_table$Sequence), isolate_table$Sample_Name, pvi_fasta, open = "a", nbchar = 60, as.string = FALSE)

for(row in 1:nrow(isolate_table)){
  
  sample=isolate_table[row,'Sample_Name']
  date=isolate_table[row,'Collection_Date']
  length=isolate_table[row,'Sequence_length']
  host=isolate_table[row,'Host_Common_Name']
  
  insert_row=c(sample,'ncov','?','?',date,'North America','USA','New York','?','?','?','?','genome',length,host,'?','?','?','?','?',
               '?','?','?','?','?','?','2020-03-16','STOOL','STOOL','?','?','?','?','?','?','?','?')
  
  nextstrain_metadata<-rbind(nextstrain_metadata,insert_row)
}

write(isolate_table$Sample_Name,
      paste(next_dir,'/defaults/include.txt',sep=''),append = T,sep='\n')

nextstrain_metadata[is.na(nextstrain_metadata)] <- '?'
nextstrain_metadata[nextstrain_metadata==''] <- '?'

write.table(nextstrain_metadata,file=paste(next_dir,'/data/metadata.tsv',sep=''),
            sep='\t',row.names = F,quote = FALSE)