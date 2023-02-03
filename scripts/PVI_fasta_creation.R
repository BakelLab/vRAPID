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
#phi_cred<-as.data.frame(t(read.table('~/.my.cnf.phidb',sep='=')))
#phi_cred<-phi_cred %>%
 # mutate_all(as.character)

#con <- dbConnect(odbc(),
 #                Driver = driver_file,
  #               Server = phi_cred$host,
   #              Database = phi_cred$database,
    #             UID = phi_cred$user,
     #            PWD = phi_cred$password,
      #           Port = phi_cred$port)


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

