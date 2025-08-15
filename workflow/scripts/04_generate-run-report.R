#!/usr/bin/env Rscript

get_time <- function() format(Sys.time(), "%H:%M:,%S")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
## Libraries
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# Check for 'connectPDB' and install if necessary
if (!requireNamespace("connectPDB")) {
  if (!requireNamespace("devtools")) {
    install.packages("devtools", repos = 'http://cran.us.r-project.org')
  }
  library(devtools)
  pdbtoken <- readLines("~/.pdbtoken")
  install_github("BakelLab/pathogendb-connect", auth_token = pdbtoken)
}
# List of libraries to load
packages <- c("getopt", "lubridate", "tidyverse","connectPDB")

# Function to check if package is installed and load it
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = 'http://cran.us.r-project.org')
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Load all required packages
lapply(packages, load_or_install)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
## Input Files & Error Messages Configuration
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# format output name
output_name <- paste0(snakemake@config[["run_id"]], "_run_report.csv")

#############
# FUNCTIONS #
#############

db_pathogendb <- dbconnect_mariadb('vanbah01_pathogens')

# -------------------------------------------------------------------------------------------------------------------------------------------------------------
## Pathogen DB Info
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

### ASSEMBLY DATA ###

cat(paste0("\n", get_time(), " [-] Obtain Assembly data\n"))

# get assembly data, filtering for assembly run ID
assembly_data <-  suppressWarnings(dbGetQuery(db_pathogendb, "SELECT * FROM tCEIRS_assemblies"))
isolates_data <- suppressWarnings(dbGetQuery(db_pathogendb, "SELECT Isolate_ID,Sample_Name,Flu_Type,Collection_Date FROM tCEIRS_Isolates"))
extract_data <- suppressWarnings(dbGetQuery(db_pathogendb, "SELECT Extract_ID,Isolate_ID,Sample_Systematic_ID,PSP_investigation_ID, Contract_Year FROM tCEIRS_Extracts"))

# standardize column names
# in assembly data
colnames(assembly_data)[c(6,9,11, 17)] <- c("Subtype Found", "Uniquely Mapped Read Percent", "Genome Completeness", "Consensus variants (15%)")
colnames(assembly_data) <- str_to_upper(gsub(" ", "_", colnames(assembly_data)))
# in isolates data
colnames(isolates_data)[c(3, 4)] <- c("Virus Type", "Sampling Date")
colnames(isolates_data) <- str_to_upper(gsub(" ", "_", colnames(isolates_data)))
# in extract data
colnames(extract_data)[3] <- "Sample ID"
colnames(extract_data) <- str_to_upper(gsub(" ", "_", colnames(extract_data)))

# join to get final assembly data like in PDB
inter_1 <- suppressMessages(full_join(isolates_data, extract_data))
inter_2 <- suppressMessages(full_join(assembly_data, inter_1))

cat(paste0("\n", get_time(), " [-] Reformat Assembly data\n"))

inter_2 <- inter_2 %>% filter(ASSEMBLY_RUN == snakemake@config[["run_id"]])

# drop isolates ID and date created columns after merging on isolates ID
inter_2 <- inter_2[-c(5, 7, 21)]

dbDisconnect(db_pathogendb)

inter_2$COLLABORATOR <- unlist(lapply(strsplit(inter_2$SAMPLE_ID, "_", fixed = TRUE),function(x) x[1]))

inter_3 <- inter_2 %>% count(COLLABORATOR, ASSEMBLY_QUALITY, ASSEMBLY_STATUS, VIRUS_TYPE, SUBTYPE_FOUND)
colnames(inter_3) <- str_to_title(colnames(inter_3))
colnames(inter_3)[6] <- "Count"

cat(paste0("\n", get_time(), " [-] Write to file\n"))

write_csv(inter_3, output_name)
write(paste0("\nTotal Samples in PDB = ", sum(inter_3$Count), "\n"), file = output_name, append = T)

cat(paste0("\n", get_time(), " [-] Done\n"))
