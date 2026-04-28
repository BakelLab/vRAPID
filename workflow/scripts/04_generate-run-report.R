#!/usr/bin/env Rscript

get_time <- function() format(Sys.time(), "%H:%M:,%S")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
## Libraries
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

suppressPackageStartupMessages(library(connectPDB))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(tidyverse))

#############
# FUNCTIONS #
#############


# PDB database config
pdb_db_group       <- snakemake@params[["pdb_db_group"]]
pdb_db_config      <- snakemake@params[["pdb_db_config"]]


# Utility function for timestamped messages
get_time <- function() format(Sys.time(), "%H:%M:%S")

# Snakemake logging function
log_smk <- function() {
    if (exists("snakemake") && length(snakemake@log) != 0) {
        log_path <- snakemake@log[[1]]
        log_con <- file(log_path, open = "wt")
        sink(log_con, split = FALSE)
        sink(log_con, type = "message")
        return(log_con)
    }
    return(NULL)
}


########
# MAIN #
########

# NOTE:  We define the main block to run as a function.
#        This is done to ensure that we can run an on.exit function if there is any error during
#        the execution of the main block. This ensures that we always clean up any stray tunnels and
#        database connections.

main <- function() {

# -------------------------------------------------------------------------------------------------
# ENSURE CONNECTION CLEANUP
# -------------------------------------------------------------------------------------------------

# Initialize resources so on.exit can clean them
db_pathogendb       <- NULL

# Define on.exit function for the main block
on.exit({
    # Print first while sinks are still active, so these messages go to the log
    cat("\n", get_time(), "[!] Starting connection cleanup...\n")

    # Close DBs (best effort)
    cat("\n", get_time(), "[!] Closing database handles...\n")
    if (!is.null(db_pathogendb)) try(dbDisconnect(db_pathogendb), silent = TRUE)

    # Unsink strictly BEFORE closing the connection
    # Handle the 'message' type sink first
    if (sink.number(type = "message") > 0) {
        sink(type = "message")
    }

    # Handle standard output sinks
    while (sink.number() > 0) {
        sink()
    }

    # Close the file connection
    if (!is.null(log_con)) {
        cat(get_time(), "[!] Log connection closed.\n") # This goes to console now
        try(close(log_con), silent = TRUE)
    }

}, add = TRUE)

# -------------------------------------------------------------------------------------------------
# START LOGGING TO SNAKEMAKE LOG FILE
# -------------------------------------------------------------------------------------------------

log_con = log_smk()

# -------------------------------------------------------------------------------------------------
# DATABASE CONNECTIONS
# -------------------------------------------------------------------------------------------------

# connectPDB functions open MariaDB and MS SQL connections through the SSH tunnels created above.
cat("[-] Connecting to PathogenDB.\n")
db_pathogendb  <- dbconnect_mariadb( db_group = pdb_db_group,    db_config = pdb_db_config)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------
## Input Files & Error Messages Configuration
# -------------------------------------------------------------------------------------------------------------------------------------------------------------

# format output name
output_name <- paste0(snakemake@config[["run_id"]], "_run_report.csv")

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
}
main()
