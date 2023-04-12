#!/usr/bin/env Rscript

# 16.05.2013 12:56:20 EDT
# Harm van Bakel <hvbakel@gmail.com>
# Divya Kriti <divya.kriti@mssm.edu>

#############
# ARGUMENTS #
#############

#  Col 1: the long flag name. A multi-character string.
#  Col 2: short flag alias of Column 1. A single-character string.
#  Col 3: 0=no argument, 1=required argument, 2=optional argument.
#  Col 4: Possible values: logical, integer, double, complex, character.
#  Col 5: A brief description of the purpose of the option.

library(getopt)
args = matrix(c('krakeninput', 'k', 1, "character", "Input file with taxonomic labels assigned to all sequenced reads",
                'output'     , 'o', 2, "character", "Plot output file prefix",
                'rows'       , 'r', 2, "numeric",   "Number of plot rows per page. Default: 2",
                'cols'       , 'c', 2, "numeric",   "Number of plot columns per page. Default: 1",
                'scale'      , 's', 2, "numeric",   "Plot scaling factor. Default: 0.8",
                'help'       , 'h', 0, "logical",   "Brief help message"
               ), ncol=5,byrow=T);
opt  = getopt(args);

# Specify default argument values
if ( is.null(opt$rows)     ) { opt$rows     = 2    }
if ( is.null(opt$cols)     ) { opt$cols     = 1    }
if ( is.null(opt$scale)    ) { opt$scale    = 0.8 }

# Help message
if ( !is.null(opt$help) || is.null(opt$output)) {
   self=sub("^.*/", "", commandArgs()[4], perl=T)
   cat("\n", getopt(args, command=self, usage=T), "\n")
   q(status=1);
}

#############
# LIBRARIES #
#############

library(ggplot2);
library(tools);
library(dplyr);

########
# MAIN #
########

# Open plot file
pdf(file=opt$output);

# Make kraken plots
kraken_table <- try(read.delim(opt$krakeninput, header = FALSE, sep = "\t", stringsAsFactors = F))
if(inherits(kraken_table, "try-error")){
    print("Kraken report input is empty!")}
    
var_nam <- unlist(kraken_table[6]);
var_per <- unlist(kraken_table[1]);

pos_fluA = grep('Influenza A virus', var_nam);
fluA_per <- trimws(var_per[pos_fluA[1]]);
fluA_per[is.na(fluA_per)] <- 0;
fluA_per <- as.numeric(fluA_per);

pos_fluB = grep('Influenza B virus', var_nam);
fluB_per <- trimws(var_per[pos_fluB[1]]);
fluB_per[is.na(fluB_per)] <- 0;
fluB_per <- as.numeric(fluB_per);

pos_host = grep('Eukaryota', var_nam);
host_per <- trimws(var_per[pos_host[1]]);
host_per[is.na(host_per)] <- 0;
host_per <- as.numeric(host_per);

pos_bact = grep('Bacteria', var_nam);
bact_per <- trimws(var_per[pos_bact[1]]);
bact_per[is.na(bact_per)] <- 0;
bact_per <- as.numeric(bact_per);

pos_cov = grep('Severe acute respiratory syndrome coronavirus 2', var_nam);
cov_per <- trimws(var_per[pos_cov[1]]);
cov_per[is.na(cov_per)] <- 0;
cov_per <- as.numeric(cov_per);

pos_229E = grep('Human coronavirus 229E', var_nam);
E229_per <- trimws(var_per[pos_229E[1]]);
E229_per[is.na(E229_per)] <- 0;
E229_per <- as.numeric(E229_per);

pos_NL63 = grep('Human coronavirus NL63', var_nam);
NL63_per <- trimws(var_per[pos_NL63[1]]);
NL63_per[is.na(NL63_per)] <- 0;
NL63_per <- as.numeric(NL63_per);

pos_HKU1 = grep('Human coronavirus HKU1', var_nam);
HKU1_per <- trimws(var_per[pos_HKU1[1]]);
HKU1_per[is.na(HKU1_per)] <- 0;
HKU1_per <- as.numeric(HKU1_per);

pos_OC43 = grep('Human coronavirus OC43', var_nam);
OC43_per <- trimws(var_per[pos_OC43[1]]);
OC43_per[is.na(OC43_per)] <- 0;
OC43_per <- as.numeric(OC43_per);

pos_MPX = grep('Monkeypox virus', var_nam);
MPX_per <- trimws(var_per[pos_MPX[1]]);
MPX_per[is.na(MPX_per)] <- 0;
MPX_per <- as.numeric(MPX_per);


mycols <- c("#868686FF", "#EFC000FF", "#0073C2FF", "#CD534CFF", "#2D964DFF", "#C18181", "#825656", "#630000", "#FC8686", "#FF8B00");
classification <- c("IAV", "IBV", "Host", "Bacteria", "SARS-CoV-2", "229E", "NL63", "HKU1", "OC43", "MPX");
n = c(fluA_per,fluB_per,host_per,bact_per, cov_per, E229_per, NL63_per, HKU1_per, OC43_per, MPX_per);
classification <- paste(classification, n, sep=" ");
classification <- paste(classification, "%", sep="");
count.data <- data.frame(classification, n);

ggplot(count.data, aes(x = 2, y = n, fill = classification)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = mycols) +
    ggtitle("Taxonomic Classification for Reads (Viruses of Interest)") +
    labs(caption = "Viruses of interest are subject to change") +
    theme_void() +
    xlim(0.5, 2.5);

# Close plot file
dev=dev.off()
