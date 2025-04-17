#!/usr/bin/env Rscript

#############
# LIBRARIES #
#############

library(ggplot2)
library(tools)
library(dplyr)

########
# MAIN #
########

# Open plot file
pdf(file = snakemake@output[["pdf"]])

# Input file path
krakeninput <- snakemake@input[["report"]]

# Read kraken report
kraken_table <- try(read.delim(krakeninput, header = FALSE, sep = "\t", stringsAsFactors = FALSE))

# Check if the input is empty
if (inherits(kraken_table, "try-error")) {
    print("Kraken report input is empty!")
}

# Extract names and percentages
var_nam <- unlist(kraken_table[6])
var_per <- unlist(kraken_table[1])

# Helper function to extract and clean percentages
get_percentage <- function(name, var_nam, var_per) {
    pos <- grep(name, var_nam)
    perc <- trimws(var_per[pos[1]])
    perc[is.na(perc)] <- 0
    return(as.numeric(perc))
}

# Extract percentages for each classification
fluA_per <- get_percentage('Influenza A virus', var_nam, var_per)
fluB_per <- get_percentage('Influenza B virus', var_nam, var_per)
host_per <- get_percentage('Eukaryota', var_nam, var_per)
bact_per <- get_percentage('Bacteria', var_nam, var_per)
cov_per <- get_percentage('Severe acute respiratory syndrome coronavirus 2', var_nam, var_per)
E229_per <- get_percentage('Human coronavirus 229E', var_nam, var_per)
NL63_per <- get_percentage('Human coronavirus NL63', var_nam, var_per)
HKU1_per <- get_percentage('Human coronavirus HKU1', var_nam, var_per)
OC43_per <- get_percentage('Human coronavirus OC43', var_nam, var_per)
MPX_per <- get_percentage('Monkeypox virus', var_nam, var_per)

# Define colors and classification labels
mycols <- c("#868686FF", "#EFC000FF", "#0073C2FF", "#CD534CFF", "#2D964DFF", 
            "#C18181", "#825656", "#630000", "#FC8686", "#FF8B00")
classification <- c("IAV", "IBV", "Host", "Bacteria", "SARS-CoV-2", "229E", 
                    "NL63", "HKU1", "OC43", "MPX")

# Combine classification and percentages
n <- c(fluA_per, fluB_per, host_per, bact_per, cov_per, E229_per, 
       NL63_per, HKU1_per, OC43_per, MPX_per)
classification <- paste(classification, n, sep = " ")
classification <- paste(classification, "%", sep = "")
count.data <- data.frame(classification, n)

# Plot the data
ggplot(count.data, aes(x = 2, y = n, fill = classification)) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = mycols) +
    ggtitle("Taxonomic Classification for Reads (Viruses of Interest)") +
    labs(caption = "Viruses of interest are subject to change") +
    theme_void() +
    xlim(0.5, 2.5)

# Close plot file
dev.off()

