#!/usr/bin/Rscript

# Aims
# Reformat snp gene matrix for CSeQTL input

library(reshape2)

# Inherit bash args
args <- commandArgs(trailingOnly = TRUE)
# Extract the variables/values
gene_id <- args[1]
out_dir <- args[2]

# Load the data without headers
data <- read.table(file = paste0(out_dir, "/", gene_id, "_temp.txt"),
                   sep = "\t", header = FALSE, col.names = c("SAMPLE", "SNP_ID", "GENOTYPE"))

# Pivot the data using reshape2 
format_data <- dcast(data, SNP_ID ~ SAMPLE, value.var = "GENOTYPE")

# Save the formatted data
write.table(format_data, file = paste0(out_dir, "/", gene_id, ".txt"),
            col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
