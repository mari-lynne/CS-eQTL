#!/usr/bin/Rscript

# Date Oct 16 2023

# AIMS:
# Using trecase output, deconvolute cell fractions using cibersort and LM22 signature matrix
# script saved in cseqtl/scripts/split_scripts.R 

## Directories and file names --------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(viridis)
library(ggplot2)
library(smarter)

# Input Directories   
in_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/resources" # cibersort script, gene length info
trec_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE" # Total Read Count

# Output Directories
out_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/ciber/lls"
plot_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/ciber/lls"
setwd(out_dir)

# Input files
sig_tpm <- fread(file.path(in_dir, "LM22.txt"))
str(sig_tpm)

# Output files
sig_fn = file.path(out_dir, "signature.txt") # TPM signature expression matrix
mix_fn = file.path(out_dir, "mixture.txt") # TPM bulk expression matrix

# 1) Gather bulk expression data -------------------------------------------------

# Get a list of all files with the "-TReC.txt" extension
file_list <- list.files(path=trec_dir, pattern = "-TReC.txt")
length(file_list) # 1325 files

# Create an empty data table to store the combined data
bulk <- data.table(ensembl_gene_id = character(0))

# Loop through the files, read them, and combine them into data frame
for (file in file_list) {
  sample_id <- gsub("-TReC.txt", "", file)  # Extract sample_id from the file name
  file_data <- fread(file.path(trec_dir, file), header = FALSE, sep = "\t", col.names = c("ensembl_gene_id", sample_id))
  bulk <- merge(bulk, file_data, by = "ensembl_gene_id", all = TRUE)
}

# Replace missing values with zeros
bulk[is.na(bulk)] <- 0

# Tidy bulk data transcript_ids
bulk$ensembl_gene_id <- str_remove_all(bulk$ensembl_gene_id, "\\.\\d")

cat("Compiled TREC data - DONE!")

## 2) Calculate TPM matrix -----------------------------------------------------

genes <- read.csv(file.path(in_dir, "gene_length_info.csv"))
# Filter study data for genes with hgcn symbols, i.e are in gene length table
filtered_genes <-
  inner_join(
    bulk %>% group_by(ensembl_gene_id) %>% mutate(id = row_number()),
    genes %>% group_by(ensembl_gene_id) %>% mutate(id = row_number()),
    by = c("ensembl_gene_id", "id")
  )

cat("Filtered TREC data - DONE!")
write.csv(filtered_genes, file.path(out_dir, "filtered_genes.csv"))

# Remake total read count table
counts <- filtered_genes[, 2:(ncol(bulk))] 
counts <- as.matrix(counts)
row.names(counts) <- filtered_genes$hgnc_symbol

# Calculate tpm:
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
bulk_tpm <- as.data.frame(tpm(counts, filtered_genes$size))
head(bulk_tpm)

# Test formats
bulk_tpm$`Gene Symbol` <- row.names(bulk_tpm)
bulk_tpm <- bulk_tpm[, c(ncol(bulk_tpm), 1:(ncol(bulk_tpm)-1))] # reorder cols

cat("Calculated TPM - DONE!")

# Write data
write.table(sig_tpm, file = sig_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(bulk_tpm, file = mix_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 3) Run cibersort --------------------------------------------------------------
source(file.path(in_dir, "cibersort_source.R")) # Obtained from CIBERSORT website

# Output 
results <- CIBERSORT(sig_matrix = sig_fn, mixture_file = mix_fn, filename = "DECON2",
                     perm = 0, QN = FALSE, absolute = FALSE, abs_method = 'sig.score')

results <- as.data.frame(results)
cat("Cibersort imputation - DONE!")

# Extract proportion of transcripts per cell type per sample
pp_bar_ciber <- results[1:22] # 22 immune subsets


#4) Combine cell types ---------------------------------------------------------

# Adjust for cell sizes 

# Cell size table (estimated and taken from EPIC paper) 
# TODO update with our own FSC/SSC data
T_cells <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "T")==TRUE]
B_cells <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "B|Plasma")==TRUE]
NK <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "NK")==TRUE]
mono <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "Mono|Macro|Den|Mast")==TRUE]
neut <- colnames(pp_bar_ciber)[str_starts(colnames(pp_bar_ciber), "Neut|Eo")==TRUE]

sizes <- c(rep(0.4, length(B_cells)),
           rep(0.4, length(T_cells)),
           rep(0.42, length(NK)),
           rep(1.4, length(mono)),
           rep(0.15, length(neut)))

pp_bar_ciber <- rbind(pp_bar_ciber, sizes)
str(pp_bar_ciber)

# Adjust expression for cell size S
pp_hat_ciber <- t(apply(pp_bar_ciber, 1, function(xx) {
  yy <- xx / sizes
  yy / sum(yy)
}))

pp_hat_ciber <- pp_hat_ciber[-nrow(pp_hat_ciber),]

write.csv(pp_hat_ciber, file.path(out_dir, "pp_ciber.csv"))

# Create a new matrix with combined values and columns
combined_matrix <- matrix(0, nrow = nrow(pp_hat_ciber), ncol = 5)  # Adjust the number of columns as needed
colnames(combined_matrix) <- c("T_cells", "B_cells", "NK", "mono", "neut")
row.names(combined_matrix) <- row.names(pp_hat_ciber)

combined_matrix[, "T_cells"] <- rowSums(pp_hat_ciber[, T_cells])
combined_matrix[, "B_cells"] <- rowSums(pp_hat_ciber[, B_cells])
combined_matrix[, "NK"] <- rowSums(pp_hat_ciber[, NK])
combined_matrix[, "mono"] <- rowSums(pp_hat_ciber[, mono])
combined_matrix[, "neut"] <- rowSums(pp_hat_ciber[, neut])

cat("Collapsed Cibersort data - DONE!")

write.csv(combined_matrix, file.path(out_dir, "pp_ciber_combo.csv"))

# 5) Plot Cibersort results ----------------------------------------------------

# 1) LM22 all cell types
ciber <- as.data.frame(pp_hat_ciber)
ciber$sample_id <- row.names(ciber)
# Reformat data
ciber <- ciber %>% group_by(sample_id)
ciber <- melt(ciber, id.vars = c("sample_id"), measure.vars = c(1:(ncol(ciber)-1)))
ciber <- ciber %>% rename(fraction = value,
                          cell_type = variable)

plot_name <- ("ciber_fractions_lm22.png")
output_file <- file.path(plot_dir, plot_name)

# Plot
p1 <- ciber %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort Imputed Cell Fractions")
ggsave(output_file, plot = p1, device = "png")

# 2) Collapsed cell types

ciber <- as.data.frame(combined_matrix)
ciber$sample_id <- row.names(ciber)
ciber <- ciber %>% group_by(sample_id)
ciber <- melt(ciber, id.vars = c("sample_id"), measure.vars = c(1:(ncol(ciber)-1)))
ciber <- ciber %>% rename(fraction = value,
                          cell_type = variable)

plot_name <- ("ciber_fractions_epic.png")
output_file <- file.path(plot_dir, plot_name)

# Plot collapsed
p2 <- ciber %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort Imputed Cell Fractions")
ggsave(output_file, plot = p2, device = "png")

cat("Plotting - Done!")