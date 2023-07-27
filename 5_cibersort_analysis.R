# Deconvolution -------------------------------------------------------------

# AIMS:
# Using trecase output from step 4; deconvolute cell fractions using cibersort and LM22 signature matrix

## Directories and file names --------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(viridis)

plot_dir <- "~/Documents/CSeQTL/data/plots"
# work_dir = "fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/ciber/"
work_dir <-  "~/Documents/CSeQTL/data/ciber_ase"
setwd(work_dir)

# Input files
sig_tpm <- fread("LM22.txt")

# Output files
sig_fn = file.path(work_dir, "signature.txt") # TPM signature expression matrix
mix_fn = file.path(work_dir, "mixture.txt") # TPM bulk expression matrix

# 1) Import bulk expression data -------------------------------------------------

# run script to make signature matrix in same wd
system("bash ~/Documents/CSeQTL/scripts/CSeQTL/5_mixture_prep.sh")

bulk <- fread("total_counts.txt")
# Tidy sample names
colnames(bulk) <- str_extract(colnames(bulk), "\\d*(?=-)")
# rename gene_id col
name <- names(bulk[,1])
bulk <- dplyr::rename(bulk, "ensembl_gene_id" = name)

# Tidy bulk data transcript_ids
bulk$ensembl_gene_id <- str_remove_all(bulk$ensembl_gene_id, "\\.\\d")


# 2) Calculate TPM matrix ---------------------------------------------------------
# function
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# Get gene/transcript length from biomart
# ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# genes <- getBM(
#   attributes=c('ensembl_gene_id', 'hgnc_symbol','start_position','end_position'),
#   mart = ensembl)
# genes$size <- c(genes$end_position - genes$start_position)
# write.csv(genes, file = "~/Documents/CSeQTL/data/ciber_ase/gene_length_info.csv", row.names = FALSE)

genes <- read.csv(file = "~/Documents/CSeQTL/data/ciber_ase/gene_length_info.csv")

# Filter study data for genes with hgcn symbols
filtered_genes <- 
  inner_join(bulk %>% group_by(ensembl_gene_id) %>% mutate(id = row_number()),
             genes %>% group_by(ensembl_gene_id) %>% mutate(id = row_number()), 
             by = c("ensembl_gene_id", "id"))


write.csv(filtered_genes, file = "~/Documents/CSeQTL/data/ciber_ase/filtered_genes_test.csv")

counts <- filtered_genes[, 2:(ncol(bulk))]

# 62,000 genes - 48,000
# Might also have to filter trecase output later/filter duplicates


### Make count matrix ----------------------------------------------------------

counts <- as.matrix(counts)
row.names(counts) <- filtered_genes$hgnc_symbol

bulk_tpm <- as.data.frame(tpm(counts, filtered_genes$size))
head(bulk_tpm)

# test formats
bulk_tpm$`Gene Symbol` <- row.names(bulk_tpm)
bulk_tpm <- bulk_tpm[, c(ncol(bulk_tpm), 1:(ncol(bulk_tpm)-1))] # reorder columns

# write tables
write.table(sig_tpm, file = sig_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(bulk_tpm, file = mix_fn, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# 3 Run cibersort --------------------------------------------------------------
source("5_cibersort_source.R") # Obtained from CIBERSORT website

# Output 
results <- CIBERSORT(sig_matrix = sig_fn, mixture_file = mix_fn, filename = "DECON",
                     perm = 0, QN = FALSE, absolute = FALSE, abs_method = 'sig.score')

# filename = "DECON"
results = as.data.frame(results)
# Extract proportion of transcripts per cell type per sample
pp_bar_ciber <- results[1:22] # 22 immune subsets

# Cell size table (estimated in EPIC paper) 
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

# Adjust expression for cell size S
pp_hat_ciber <- t(apply(pp_bar_ciber, 1, function(xx) {
  yy <- xx / sizes
  yy / sum(yy)
}))

pp_hat_ciber <- pp_hat_ciber[-nrow(pp_hat_ciber),]

write.csv(pp_hat_ciber, file = "~/Documents/CSeQTL/data/ciber_ase/pp_hat_ciber_test.csv")

## Combine cell types -----------------------------------------------------------

# Create a new matrix with combined values and columns
combined_matrix <- matrix(0, nrow = nrow(pp_hat_ciber), ncol = 5)  # Adjust the number of columns as needed
colnames(combined_matrix) <- c("T_cells", "B_cells", "NK", "mono", "neut")
row.names(combined_matrix) <- row.names(pp_hat_ciber)

combined_matrix[, "T_cells"] <- rowSums(pp_hat_ciber[, T_cells])
combined_matrix[, "B_cells"] <- rowSums(pp_hat_ciber[, B_cells])
combined_matrix[, "NK"] <- rowSums(pp_hat_ciber[, NK])
combined_matrix[, "mono"] <- rowSums(pp_hat_ciber[, mono])
combined_matrix[, "neut"] <- rowSums(pp_hat_ciber[, neut])

write.csv(combined_matrix, file = "~/Documents/CSeQTL/data/ciber_ase/pp_hat_ciber_test.csv")

## Plotting ------------------------------------------------------------------

ciber <- as.data.frame(pp_hat_ciber)
# ciber <- as.data.frame(combined_matrix)
ciber$sample_id <- row.names(ciber)

ciber <- ciber %>% group_by(sample_id)
ciber <- melt(ciber, id.vars = c("sample_id"), measure.vars = c(1:(ncol(ciber)-1)))
ciber <- ciber %>% rename(fraction = value,
                        cell_type = variable)


plot_name <- ("ciber_fractions.png")
# Set the output file path
output_file <- file.path(plot_dir, plot_name)

# Start capturing the plot to a PNG file
png(output_file)

ciber %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort Output")

dev.off()


