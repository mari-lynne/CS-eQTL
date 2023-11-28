# Deconvolution -------------------------------------------------------------

# AIMS:
# Using trecase output from step 4; deconvolute cell fractions using cibersort and LM22 signature matrix

## Directories and file names --------------------------------------------------

library(data.table)
library(dplyr)
library(stringr)
library(viridis)


# work_dir = "fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/ciber/"
work_dir <-  "~/Documents/CSeQTL/data/ciber_ase"
plot_dir <- "~/Documents/CSeQTL/data/plots"
setwd(work_dir)

# Input files
sig_tpm <- fread("LM22.txt")

# Output files
sig_fn = file.path(work_dir, "signature.txt") # TPM signature expression matrix
mix_fn = file.path(work_dir, "mixture.txt") # TPM bulk expression matrix

# 1) Gather bulk expression data -------------------------------------------------

# Get a list of all files with the "-TReC.txt" extension
file_list <- list.files(pattern = "-TReC.txt")

# Create an empty data table to store the combined data
bulk <- data.table(ensembl_gene_id = character(0))

# Loop through the files, read them, and combine them into data frame
for (file in file_list) {
  sample_id <- gsub("-TReC.txt", "", file)  # Extract sample_id from the file name
  file_data <- fread(file, header = FALSE, sep = "\t", col.names = c("ensembl_gene_id", sample_id))
  bulk <- merge(bulk, file_data, by = "ensembl_gene_id", all = TRUE)
}

# Replace missing values with zeros
bulk[is.na(bulk)] <- 0

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


write.csv(filtered_genes, file = "~/Documents/CSeQTL/data/ciber_ase/filtered_genes_test_oct.csv")

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

# 3 Run/Load cibersort ---------------------------------------------------------
source(file.path(in_dir, "cibersort_source.R")) # Obtained from CIBERSORT website

# Output 
results <- CIBERSORT(sig_matrix = sig_fn, mixture_file = mix_fn, filename = "DECON_QC",
                     perm = 0, QN = FALSE, absolute = FALSE, abs_method = 'sig.score')

results <- as.data.frame(results)
cat("Cibersort imputation - DONE!")

# Or upload from local docker folder
results <- fread(file = paste0(work_dir, "/docker/output/CIBERSORTx_Adjusted.txt"))

# Extract proportion of transcripts per cell type per sample
ciber_fractions <- as.data.table(results[, 2:23]) # 22 immune subsets
row.names(ciber_fractions) <- results$Mixture

# Update cell sizes ------------------------------------------------------------
# Cell size table (estimated in EPIC paper) 
# TODO update with our own FSC/SSC data
T_cells <- colnames(ciber_fractions)[str_starts(colnames(ciber_fractions), "T")==TRUE]
B_cells <- colnames(ciber_fractions)[str_starts(colnames(ciber_fractions), "B|Plasma")==TRUE]
NK <- colnames(ciber_fractions)[str_starts(colnames(ciber_fractions), "NK")==TRUE]
mono <- colnames(ciber_fractions)[str_starts(colnames(ciber_fractions), "Mono|Macro|Den|Mast")==TRUE]
neut <- colnames(ciber_fractions)[str_starts(colnames(ciber_fractions), "Neut|Eo")==TRUE]

cell_sizes <- c(rep(0.4, length(T_cells)),
               rep(0.4, length(B_cells)),
               rep(0.42, length(NK)),
               rep(1.4, length(mono)),
               rep(0.15, length(neut)))

# Calculate proportion of cell types per sample
pp_hat_ciber = t(apply(ciber_fractions,1,function(xx){
  yy = xx / cell_sizes; yy / sum(yy)
}))


## Combine cell types -----------------------------------------------------------

# Create a new matrix with combined values and columns
cell_combo <- matrix(0, nrow = nrow(pp_hat_ciber), ncol = 5)  # Adjust the number of columns as needed
colnames(cell_combo) <- c("T_cells", "B_cells", "NK", "Mono", "Neut")
row.names(cell_combo) <- row.names(ciber_fractions)

input <- pp_hat_ciber # pp_hat_ciber

cell_combo[, "T_cells"] <- rowSums(input[, T_cells])
cell_combo[, "B_cells"] <- rowSums(input[, B_cells])
cell_combo[, "NK"] <- rowSums(input[, NK])
cell_combo[, "Mono"] <- rowSums(input[, mono])
cell_combo[, "Neut"] <- rowSums(input[, neut])

# write.csv(cell_combo, file = "~/Documents/CSeQTL/data/ciber_ase/pp_hat_ciber_test.csv")

## Plotting ------------------------------------------------------------------

ciber <- as.data.frame(pp_hat_ciber)
ciber <- as.data.frame(cell_combo)
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


# Compare overlapping samples --------------------------------------------------

# 20 samples overlap between sct and lls cohort
# get meta files for each, use intersect topmed_nwdid, then get the rnaseq ids from sct, use that to subset first data set and plot histos



