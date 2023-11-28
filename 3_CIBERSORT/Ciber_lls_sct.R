# Deconvolution -------------------------------------------------------------

# AIMS:
# Using trecase output from step 4; deconvolute cell fractions using cibersort and LM22 signature matrix

## Directories and file names --------------------------------------------------

library(data.table)
library(dplyr)
library(tidylog)
library(stringr)
library(viridis)
library(ggplot2)
library(janitor)
library(smarter)
library(patchwork)

# work_dir = "fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/ciber/"
# sct dir = "~/Documents/CSeQTL/data/ciber_ase"
plot_dir <- "~/Documents/CSeQTL/data/plots"

# Local ------------------------------------------------------------------------
# Input files
# sig_tpm <- fread("~/Documents/CSeQTL/data/ciber_ase/LM22.txt")
# bulk <- fread(file.path(work_dir, "/sct_tpm_rnaseq.txt"))

# LLS ==========================================================================

work_dir = "~/Documents/CSeQTL/data/ciber_ase"
sig_fn = file.path(work_dir, "LM22.txt") # TPM signature expression matrix
mix_fn = file.path(work_dir, "mixture_qc.txt") # TPM bulk expression matrix
setwd(work_dir)

# SCT ==========================================================================

work_dir = "~/Documents/whi_sct/rna/results/ciber/R"
# Output files
sig_fn = file.path(work_dir, "signature.txt") # TPM signature expression matrix
mix_fn = file.path(work_dir, "mixture.txt") # TPM bulk expression matrix


# Run Deconvolution ============================================================ 

## Run cibersort
source("~/Documents/CSeQTL/scripts/CSeQTL/5_cibersort_source.R")

results <- CIBERSORT(sig_matrix = sig_fn,
                     mixture_file = mix_fn,
                     filename = "DECON",
                     perm = 0,
                     QN = FALSE,
                     absolute = FALSE,
                     abs_method = 'sig.score')

results = as.data.frame(results)
# Extract proportion of transcripts per cell type per sample
pp_bar_ciber <- results[1:22] # 22 immune subsets

# Cell size table (estimated in EPIC paper from total mRNA)

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

# Update without scaling
sizes <- c(rep(0.4, length(B_cells)),
           rep(0.4, length(T_cells)),
           rep(0.4, length(NK)),
           rep(0.4, length(mono)),
           rep(0.4, length(neut)))


pp_bar_ciber <- rbind(pp_bar_ciber, sizes)

# Adjust expression for cell size S
pp_hat_ciber <- t(apply(pp_bar_ciber, 1, function(xx) {
  yy <- xx / sizes
  yy / sum(yy)
}))

pp_hat_ciber <- pp_hat_ciber[-nrow(pp_hat_ciber),]

# Save interim data 
pp_hat_lls_scale <- pp_hat_ciber
pp_hat_lls_noscale <- pp_hat_ciber

pp_hat_sct_scale <- pp_hat_ciber
pp_hat_sct_ns  <- pp_hat_ciber

write.csv(pp_hat_ciber, file = file.path(work_dir, "pp_hat_ciber_sct_ns.csv"))

# Merge cell types -------------------------------------------------------------
pp_hat_ciber <- read.csv(file = file.path(work_dir, "pp_hat_ciber_sct_ns.csv"))

# Merged column names:
T_cells <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "T")==TRUE]
B_cells <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "B|Plasma")==TRUE]
NK <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "NK")==TRUE]
mono <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "Mono|Macro|Den|Mast")==TRUE]
neut <- colnames(pp_hat_ciber)[str_starts(colnames(pp_hat_ciber), "Neut|Eo")==TRUE]

# Create a new matrix with combined values and columns
combined_matrix <- matrix(0, nrow = nrow(pp_hat_ciber), ncol = 5)  # Adjust the number of columns as needed
colnames(combined_matrix) <- c("T_cells", "B_cells", "NK", "mono", "neut")
row.names(combined_matrix) <- row.names(pp_hat_ciber)
row.names(combined_matrix) <- (pp_hat_ciber$X)

combined_matrix[, "T_cells"] <- rowSums(pp_hat_ciber[, T_cells])
combined_matrix[, "B_cells"] <- rowSums(pp_hat_ciber[, B_cells])
combined_matrix[, "NK"] <- rowSums(pp_hat_ciber[, NK])
combined_matrix[, "mono"] <- rowSums(pp_hat_ciber[, mono])
combined_matrix[, "neut"] <- rowSums(pp_hat_ciber[, neut])

write.csv(combined_matrix, file =  file.path(work_dir, "pp_hat_ciber_merged_sct_ns.csv"))


# Plotting ---------------------------------------------------------------------

# ciber <- as.data.frame(pp_hat_ciber)
ciber <- as.data.frame(combined_matrix)
ciber$sample_id <- row.names(ciber)

ciber <- ciber %>% group_by(sample_id)
ciber <- melt(ciber, id.vars = c("sample_id"), measure.vars = c(1:(ncol(ciber)-1)))
ciber <- ciber %>% rename(fraction = value,
                          cell_type = variable)
# Plot
ciber %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort SCT - Unscaled")

# ciber_lls_scaled <- ciber
# ciber_lls_ns <- ciber
# ciber_sct_scaled <- ciber
# ciber_sct_ns <- ciber

# Check LLS v SCT --------------------------------------------------------------

# Get IDs
sct_ids <- row.names(pp_hat_sct_scale) # rnaseq_ids
lls_ids <- row.names(pp_hat_lls_noscale) # TOR_ids
# Read in ID link file
link <- read.csv(file = "~/Documents/whi_sct/rna/meta/sct_all_covars_Jul23.csv")
link_lls <- read.csv(file = "~/Documents/CSeQTL/data/meta/lls/LLS_ID_pheno_08-23.csv")


# Check Overlap
share <- intersect(link$subject_id, link_lls$subject_id)
# Convert to same ID
share_sct <- link[link$subject_id %in% share, c("subject_id","rnaseq_ids")]
share_lls <- link_lls[link_lls$subject_id %in% share,c("subject_id","lls_torid")]

# Filter tables and Add subject ids for later merging
sub_sct <- ciber_sct_scaled %>% filter(sample_id %in% share_sct$rnaseq_ids) %>%
  mutate(rnaseq_ids = sample_id)
sub_sct <- left_join(sub_sct, share_sct, by = "rnaseq_ids")

sub_lls <- ciber_lls_scaled %>% filter(sample_id %in% share_lls$lls_torid) %>%
  mutate(lls_torid = sample_id)
sub_lls <- left_join(sub_lls, share_lls, by = "lls_torid")

### Save/Load --------------------------------------------------------------------

# save.image(file = "~/Documents/CSeQTL/data/ciber_ase/LLS_SCT_ciber.RData")
load(file = "~/Documents/CSeQTL/data/ciber_ase/LLS_SCT_ciber.RData")

# Subset plotting ==============================================================
a <- sub_sct %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort SCT - Scaled")

b <- sub_lls %>% ggplot(aes(x=sample_id, y=fraction, fill=cell_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete=TRUE, option = "H", direction = -1, alpha = 0.95) +
  theme_light() +
  labs(y = "Fraction \n", x= "Sample ID", title = "Cibersort LLS - Scaled")

(a|b) + plot_layout(guides = "collect")

# Correlations ====================================================

# Correlate fraction values
# Merge tables by subject id

cor_table <- left_join(sub_sct, sub_lls, by = c("subject_id", "cell_type"))
cor_table <- rename(cor_table, sct_fractions = fraction.x,
                    lls_fractions = fraction.y)

# Plot 
cor_table %>% ggplot(aes(x = sct_fractions, y = lls_fractions, color = cell_type)) + geom_point()

ggsave(file = file.path(plot_dir, "sct_lls_cibersort_allcors.png"))

cor_table %>% filter(cell_type == "NK") %>% ggplot(aes(x = sct_fractions, y = lls_fractions)) + geom_point()
# worse cors with monocytes and NK cells

nk <- cor_table %>% filter(cell_type == "NK") %>% ggplot(aes(x = sct_fractions, y = lls_fractions)) + geom_point(color = "#00BF7D") + ggtitle("NK Cells")
mono <- cor_table %>% filter(cell_type == "mono") %>% ggplot(aes(x = sct_fractions, y = lls_fractions)) + geom_point(color = "#00B0F6") + ggtitle("Monocytes")

(nk | mono)
ggsave(file = file.path(plot_dir, "sct_lls_nk_v_mono.png"))

library(scales)
hex <- hue_pal()(length(levels(cor_table$cell_type)))
show_col(hex)
