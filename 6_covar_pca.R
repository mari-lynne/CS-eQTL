# Set-up -----------------------------------------------------------------------
# Libraries
library(stringi)
library(stringr)
library(dplyr)
library(janitor)
library(tidylog)
source("~/scripts/functions.R")

# Directories and files --------------------------------------------------------
wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
setwd(wd)

# Set up ID file/manifest ------------------------------------------------------
manifest <- read.csv(file = "~/scripts/slurm_22_05/manifest.csv")

regexp <- "[[:digit:]]+"
manifest <- manifest %>%
  mutate(nwgc_sample_id = unlist(str_extract_all(manifest$bam_file, regexp)))
# these match up rnaseq_ids of dataset :) 
# write.csv(manifest, file = "~/Documents/whi_sca/rna/ids/manifest2.csv", row.names = FALSE)

# Link nwgc sample id to TOR id
id_link <- clean_names(read.csv(file = "~/Documents/whi_sca/rna/ids/lookup_whi_topmed_to6_rnaseq_1_final.csv"))
id_link$nwgc_sample_id <- as.character(id_link$nwgc_sample_id)
id_link <- inner_join(manifest, id_link, by = "nwgc_sample_id")

# Rename vars
id_link <- id_link %>% rename(flagged_seq_qc = flagged_seq_qc_metric,
                                  pick_plate1 = pick_plate_1_indicates_which_samples_were_prepped_together)
# Make new plate var
id_link <- id_link %>%
  mutate(plate = coalesce(id_link$pick_plate1, id_link$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))
# Update factors
id_link$plate <- as.factor(id_link$plate)
# Filter flagged QCC
id_link <- filter(id_link, flagged_seq_qc != "Yes")
id_link <- select(id_link, bam_file, geno_file, nwgc_sample_id, torid, investigator_id, topmed_nwdid, sct_nwdid, plate, assay_sex)

# Link TOR id to topmed nwdid
geno_link <- clean_names(read.csv(file = paste(id_dir, "WHI_all_subjects_PAGE_final_map_genotypes_phenotypes_2023-04-22.csv", sep = "/")))
geno_link <- geno_link %>% select(-c("mz_pair", "as311", "as315", "baa23", "ethnic", "aim1", "aim2", "sct", "sample_match"))
geno_link$torid <- geno_link$sct_torid
all_ids <- inner_join(id_link, geno_link, by = "torid") #1029 samples
all_ids <- mutate(all_ids, geno_file = str_c(all_ids$topmed_nwdid, ".bcf"))
all_ids <- all_ids %>%
  select(bam_file, nwgc_sample_id, geno_file, topmed_nwdid, torid, investigator_id, sct_nwdid, subject_id, plate, assay_sex)

# Extra phenotype data ---------------------------------------------------------

# Phenotype file (made in sct_id_pheno_links.R)
pheno <- fread(paste0(meta_dir, "/rnaseq_pheno_all.txt"))
pheno <- select(pheno, rnaseq_ids, subject_id, bmi_t0, age)

geno_samples <- filter(all_ids, !is.na(topmed_nwdid))
geno_samples <- left_join(geno_samples, pheno, by = "subject_id")

# RNAseq PCA -------------------------------------------------------------------

# See ~/Documents/whi_sca/rna/results/pca.R
pc_top <- data.frame(fread(file = "~/Documents/whi_sca/rna/results/pc_top.txt"))
pc_top <- pc_top %>% dplyr::rename(rnaseq_ids = V1)
covar <- left_join(geno_samples, pc_top, by = "rnaseq_ids")

write.csv(covar, file = "~/Documents/whi_sca/rna/meta/sct_cseqtl_test_covars.csv", row.names = F)


# Library size -----------------------------------------------------------------

# One of the non-intercept columns should correspond to centered log-transformed library size.
# The rownames(XX) needs to be specified.

load(file = "~/Documents/whi_sca/rna/results/sct_dge.RData")
data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 

# RNA seq
data <- read_omic(name = data_name, wd = wd)
pheno <- covar[covar$nwgc_sample_id %in% c("941357", "941437") , 12:27]

# Filter samples
sub <- data[,colnames(data)%in% covar$rnaseq_ids]

# Gene_info
gene_info <- data[, 1:2]

# Calculate library sizes
counts <- as.matrix(sub)
dge <- DGEList(counts=counts, samples=covar, genes=gene_info)

# poss redo manually as colsums(mat) > centre + log transform lsizes, save as vector
dge <- calcNormFactors(dge)
dge$samples$lib.size

write.csv(dge$samples, file = "~/Documents/whi_sca/rna/meta/sct_cseqtl_test_covars.csv", row.names = T)

# GWAS PCs ---------------------------------------------------------------------
pc_geno <- read.csv(file = "sct_wgs_PCs1-12.csv")
pc_geno <- pc_geno %>% rename(topmed_nwdid = sample.id)
test <- left_join(covar, pc_geno, by = "topmed_nwdid")
