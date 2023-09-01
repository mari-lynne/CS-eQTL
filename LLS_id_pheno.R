setwd("~/Documents/CSeQTL/data/meta/lls")

library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(tidylog)

# Date 30/08/23

# Aims:
# Match IDs from wgs bcf files (wgs), with rnaseq ids so we can filter out extra samples in the bcf files
# Also get pheontype data

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LLS data input

# Geno bcf ids
# From /fh/scratch/delete90/kooperberg_c/topmed_freeze10/phased/whi_only. bcf file query
nw <- fread(file = "whi_geno_nwids.txt", header = FALSE)

# RNA bam ids
# /fh/scratch/delete90/kooperberg_c/lls_rna/bam_files/bam_files
tor <- fread("whi_rna_torids.txt", header = FALSE)

# Link file 
# /fh/scratch/delete90/kooperberg_c/PAGE_phenotype/WHI_updated_phenotypes
link <- clean_names(read.csv(file = "WHI_all_subjects_PAGE_final_map_genotypes_phenotypes_2023-06-29.csv"))

# Joining

# add joining cols suspected to match names in link, then join
nw$topmed_nwdid <- nw$V1
geno_join <- inner_join(nw, link, by = "topmed_nwdid") # most all match

tor$lls_torid <- tor$V1
rna_join <- inner_join(tor, link, by = "lls_torid") # matched 1,362 samples

# Final join
all_match <- as.data.frame(inner_join(rna_join, geno_join, by = "subject_id"))
# topmed_nwdid is the genotype id, and lls_torid is the bam ID

# Clean Cols
all_match <- all_match[ , !str_detect(colnames(all_match), ".y")]
colnames(all_match) <- str_remove(colnames(all_match), ".x")

# Make csv array:
array <- all_match %>% select(topmed_nwdid)
array$index <- seq(1:nrow(array))
array <- array %>% select(index, topmed_nwdid)
write.table(all_match$topmed_nwdid, file = "lls_geno_filt.txt", row.names = F, col.names = F, quote = F)
write.csv(array, file = "lls_sample_array.csv", quote = FALSE, row.names = FALSE)

# Phenotype ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pheno <- clean_names(read.csv(file = "LLS_updated_phenotypes_2023-05-19.csv"))
pheno <- left_join(all_match, pheno, by = "subject_id")

write.csv(pheno, file = "LLS_ID_pheno_08-23.csv", quote = F, row.names = F)
