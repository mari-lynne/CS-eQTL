
library(stringi)
library(stringr)
library(dplyr)
library(janitor)
library(tidylog)

# RNA seq

source("~/scripts/functions.R")
wd <- c("~/Documents/whi_sca/rna") # Where main data is saved
id_dir <- ("~/Documents/whi_sca/rna/ids")

manifest <- read.csv(file = "~/scripts/slurm_22_05/manifest.csv")

regexp <- "[[:digit:]]+"
manifest <- manifest %>%
  mutate(nwgc_sample_id = unlist(str_extract_all(manifest$bam_file, regexp)))
# these match up rnaseq_ids of dataset :) 
write.csv(manifest, file = "~/Documents/whi_sca/rna/ids/manifest2.csv", row.names = FALSE)

# Link nwgc sample id to TOR id
id_link <- clean_names(read.csv(file = "~/Documents/whi_sca/rna/ids/lookup_whi_topmed_to6_rnaseq_1_final.csv"))
id_link$nwgc_sample_id <- as.character(id_link$nwgc_sample_id)
id_link <- inner_join(manifest, id_link, by = "nwgc_sample_id")

id_link <- select(id_link, bam_file, nwgc_sample_id, torid, investigator_id)

# Link TOR id to topmed nwdid
geno_link <- clean_names(read.csv(file = paste(id_dir, "WHI_all_subjects_PAGE_final_map_genotypes_phenotypes_2023-04-22.csv", sep = "/")))
geno_link <- geno_link %>% select(-c("mz_pair", "as311", "as315", "baa23", "ethnic", "aim1", "aim2", "sct", "sample_match"))
geno_link$torid <- geno_link$sct_torid
all_ids <- inner_join(id_link, geno_link, by = "torid") #1029 samples
all_ids <- all_ids %>% select(bam_file, nwgc_sample_id, torid, investigator_id, topmed_nwdid, sct_nwdid)

geno_link$torid <- geno_link$sct_torid
rm(data, rna_ids)
geno_samples <- filter(all_ids, !is.na(topmed_nwdid))

bcf_ids <- fread(file = "~/Documents/CSeQTL/data/meta/geno_ids.txt", header = F)
bcf_ids <- rename(bcf_ids, topmed_nwdid = V1)

manifest <- geno_samples %>% select(bam_file, topmed_nwdid)
geno_samples <- rename(geno_samples, geno_file = topmed_nwdid)
geno_samples <- mutate(geno_samples, geno_file = str_c(geno_samples$geno_file, ".bcf"))


write.csv(manifest, file = "~/Documents/whi_sca/rna/ids/manifest_geno.csv", row.names = FALSE)
write.csv(manifest, file = "~/scripts/slurm_22_05/manifest_geno.csv", row.names = FALSE)

table <- select(geno_samples, geno_file)
write.table(table, file = "~/scripts/slurm_22_05/geno_files_filt.txt", row.names = F, quote = F, sep = "\t")

# check <- inner_join(bcf_ids, geno_samples, by = "topmed_nwdid")
