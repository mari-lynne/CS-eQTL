library("smarter")
library("CSeQTL")
library("stringr")
library("dplyr")
library("tidylog")
library("caret")

# Data Input Format: -----------------------------------------------------------

## Per gene, and per SNP all samples data: 

# SNP	# phased genotype vector

# TREC 	# TReC vector (Per gene!)
# hap2	# 2nd haplotype counts
# ASREC 	# total haplotype counts = hap1 + hap2
# PHASE 	# Indicator vector of whether or not to use haplotype counts

# SNP phased genotype vector/matrix made in 5_cseqtl_format
# TREC, hap2, ASREC, PHASE made in 3_aseq_format, and aggregate_genes.sh?

## Sample-specific variables 
# RHO 			# cell type proportions matrix
# XX_base 	# baseline covariates, continuous variables centered
# XX_genoPC 	# genotype PCs, centered
# XX_trecPC 	# residual TReC PCs, centered
# XX = cbind(Int = 1,XX_base,XX_genoPC,XX_trecPC)


# Data Directories: ------------------------------------------------------------
script_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/split_scripts"
chrom <- "7" # S{SLURM_TASK_ARRAY}
base_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results"

snp_dir <- paste0(base_dir, "/genotype/LLS/Sept/LLS_Hap/gene_data/chr", chrom)
aseq_dir <- file.path(base_dir, "ASE/per_gene")
ciber_dir <- file.path(base_dir, "rnaseq/ciber/lls") # pp_ciber_combo.csv

meta_dir <- file.path(base_dir, "metadata")
rna_pca_dir <- file.path(base_dir, "metadata")
wgs_pca_dir <- file.path(base_dir, "genotype/LLS/Sept/LLS_Hap/study_pca")

# Note (needed to remove the first two samples as not in aseq data (diff to the outliers), update eventually)

# Load Covariates --------------------------------------------------------------

# XX_base
XX_base <- read.csv(file.path(meta_dir, file = "lls_covar_outlierRM.csv"))

# XX_genoPC (only first 10)	
# genotype PCs, centered
XX_genoPC <- fread(file.path(wgs_pca_dir, "whi_lls.LD.eigenvec"))
# remove outliers
XX_genoPC <- XX_genoPC[XX_genoPC$`#IID` %in% XX_base$topmed_nwdid, ]
# convert ID col to rownames
XX_genoPC <- column_to_rownames(XX_genoPC, var = "#IID")

# XX_trecPC
# residual TReC PCs, centered
# NOTE: Currently just PCs done on raw data, not residuals
XX_trecPC <- fread(file.path(rna_pca_dir, "lls_rnaseq_pcs_outlierRM.csv"))
# Update IDs/order to match geno/covars
XX_trecPC$lls_torid <- XX_trecPC$V1
XX_trecPC <- left_join(XX_base, XX_trecPC, by = "lls_torid") %>%
  column_to_rownames(., var = "topmed_nwdid") %>%
  dplyr::select(., starts_with("PC"))
# Select num of PCs = LLS use 12
XX_trecPC <- (XX_trecPC[,1:12])
colnames(XX_trecPC) <- paste0(colnames(XX_trecPC), "_rna")

# Select relevant covars must include an intercept column, and log-transformed library sizes (gonna test it with the norm multiplied version first)
#TODO: Categorical covars to consider later, smoking, sct (35 individuals), ethnicity
XX <- XX_base %>% dplyr::select(lib.size.norm, neut, bmi)

# Add PCs (checked row order)
XX <- cbind(XX, XX_genoPC, XX_trecPC)

# Impute missing values and centre data
XX_preprocess <- preProcess(XX, method = c("scale", "knnImpute"))
XX <- predict(XX_preprocess, XX)
# Add intercept column
XX <- data.frame(intercept = 1, XX)

# Load RHO (Cell type proportions) ---------------------------------------------

RHO <- read.csv(file.path(ciber_dir, "pp_ciber_combo.csv"))
RHO <- column_to_rownames(RHO, var = "X")
# save.image(file = "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/cseqtl_setup_nov3.RData")

# Load ASeq and SNP data -------------------------------------------------------

# Aseq data made in Aseq format > run_aggregate.sh
# Genotype data made in aggregate, gathers SNPs per gene

# ASEQ data
test_gene <- "ENSG00000002933"
aseq_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE/per_gene/test"
ase_file <- list.files(aseq_dir, pattern = test_gene)
aseq_dat <- fread(file.path(aseq_dir, ase_file))
# ORDER by covar matrix
index <- match(row.names(XX), aseq_dat$sample_id)
# Reorder gene_dat based on the indices
aseq_dat <- aseq_dat[index, ]
# Remove nas
aseq_dat <- aseq_dat %>% filter(!is.na(sample_id))

# Remove missing samples
XX <- XX[rownames(XX) %in% aseq_dat$sample_id, ]

# Make vectors for model
TREC <-  aseq_dat$total
names(TREC) <- aseq_dat$sample_id

ASREC <- aseq_dat$ASREC
names(ASREC) <- aseq_dat$sample_id

PHASE <- rep(1, nrow(aseq_dat))
names(PHASE) <- aseq_dat$sample_id

hap2 <- aseq_dat$hap2
names(hap2) <- aseq_dat$sample_id

# Load genotype data
SNP <- fread(paste0(snp_dir,"/", test_gene, ".txt"))
SNP <- column_to_rownames(SNP, var = "SNP_ID")

## Check orders -----------------------------------------------------------------
# Filter for matching samples and order
SNP <- SNP[, colnames(SNP) %in% row.names(XX)]
index <- match(row.names(XX), colnames(SNP))
SNP <- as.matrix(SNP[ , index])

RHO <- RHO[row.names(RHO) %in% row.names(XX), ]
index <- match(row.names(XX), row.names(RHO))
RHO <- as.matrix(RHO[index, ])

XX <- as.matrix(XX)


CSeQTL_GS(XX=XX, SNP = SNP, PHASE = PHASE,
          TREC = TREC, hap2 = hap2, ASREC = ASREC, RHO = RHO)