# devtools::install_github("pllittle/smarter")
# devtools::install_github("pllittle/CSeQTL")
library("smarter")
library("CSeQTL")

# Get genotype PCs for these participants
# Calculate TREC Pcs in cibersort script

# 1) Sample-specific variables -------------------------------------------------
# See 6_covar_pca.R

# XX_base baseline covariates, continuous variables centered
# Also contains XX_genoPC; genotype PCs, centered, XX_trecPC 	# residual TReC PCs, centered
covar <- read.csv(file = "~/Documents/whi_sca/rna/meta/sct_cseqtl_test_covars.csv")
# filter for test samples and relevant covar columns 
covar2 <- covar[covar$nwgc_sample_id %in% c("941357", "941437") , c(2,3,16:31)]  # check order is the same
# group col accting as y intercept

XX <- as.matrix(covar2)

# RHO cell type proportions matrix
RHO <- read.csv(file = "~/Documents/CSeQTL/data/ciber_ase/pp_hat_ciber_test.csv")
row.names(RHO) <- RHO$X
RHO <- as.matrix(RHO[,-1]) 
 	

# 2) Gene and sample variables -------------------------------------------------
 	
### TREC 
# Vector containing total read count of a gene per sample
TREC <- read.csv(file = "~/Documents/CSeQTL/data/ciber_ase/filtered_genes_test.csv")
# Pick test gene and samples
TREC <- TREC[TREC$hgnc_symbol == "SLC25A5", colnames(TREC) %in% c("X941357", "X941437")]

## SNP, hap2, ASREC, PHASE -----------------------------------------------------
# See 4_test_ase.R
# cts = cbind(ct,ct1,ct2,ctN) # trec, hap1, hap2, hapN 
# saved as thisRun$topmed_nwdid, "-output.trecase.txt" in ase dir
# SNP	# phased genotype vector in trecase dir (see step 4)

# sample = thisRun$topmed_nwdid test = "941357"

cts <- fread("~/Documents/CSeQTL/data/ciber_ase/NWD543548-output.trecase.txt")

SNP <-
hap2	# 2nd haplotype counts
ASREC 	# total haplotype counts = hap1 + hap2
PHASE 	# Indicator vector of whether or not to use haplotype counts

# 3) Tuning arguments ----------------------------------------------------------
trim 		# TRUE for trimmed analysis, FALSE for untrimmed
thres_TRIM 	# if trim = TRUE, the Cooks' distance cutoff to trim sample TReCs
ncores 	# number of threads/cores to parallelize the loop, improve runtime

# 4) ASREC-related cutoffs to satisfy to use allele-specific read counts -------
numAS 		# cutoff to determine if a sample has sufficient read counts
numASn 	# cutoff of samples having ASREC >= numAS
numAS_het 	# minimum number for sum(ASREC >= numAS & SNP %in% c(1,2))
cistrans 	# p-value cutoff on cis/trans testing to determine which p-value 
# (TReC-only or TReC+ASReC) to report


# Run eQTL mapping -------------------------------------------------------------
gs_out = CSeQTL_GS(XX = XX,TREC = TREC,SNP = SNP,hap2 = hap2,
                   ASREC = ASREC, PHASE = PHASE, RHO = RHO, trim = trim,
                   thres_TRIM = 20, numAS = 5, numASn = 5, numAS_het = 5,
                   cistrans = 0.01, ncores = ncores, show = TRUE)
