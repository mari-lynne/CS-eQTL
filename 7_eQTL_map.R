# devtools::install_github("pllittle/smarter")
# devtools::install_github("pllittle/CSeQTL")
library("smarter")
library("CSeQTL")

wd <- "/home/mari/Documents/CSeQTL/data/ciber_ase"

# Input ------------------------------------------------------------------------

# TReC data
# ASReC data (haplotype 1 and 2 read counts)
# Phased genotypes
# Sample covariates:
#   observed confounders,
# genotype principal components (PCs),
# latent batch effects (derived from residual TReC PCs)
# Cell type proportions

# Get genotype PCs for these participants
# Calculate TREC Pcs in cibersort script

# Sample-specific variables
RHO 			# cell type proportions matrix
XX_base 	# baseline covariates, continuous variables centered
XX_genoPC 	# genotype PCs, centered
XX_trecPC 	# residual TReC PCs, centered
XX = cbind(Int = 1,XX_base,XX_genoPC,XX_trecPC)

# Gene and sample variables
TREC 	# TReC vector (Per gene!)
SNP	# phased genotype vector
hap2	# 2nd haplotype counts
ASREC 	# total haplotype counts = hap1 + hap2
PHASE 	# Indicator vector of whether or not to use haplotype counts

# Tuning arguments
trim 		# TRUE for trimmed analysis, FALSE for untrimmed
thres_TRIM 	# if trim = TRUE, the Cooks' distance cutoff to trim sample TReCs
ncores 	# number of threads/cores to parallelize the loop, improve runtime

# ASREC-related cutoffs to satisfy to use allele-specific read counts
numAS 		# cutoff to determine if a sample has sufficient read counts
numASn 	# cutoff of samples having ASREC >= numAS
numAS_het 	# minimum number for sum(ASREC >= numAS & SNP %in% c(1,2))
cistrans 	# p-value cutoff on cis/trans testing to determine which p-value 
# (TReC-only or TReC+ASReC) to report



# 1) Sample-specific variables -------------------------------------------------
# See 6_covar_pca.R

# XX_base baseline covariates, continuous variables centered
# XX_genoPC; genotype PCs, centered
# XX_trecPC 	# residual TReC PCs, centered

covar <- read.csv(file = "~/Documents/whi_sca/rna/meta/sct_cseqtl_test_covars.csv")
row.names(covar) <- covar$nwgc_sample_id
# Filter for test samples and relevant covariate columns
samples <- c("941357", "941437")
covar2 <- covar[covar$nwgc_sample_id %in% samples , c(2,3,16:31)]  # check order is the same
# group col acting as y intercept
XX <- as.matrix(covar2)
row.names(XX) <- covar$X

# RHO cell type proportions matrix
RHO <- read.csv(file = "~/Documents/CSeQTL/data/ciber_ase/pp_hat_ciber_test.csv")
row.names(RHO) <- RHO$X
RHO <- as.matrix(RHO[,-1]) 
 	

# 2) Gene and sample variables -------------------------------------------------
# Eventually run for manifest in SLURM
# args = commandArgs(trailingOnly=TRUE)
# manifest <- read.csv(args[1] , stringsAsFactors = F) 
# thisRun <- manifest[Sys.getenv('SLURM_ARRAY_TASK_ID'),] 

thisRun <- c("/NWD543548")

### TREC 
# Vector containing total read count of a gene per sample
# TREC <- read.csv(file = "~/Documents/CSeQTL/data/ciber_ase/filtered_genes_test.csv")
trec_ase <-fread(paste0(wd, thisRun, "-output.trecase.txt"))
TREC <- trec_ase[,2]

# hap2, ASREC, PHASE
hap2 <- trec_ase[,4]
ASREC <- trec_ase[,3] + trec_ase[,4]

PHASE <- rep(1, nrow(trec_ase))

# Pick test gene and sample

### SNP,  -----------------------------------------------------
# See 4_test_ase.R
# cts = cbind(ct,ct1,ct2,ctN) # trec, hap1, hap2, hapN 
# saved as thisRun$topmed_nwdid, "-output.trecase.txt" in ase dir
# SNP	# phased genotype vector in trecase dir (see step 4)



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
gs_out = CSeQTL_GS(XX = XX,TREC = TREC, SNP = SNP,hap2 = hap2,
                   ASREC = ASREC, PHASE = PHASE, RHO = RHO, trim = trim,
                   thres_TRIM = 20, numAS = 5, numASn = 5, numAS_het = 5,
                   cistrans = 0.01, ncores = ncores, show = TRUE)
