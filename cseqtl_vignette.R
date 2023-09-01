library(CSeQTL)
library(data.table)
library(dplyr)
library(smarter)
library(ggplot2)

# An Introduction
vignette(topic = "intro",package = "CSeQTL")

# Directories and input

gen_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/test_aug/ASE/genes")
test_gene <- c("ENSG00000000457.14_counts.txt")

# Per gene, and per SNP all samples data: --------------------------------------
TREC 	# TReC vector (Per gene!)
SNP	# phased genotype vector
hap2	# 2nd haplotype counts
ASREC 	# total haplotype counts = hap1 + hap2
PHASE 	# Indicator vector of whether or not to use haplotype counts

# make in aggregate_genes.sh
gene_dat <- fread(file.path(gen_dir, test_gene))
TREC <- gene_dat$total
SNP <- rep(1, nrow(gene_dat)) # Replace with SNP of interest phase info
ASREC <- gene_dat$ASREC
PHASE <- rep(1, nrow(gene_dat))

# Get haplotype info for 
# Simulate a data object for a gene and SNP
SNP <- sim$true_SNP

# Sample-specific variables -------------------------------------------------
# See 6_covar_pca.R

RHO 			# cell type proportions matrix
XX_base 	# baseline covariates, continuous variables centered
XX_genoPC 	# genotype PCs, centered
XX_trecPC 	# residual TReC PCs, centered

covar <- read.csv(file = "~/Documents/whi_sca/rna/meta/sct_cseqtl_test_covars.csv")
row.names(covar) <- covar$nwgc_sample_id
covar <- covar[1:20, ]
# Sample covars
XX_base <- covar %>% select(age, bmi_t0)
XX_base <- scale(XX_base)
XX_base <- cbind(covar$group, XX_base) # group col acting as y intercept

# RNAseq covars
XX_trecPC <- covar %>% select(starts_with("PC"))

# Geno covars
geno_pc <- fread(file = "~/Documents/CSeQTL/data/ciber_ase/ase/whi_sub_concat.LD.eigenvec")
XX_genoPC <- geno_pc[1:20, 3:10] # use first 7 for now, prob with 20 samples will reduce rnaseq covars too
XX_genoPC <- scale(XX_genoPC)

# RHO cell type proportions matrix
RHO <- read.csv(file = "~/Documents/whi_sca/rna/results/ciber/R/pp_hat_ciber_merged.csv")
RHO <- RHO[RHO$X %in% covar$X, ]
RHO <- RHO[match(covar$X, RHO$X),] # Match up row orders
row.names(RHO) <- covar$nwgc_sample_id
RHO <- RHO[, -1] # remove ID col
XX = cbind(XX_base,XX_genoPC,XX_trecPC)

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
gs_out = CSeQTL_GS(XX = XX,TREC = TREC, SNP = SNP, hap2 = hap2,
                   ASREC = ASREC, PHASE = PHASE, RHO = RHO, trim = trim,
                   thres_TRIM = 20, numAS = 5, numASn = 5, numAS_het = 5,
                   cistrans = 0.01, ncores = ncores, show = TRUE)



# Vignette ---------------------------------------------------------------------
vignette(topic = "intro",package = "CSeQTL")

req_packs = c("ggplot2","smarter","CSeQTL")
for(pack in req_packs){
  library(package = pack,character.only = TRUE)
}

# List package's exported functions
ls("package:CSeQTL")
# Fix seed
set.seed(2)

# Simulate cell type expression/count data -------------------------------------
# sample size
NN = 300

# fold-change between q-th and 1st cell type
true_KAPPA  = c(1,3,1)

# eQTL effect size per cell type, 
#   fold change between B and A allele
true_ETA    = c(1,1,1)

# cis/trans effect size
true_ALPHA  = c(1,1,1)

# count number of cell types
QQ = length(true_KAPPA)

# TReC model overdispersion
true_PHI = 0.1

# ASReC model overdispersion
true_PSI = 0.05

# Simulate cell type proportions
tRHO = gen_true_RHO(wRHO = 1,NN = NN,QQ = QQ)
plot_RHO(RHO = tRHO)


# Simulate a data object for a gene and SNP ------------------------------------
sim = CSeQTL_dataGen(NN = NN,MAF = 0.3,true_BETA0 = log(1000),
                     true_KAPPA = true_KAPPA,true_ETA = true_ETA,true_PHI = true_PHI,
                     true_PSI = true_PSI,prob_phased = 0.05,true_ALPHA = true_ALPHA,
                     RHO = tRHO,cnfSNP = TRUE,show = FALSE)
# List contains input vars and datasets for downstream CSeQTL functions

names(sim)
#>  [1] "true_PARS" "true_SNP"  "dat"       "XX"        "true_RHO"  "QQ"       
#>  [7] "GI"        "np"        "I_np"      "iBETA"     "iPHI"      "MU"       
#> [13] "true_OF"   "vname"

data <- sim$dat

# TReC, ASReC, haplotype 2 counts
sim$dat[1:5,]

#>   total     LGX1 total_phased hap2      LBC phased   log_mu        mu  pp
#> 1   456 2339.837           27   13 16.81415      1 6.209945  497.6740 0.5
#> 2  1080 6467.905           47   17 28.63941      1 6.973643 1068.1067 0.5
#> 3   702 3903.057           30   10 17.21821      1 6.566389  710.7982 0.5
#>   phase0 log_lib_size  geno_col
#> 1      0     5.538577 #FF0000BF
#> 2      0     5.850341 #0000FFBF
#> 3      0     5.751526 #00FF00BF

# Batch covariates including the intercept
sim$XX[1:5,]
#>      X1         X2 X3         X4         X5
#> [1,]  1 -1.7890676  1 -0.5523594 1.35690857
#> [2,]  1 -0.6066932  0  1.5950620 1.06844120
#> [3,]  1 -0.9814495  0 -0.9546107 0.06637533

sim$dat$SNP = sim$true_SNP
sim$dat$SNP = factor(sim$dat$SNP,
                     levels = sort(unique(sim$dat$SNP)),
                     labels = c("AA","AB","BA","BB"))

# TReC vs SNP
ggplot(data = sim$dat,
       mapping = aes(x = SNP,y = log10(total + 1))) +
  geom_violin(aes(fill = SNP)) + geom_jitter(width = 0.25) +
  geom_boxplot(width = 0.1) +
  xlab("Phased Genotype") + ylab("log10(TReC + 1)") +
  theme(legend.position = "right",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))

# Cell specific testing --------------------------------------------------------

cistrans_thres = 0.01
res = c()

for(trec_only in c(TRUE,FALSE)){
  for(neg_binom in c(TRUE,FALSE)){
    for(beta_binom in c(TRUE,FALSE)){
      
      cat(sprintf("%s: trec_only = %s, neg_binom = %s, beta_binom = %s ...\n",
                  date(),trec_only,neg_binom,beta_binom))
      
      PHASE = sim$dat$phased * ifelse(trec_only,0,1)
      upPHI = ifelse(neg_binom,1,0)
      upPSI = ifelse(beta_binom,1,0)
      upKAPPA = rep(1,QQ)
      upETA = upKAPPA
      upALPHA = upETA
      
      sout = CSeQTL_smart(TREC = sim$dat$total,hap2 = sim$dat$hap2,
                          ASREC = sim$dat$total_phased,PHASE = PHASE,
                          SNP = sim$true_SNP,RHO = sim$true_RHO,XX = sim$XX,upPHI = upPHI,
                          upKAPPA = upKAPPA,upETA = upETA,upPSI = upPSI,upALPHA = upALPHA,
                          iFullModel = FALSE,trim = FALSE,thres_TRIM = 10,
                          hypotest = TRUE,swap = FALSE,numAS = 5,numASn = 5,
                          numAS_het = 5,cistrans_thres = cistrans_thres)
      
      res = rbind(res,smart_df(Model = ifelse(trec_only,"TReC-only","TReCASE"),
                               TReC_Dist = ifelse(neg_binom,"Negative Binomial","Poisson"),
                               ASReC_Dist = ifelse(beta_binom,"Beta-Binomial","Binomial"),
                               sout$res))
      rm(sout)
      
    }}}

# Check

sout = CSeQTL_smart(TREC = sim$dat$total,hap2 = sim$dat$hap2,
                    ASREC = sim$dat$total_phased,PHASE = PHASE,
                    SNP = sim$true_SNP,RHO = sim$true_RHO,XX = sim$XX,upPHI = upPHI,
                    upKAPPA = upKAPPA,upETA = upETA,upPSI = upPSI,upALPHA = upALPHA,
                    iFullModel = FALSE,trim = FALSE,thres_TRIM = 10,
                    hypotest = TRUE,swap = FALSE,numAS = 5,numASn = 5,
                    numAS_het = 5,cistrans_thres = cistrans_thres)
