# Set up -----------------------------------------------------------------------
library(CSeQTL)
library(data.table)
library(dplyr)
library(smarter)
library(ggplot2)

# An Introduction
vignette(topic = "intro",package = "CSeQTL")

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

### Simulate cell type expression/count data -----------------------------------
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

### Simulate a data object for a gene and SNP ----------------------------------
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
sim$true_SNP # Intitally in 0 1 2 3 code
# Make new SNP col for plotting - recoded as AA AB BA BB
sim$dat$SNP = sim$true_SNP
sim$dat$SNP = factor(sim$dat$SNP,
                     levels = sort(unique(sim$dat$SNP)),
                     labels = c("AA","AB","BA","BB"))
sim$dat$SNP

# TReC vs SNP
ggplot(data = sim$dat,
       mapping = aes(x = SNP,y = log10(total + 1))) +
  geom_violin(aes(fill = SNP)) + geom_jitter(width = 0.25) +
  geom_boxplot(width = 0.1) +
  xlab("Phased Genotype") + ylab("log10(TReC + 1)") +
  theme(legend.position = "right",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12))

## Cell specific testing -------------------------------------------------------

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
      
      sout = CSeQTL_smart(TREC = sim$dat$total, hap2 = sim$dat$hap2,
                          ASREC = sim$dat$total_phased, PHASE = PHASE,
                          SNP = sim$true_SNP, RHO = sim$true_RHO,XX = sim$XX,upPHI = upPHI,
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

### Check input variables: ---------------------------------------

head(sim$dat)
sim$dat$SNP
sim$true_SNP

ols_out = CSeQTL_linearTest(input = sim$dat,XX = sim$XX,
                            RHO = sim$true_RHO,SNP = sim$true_SNP,MARG = TRUE)

str(sim$dat)
# 'data.frame':	300 obs. of  13 variables:
$ total       : num  456 1080 702 797 434 ...
$ LGX1        : num  2340 6468 3903 4532 2206 ...
$ total_phased: int  27 47 30 38 26 51 75 210 33 89 ...
$ hap2        : num  13 17 10 25 16 20 39 113 14 44 ...
$ LBC         : num  16.8 28.6 17.2 22.4 15.5 ...
$ phased      : num  1 1 1 1 1 1 1 1 1 1 ...
$ log_mu      : num  6.21 6.97 6.57 7.34 6.46 ...
$ mu          : num  498 1068 711 1543 637 ...
$ pp          : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
$ phase0      : num  0 0 0 0 0 0 0 0 0 0 ...
$ log_lib_size: num  5.54 5.85 5.75 5.92 5.97 ...
$ geno_col    : chr  "#FF0000BF" "#0000FFBF" "#00FF00BF" "#FF0000BF" ...
$ SNP         : Factor w/ 4 levels "AA","AB","BA",..: 1 4 2 1 4 4 3 3 2 4 ...





