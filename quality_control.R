# PCA and latent batch effect:

# Aims
# Run PCA of expression data to account for population structure/variance in datra when modelling eQTLs
# Check expression PCs against genotype PCs, cell type proprotions and age + bmi covariates

library(data.table)
library("FactoMineR")
library("factoextra")
library("PCAForQTL")
library(ggplot2)
library(viridis)
library(patchwork)
library(dplyr)
library(tidylog)
library(colorspace)
library(forcats)
library(viridis)

setwd("~/Documents/CSeQTL/data/pca")
dge <- readRDS(file = "dge_black.rds")

# Calculate PCA of RNAseq TReC data --------------------------------------------

exprs <- t(dge$counts)

prcompResult <-prcomp(exprs,center=TRUE,scale.=TRUE)
PCs <-prcompResult$x
dim(PCs)
resultRunElbow <-PCAForQTL::runElbow(prcompResult=prcompResult)
print(resultRunElbow)

# Scree Plots ------------------------------------------------------------------
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow"),values=c(K_elbow),
                         titleText="RNA-Seq PCA")
ggsave(filename = "scree.png", device = "png")

# Filter PCs
eigenval <- prcompResult$sdev[1:30]^2

pve <- data.frame(
  PC = 1:length(eigenval),
  pve = (eigenval / sum(eigenval) * 100)
)
pve$PC <- as.factor(pve$PC)

PC_Palette <- viridis(30)

pve_plot <- ggplot(pve, aes(PC, pve, fill = PC)) +
  geom_bar(stat = "identity") +
  ylab("Proportion of Variance (%)\n") +
  scale_fill_manual(values = PC_Palette) +
  theme(legend.position = "none") + geom_vline(xintercept = 11.5, linetype = "dashed", color = "black") +
  ggtitle("RNAseq TReC PCA")

pc_top <-PCs[,1:11] 

write.table(pc_top, file = "pc_top_rnaseq.txt", quote = F, col.names = T, row.names = T)

# Extra Covariates -------------------------------------------------------------
# Sample covars: age, bmi, ethnicity when required
# Cell proportion covars: cibersort output
# Genotype covars:

meta_dir <- c("~/Documents/whi_sca/rna/meta")
samp_covar <- read.csv(file.path(meta_dir, "sct_all_covars_Jul23.csv"))
ciber <- read.csv(file = "~/Documents/whi_sca/rna/results/ciber/R/pp_hat_ciber_merged.csv")
geno <-   

  # Join
samp_covar <- select(samp_covar, group, lib.size, rnaseq_ids, nwgc_sample_id, plate, age, bmi_t0) 
ciber <- rename(ciber, rnaseq_ids = X)
covar_all <- left_join(samp_covar, ciber, by = "rnaseq_ids")
pc_top <- as.data.frame(pc_top)
pc_top$rnaseq_ids <- row.names(pc_top)
covar_all <- left_join(covar_all, pc_top)
covar_sub <- select(covar_all, -c("group", "lib.size", "nwgc_sample_id", "plate", "rnaseq_ids"))

# Check latent batch effects ---------------------------------------------------

# Fit a linear model to expression data (TreC) + all covariates
# Step 1: Log-transformed TReC as response variable # where rows represent genes and columns represent samples
# Log-transformed TReC calculation
# log(TReC of gene j of sample i/read-depth of sample i) as response variable,
log_trec <- log(dge$counts / colSums(dge$counts))
log_trec <- t(log_trec)

# Step 2: eQTL mapping with linear model
eqtl_model <- lm(log_trec ~ ., data = covar_sub)

# Calculate the residuals
residuals <- residuals(eqtl_model)

# Perform PCA on the residuals
res_pca <- prcomp(residuals)

# Plot the PCA of residuals
plot(res_pca$x[,1], res_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of Residuals")

# Plot original PCA:
plot(prcompResult$x[,1], res_pca$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of TReC")

# GGPLOT --------
# Plot the PCA of residuals
residuals_pca <- data.frame(PC1 = res_pca$x[, 1], PC2 = res_pca$x[, 2], plate = as.factor(covar_all$plate))

pca_resid <-
  ggplot(residuals_pca, aes(x = PC1, y = PC2, color = plate)) +
  geom_point() +
  labs(x = "PC1", y = "PC2", title = "PCA of Residuals")

# Plot the original PCA with color-coded points by covar_sub$plate
original_pca <- data.frame(PC1 = prcompResult$x[, 1], PC2 = prcompResult$x[, 2], plate = as.factor(covar_all$plate))
pca_og <- 
  ggplot(original_pca, aes(x = PC1, y = PC2, color = plate)) +
  geom_point() +
  labs(x = "PC1", y = "PC2", title = "PCA of TReC")

# PVE plot resid ----------------
eigenval <- res_pca$sdev[1:30]^2

pve <- data.frame(
  PC = 1:length(eigenval),
  pve = (eigenval / sum(eigenval) * 100)
)
pve$PC <- as.factor(pve$PC)

PC_Palette <- viridis(length(eigenval))

pve_resid <- ggplot(pve, aes(PC, pve, fill = PC)) +
  geom_bar(stat = "identity") +
  ylab("Proportion of Variance (%)\n") +
  scale_fill_manual(values = PC_Palette) +
  theme(legend.position = "none") +
  ggtitle("RNAseq TReC - Residual PCA")

# Collate ----------------------------
(pca_og|pca_resid)/(pve_plot|pve_resid) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")

og_pca <- data.frame(PC10 = prcompResult$x[, 10], PC11 = prcompResult$x[, 11], plate = as.factor(covar_all$plate))
ggplot(og_pca, aes(x = PC10, y = PC11, color = plate)) +
  geom_point() +
  labs(x = "PC1", y = "PC2", title = "PCA of TReC")

# Pvale enrich ------------------
# Step 3: Assess covariate significance using p-value histogram
baseline_model <- lm(log_trec ~ 1, data = covar_sub)  # Baseline model without SNP genotype

# Function to compute p-value histogram
compute_pvalue_histogram <- function(model, covariate) {
  p_values <- summary(model)$coefficients[rownames(summary(model)$coefficients) == covariate, "Pr(>|t|)"]
  hist(p_values, breaks = 20, main = paste("Histogram of p-values for", covariate), xlab = "p-values")
}

# Compute and plot p-value histogram for each covariate
covariates <- colnames(covar_sub)[-1]  # Exclude the response variable column
for (covariate in covariates) {
  compute_pvalue_histogram(baseline_model, covariate)
}

# Function to compute p-value histogram
compute_pvalue_histogram <- function(model, covariate) {
  p_values <- as.numeric(summary(model)$coefficients[rownames(summary(model)$coefficients) == covariate, "Pr(>|t|)"])
  hist(p_values, breaks = 20, main = paste("Histogram of p-values for", covariate), xlab = "p-values")
}

# Compute and plot p-value histogram for each covariate
covariates <- colnames(covar_sub)[-1]  # Exclude the response variable column
for (covariate in covariates) {
  compute_pvalue_histogram(baseline_model, covariate)
}

# Compute and plot p-value histogram for each covariate
covariates <- colnames(covar_sub)[-1]  # Exclude the response variable column
par(mfrow = c(ceiling(length(covariates) / 2), 2))  # Set up multiple plots in a grid

for (i in 1:length(covariates)) {
  compute_pvalue_histogram(baseline_model, covariates[i])
}

par(mfrow = c(1, 1)) 


# Known covars correlation ---------------------------------------
# Add genotype PCs

knownCovariates<- select(dge$samples, "age", "bmi_t0")
dataCovariates[,c(1:5,66:68)] #368*8. 368 samples, 8 known covariates.
identical(rownames(knownCovariates),rownames(expr))
covars <- cbind(pheno, pc_top)



