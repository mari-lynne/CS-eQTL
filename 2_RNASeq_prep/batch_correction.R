# COMBATSeq:

# Aims:
# To batch correct and merge RNAseq data from the SCT cohort and LLS for input into CSeQTL pipeline

# Steps:

# Set up =======================================================================

library(data.table)
library(dplyr)
library(tidylog)
library(stringr)
library(viridis)
library(ggplot2)
library(janitor)
library(smarter)
library(patchwork)

plot_dir <- "~/Documents/CSeQTL/data/plots"
work_dir <- "~/Documents/CSeQTL/data/ciber_ase"
meta_dir <- "~/Documents/whi_sct/rna/meta"

setwd(work_dir)
link_sct <- read.csv(file = file.path(meta_dir,"sct_all_covars_Jul23.csv"))
link_lls <- read.csv(file = "~/Documents/CSeQTL/data/meta/lls/LLS_ID_pheno_08-23.csv")

## Load TPM data ---------------------------------------------------------------

sct_tpm <- fread("~/Documents/whi_sct/rna/results/ciber/R/mixture.txt")
# Significantly more genes in sct_tpm, had not filtered for low expressing genes
# which these might be low in mono/neuts, makes sense that they were exluded.

lls_tpm <- fread(file.path(work_dir, "mixture_qc.txt"))

## Update both to use subject ID from link file as ids -------------------------
# Run/plot PCA of data merged without BC =======================================
# Run COMBATSeq ================================================================
# Run/plot PCA of merged data ==================================================
# Check individual genes, see if there is a way of using shared samples

