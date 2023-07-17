#!/usr/bin/env Rscript

## Args and manifest -----------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)

# Test if there is not one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied, the manifest file name.n", call.=FALSE)
}

# Read in the specific manifest you indicated, then select the 
# row you want to use as arguments.
manifest <- read.csv(args[1] , stringsAsFactors = F) 
thisRun <- manifest[Sys.getenv('SLURM_ARRAY_TASK_ID'),] # subset array by index of slurm tasks, e.g 3rd task subsets 3rd row of the manifest

### Packages ------------------------------------------------------------------

# Installed extra packages to rhino home dir https://sciwiki.fredhutch.org/rModules/
library(stringr)
library(data.table)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(CSeQTL)
library(asSeq)

##  1 Input output files -------------------------------------------------------

# Input Directories
# main_dir = c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/")
main_dir = c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/test_june1/")
bam_dir = c("/fh/scratch/delete90/kooperberg_c/sct_rnaseq/release_files/")

# # Output Dirs

# rna_dir = paste0(main_dir, "rnaseq/")
# gen_dir =  paste0(main_dir, "genotype/split/sub/")
# ase_dir = paste0(main_dir, "ASE/")

rna_dir = main_dir
gen_dir =  main_dir
ase_dir = main_dir

# File names
bam_file = paste0(bam_dir, thisRun$bam_file)
bam_id = thisRun$nwgc_sample_id
bam_filt_fn = paste0(rna_dir, bam_id, "-output.filtered.asSeq.bam")
sortQ_bam_fn = paste0(rna_dir, bam_id, "-output.filtered.asSeq.sortQ")

gtf_rds_fn = c("exon_by_genes.rds") # Made in 1_get-exon-info
genes = readRDS(gtf_rds_fn)
het_snp_fn = paste0(gen_dir, thisRun$topmed_nwdid, "_hap.txt") 

PE = TRUE # Paired end reads = TRUE

# 2) Get total read count per gene ---------------------------------------------

bamfile = Rsamtools::BamFileList(bam_filt_fn, yieldSize = 1000000)

se = GenomicAlignments::summarizeOverlaps(features = genes,
                                          reads = bamfile, mode = "Union",singleEnd = !PE,
                                          ignore.strand = TRUE, fragments = PE)

ct = as.data.frame(SummarizedExperiment::assay(se)) # reads

print("Total Read Count DONE")
date()

# Filter reads by Qname
sortBam(file = bam_filt_fn, destination = sortQ_bam_fn, byQname=TRUE)

print("Sorted BAMs by Qname DONE")
date()

# 3) Get allele-specific read counts -------------------------------------------

# Requires tab delim file of heterozygous snps per sample: het_snp_fn
# - No header - Columns: chr, position, hap1 allele, hap2 allele


## 3a) Extract ASE reads --------------------------------------------------------

# Function extracts reads determined to have allele-specific expression
# A sequence read is counted as allele-specific only if:
# It harbors one or more hetrozygous genetic markers
# 2) Genotypes are consistent with one (not both) of the haplotypes specified in snpList

asSeq::extractAsReads(input = sortQ_bam_fn,
                      snpList = het_snp_fn, min.avgQ = 20, min.snpQ = 20)

print("ASE counting DONE")
date()

## 3b) Count allele specific read counts ASRC ----------------------------------

se1 = GenomicAlignments::summarizeOverlaps(features = genes,
                                           reads = sprintf("%s_hap1.bam", sortQ_bam_fn), mode = "Union",
                                           singleEnd = !PE,ignore.strand = TRUE,fragments = PE)
se2 = GenomicAlignments::summarizeOverlaps(features = genes,
                                           reads = sprintf("%s_hap2.bam", sortQ_bam_fn), mode = "Union",
                                           singleEnd = !PE,ignore.strand = TRUE,fragments = PE)
seN = GenomicAlignments::summarizeOverlaps(features = genes,
                                           reads = sprintf("%s_hapN.bam", sortQ_bam_fn), mode = "Union",
                                           singleEnd = !PE,ignore.strand = TRUE,fragments = PE)

ct1 = as.data.frame(SummarizedExperiment::assay(se1))
ct2 = as.data.frame(SummarizedExperiment::assay(se2))
ctN = as.data.frame(SummarizedExperiment::assay(seN))
cts = cbind(ct,ct1,ct2,ctN) # trec, hap1, hap2, hapN
dim(cts); cts[1:2,]

out_fn = paste0(thisRun$topmed_nwdid, "-output.trecase.txt") # output to main results dir
write.table(cts, file = paste0(ase_dir, out_fn), quote = FALSE,
            sep = "\t", eol = "\n")

print("ASE mapping DONE!")
date()


# Plot summary histograms ------------------------------------------------------

# load test data (2 indiviuals for now)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

cts <- fread("~/Documents/CSeQTL/data/ciber_ase/NWD543548-output.trecase.txt")
cts2 <- fread("~/Documents/CSeQTL/data/ciber_ase/NWD543548-output.trecase.txt")

# List all the files in your directory with the pattern "-output.trecase.txt"
setwd("~/Documents/CSeQTL/data/ciber_ase/ase")

# Create an empty data frame to store the results
summary_data <- data.frame()

# Get a list of files in the directory with the specified pattern
file_list <- list.files(pattern = "-output.trecase.txt")
for (file in file_list) {
  # Extract the sample ID from the file name
  sample_id <- str_extract(file, "(?<=^)[^\\.]+")
  
  # Read the file into a data frame
  data <- fread(file, header = TRUE, col.names = c("V1", paste0(sample_id, c("-output.filtered.asSeq.bam",
                                                                             "-output.filtered.asSeq.sortQ_hap1.bam",
                                                                             "-output.filtered.asSeq.sortQ_hap2.bam",
                                                                             "-output.filtered.asSeq.sortQ_hapN.bam"))))
  # Rename the columns
  colnames(data) <- c("gene_id", "trec", "hap1", "hap2", "hapN")
  
  # Calculate column sums
  colsums <- colSums(data[, c("trec", "hap1", "hap2", "hapN")], na.rm = TRUE)
  
  # Create a summary row with the sample ID and colsum values
  summary_row <- data.frame(sample_ID = sample_id,
                            trec_colsum = colsums["trec"],
                            hap1_colsum = colsums["hap1"],
                            hap2_colsum = colsums["hap2"],
                            hapN_colsum = colsums["hapN"])
  
  # Append the summary row to the summary_data data frame
  summary_data <- bind_rows(summary_data, summary_row)
}


# Step 1: Calculate mean and standard deviation
# Step 1: Calculate mean and standard deviation
col_means <- colMeans(summary_data[, -1])  # Exclude the first column (sample_ID)
col_sds <- apply(summary_data[, -1], 2, sd)  # Standard deviation of each column

# Step 2: Generate random values within 1 SD of the mean for each column
dummy_samples <- data.frame(sample_ID = paste0("Dummy", 1:10))  # Create a data frame for dummy samples

for (col_name in colnames(summary_data)[-1]) {
  mean_val <- col_means[col_name]
  sd_val <- col_sds[col_name]
  
  lower_limit <- mean_val - sd_val
  upper_limit <- mean_val + sd_val
  
  dummy_values <- rnorm(n = 10, mean = mean_val, sd = sd_val)
  dummy_values[dummy_values < lower_limit] <- lower_limit
  dummy_values[dummy_values > upper_limit] <- upper_limit
  
  dummy_samples[[col_name]] <- dummy_values
}

# Step 3: Combine summary_data and dummy_samples
combined_data <- rbind(summary_data, dummy_samples)

# Print the combined data
print(combined_data)


# Print the combined data
print(combined_data)

og <- summary_data
summary_data <- read.csv(file = "~/Documents/CSeQTL/data/ciber_ase/ase/hap_plots2.csv")
library(scales)

hap1 <-   ggplot(summary_data, aes(x = trec_colsum, y = hap1_colsum)) +
  geom_point() +
  labs(x = "trec", y = "hap1") +
  ggtitle("Scatter Plot: hap1 vs trec") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) + 
  theme_minimal()

hap2 <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hap2_colsum)) +
  geom_point() +
  labs(x = "trec", y = "hap2") +
  ggtitle("Scatter Plot: hap2 vs trec") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) + 
  theme_minimal()

hapN <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hapN_colsum)) +
  geom_point() +
  labs(x = "trec", y = "hapN") +
  ggtitle("Scatter Plot: hapN vs trec") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) + 
  theme_minimal()

library(patchwork)

(hap1 | hap2 |hapN) + plot_annotation(tag_levels = "a")

write.csv(summary_data, file = "hap_plots.csv", row.names = F)

read.csv(file = "hap_plots.csv")


# Test -----

data <- fread("~/Documents/CSeQTL/data/ciber_ase/NWD543548-output.trecase.txt")

data <- data.frame(gene_id = c(1, 2, 3),
                   sample_ID_trec = c(10, 20, 30),
                   sample_ID_hap1 = c(5, 15, 25),
                   sample_ID_hap2 = c(7, 17, 27),
                   sample_ID_hapN = c(8, 18, 28))

# Rename the columns
colnames(data) <- c("gene_id", "trec", "hap1", "hap2", "hapN")

# Calculate column sums
colsums <- colSums(data[, c("trec", "hap1", "hap2", "hapN")], na.rm = TRUE)

# Create a summary data frame
summary_data <- data.frame(sample_ID = "sample_ID",
                           trec_colsum = colsums["trec"],
                           hap1_colsum = colsums["hap1"],
                           hap2_colsum = colsums["hap2"],
                           hapN_colsum = colsums["hapN"]



ggplot(combined_data, aes(x = trec, y = sample_ID-hap1)) +
  geom_point() +
  labs(title = "Haplotype 1",
       x = "Total RNA-seq fragments",
       y = "Haplotype 1 fragments")
