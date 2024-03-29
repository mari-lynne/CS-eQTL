# Aims:

# Examine haplotype specific and total read counts generated by previous ASeq2 code

# Set up: ----------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(patchwork)

source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/functions.R")
covar_name <- "LLS_ID_pheno_08-23.csv"
meta_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata"
wd <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE"
plot_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots/LLS/ASeq"

## Format aseq data ------------------------------------------------------------

# List all the files in your directory with the pattern "-output.trecase.txt"
setwd(wd)

# Create an empty data frame to store the results
summary_data <- data.frame()

# Get a list of files in the directory with the specified pattern
file_list <- list.files(pattern = "-output.trecase.txt")
for (file in file_list) {
  # Extract the sample ID from the file name
  sample_id <- str_extract(file, "(?<=^)[^\\.]+")
  
  # Read the file into a data frame
  data <-
    fread(file, header = TRUE, col.names = c("V1", paste0(
      sample_id,
      c(
        "-output.filtered.asSeq.bam",
        "-output.filtered.asSeq.sortQ_hap1.bam",
        "-output.filtered.asSeq.sortQ_hap2.bam",
        "-output.filtered.asSeq.sortQ_hapN.bam"
      )
    )))
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
col_means <- colMeans(summary_data[, -1])  # Exclude the first column (sample_ID)
col_sds <- apply(summary_data[, -1], 2, sd)  # Standard deviation of each column

## Add in phenotype data -------------------------------------------------------

covar <- read_omic(name = covar_name, dir = meta_dir)
covar2 <- select(covar, topmed_nwdid, ethnic.x)

summary_data$topmed_nwdid <- sub("-output$", "", summary_data$sample_ID)
summary_data <- inner_join(summary_data, covar2, by = "topmed_nwdid")
summary_data$ethnicity <- as.factor(summary_data$ethnic.x)

summary_data$ethnicity <- recode(summary_data$ethnicity,
                                 '3' = 'A.A',
                                 '4' = 'His',
                                 '5' = 'E.A')

# change order so hispanics are plotted
summary_data$ethnicity <- fct_relevel(summary_data$ethnicity, "His", after = Inf)
summary_data <- summary_data[order(summary_data$ethnicity), ]

pal <- c("#CC4678D9", "#FCCE25D9",  "#900DA4D9")

# Filter outliers both from rnaseq analysis and trec plots
outliers <- c("NWD196624", "NWD774712", "NWD148594")
summary_data <- filter(summary_data, topmed_nwdid %!in% outliers)

# Plot summary histograms ------------------------------------------------------

hap1 <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hap1_colsum, color=ethnicity)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  labs(x = "TReC", y = "Hap1") +
  ggtitle("Hap1 vs TReC") +
  theme_minimal()

hap2 <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hap2_colsum, color=ethnicity)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  labs(x = "TReC", y = "Hap2") +
  ggtitle("Hap2 vs TReC") + 
  theme_minimal() 

hapN <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hapN_colsum, color=ethnicity)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = pal) +
  scale_y_continuous(trans='log10') +
  scale_x_continuous(trans='log10') +
  labs(x = "TReC", y = "HapN") +
  ggtitle("HapN vs TReC") +
  theme_minimal()

### Main plot ------------------------------------------------------------------

(hap1 | hap2 |hapN) + plot_annotation(tag_levels = "a") + plot_layout(guides = "collect")

ggsave(filename = "hap_trec_ethnic_scatters.png", path = plot_dir, width = 14.5, height = 5.7, units = c("in"))

### HAP1 v HAP2 ----------------------------------------------------------------
# Use the sample() function to randomly select rows
samp_data <- summary_data[sample(nrow(summary_data), 
                                 size = nrow(summary_data) * 0.4), ]

samp_data %>%
  ggplot(aes(x = hap1_colsum, y = hap2_colsum, colour=hap1_colsum)) +
  geom_point(alpha = 0.9) +
  labs(x = "hap1", y = "hap2") +
  ggtitle("hap1 vs hap2") +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) + 
  theme_minimal() + theme(legend.position = "none") 


