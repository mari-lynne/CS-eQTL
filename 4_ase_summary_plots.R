# Plot summary histograms ------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(patchwork)

# List all the files in your directory with the pattern "-output.trecase.txt"
# setwd("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/ASE")

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

hap1 <-   ggplot(summary_data, aes(x = trec_colsum, y = hap1_colsum)) +
  geom_point() +
  labs(x = "trec", y = "hap1") +
  ggtitle("Scatter Plot: hap1 vs TReC") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) + 
  theme_minimal()

hap2 <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hap2_colsum)) +
  geom_point() +
  labs(x = "trec", y = "hap2") +
  ggtitle("hap2 vs TReC") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) + 
  theme_minimal()

hapN <- 
  ggplot(summary_data, aes(x = trec_colsum, y = hapN_colsum)) +
  geom_point() +
  labs(x = "trec", y = "hapN") +
  ggtitle("hapN vs TReC") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) + 
  theme_minimal()

(hap1 | hap2 |hapN) + plot_annotation(tag_levels = "a")

# HAP1 v HAP2 
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
  


