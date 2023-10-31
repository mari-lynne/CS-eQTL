library(data.table)
library(dplyr)
library(stringr)
library(viridis)
library(ggplot2)
library(edgeR)
library(limma)
library(patchwork)
library(PCAForQTL)
source("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/scripts/functions.R")

# Dirs and vars ----------------------------------------------------------------

# Input Directories
in_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq"
meta_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/metadata"
# Output Directories
out_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq"
plot_dir <- "/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/rnaseq/plots"
setwd(out_dir)

# Input files
data_name <- c("TOPMed_Kooperberg_P4.RNASeQCv2.3.4_gene_reads.gct.gz") 
covar_name <- c("LLS_ID_pheno_08-23.csv")

# Functions:
# Calculate tpm:
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

# 1) Load data -----------------------------------------------------------------
# RNA seq
data <- read_omic(name = data_name, dir = in_dir)
covar <- read_omic(name = covar_name, dir = meta_dir)

# Filter samples
sub <- data[ , colnames(data) %in% covar$lls_torid]

# Gene_info
gene_info <- data[, 1:2]
gene_info$ensembl_gene_id <- str_remove_all(gene_info$Name, "\\.\\d")

# 2) Preprocess ----------------------------------------------------------------

counts <- as.matrix(sub)
dge <- DGEList(counts=counts, samples=covar, genes=gene_info)
dim(dge) # 58103  1327

# Remove zero counts
dge  <- dge[rowSums(dge$counts) > 0, ]

# Keep samples that have above 10 counts, in at least n samples
# n samples determined by the smallest group sample size, default = 70%, so in our case would be 
# above 10 counts in at least 70% of 96 (hispanics) 67 samples, or do of whole data (leave out), nums match the sct more that way
# Plot before and after filtering low counts
keep.exprs <- filterByExpr(dge)
dge2 <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge2)

# Calculate effective library sizes by normalising with TMM
dge2 <- calcNormFactors(dge2, method = "TMM")
head(dge2$samples$norm.factors)

# Normalize counts data using edgeR cpm function, log space minimizes the effects of outliers
cpms <- cpm(dge2, log=TRUE) # USE CPM unscalled for PCA


# 3) Outlier detection ---------------------------------------------------------

# Two approaches: Relative Log Expression analysis, and clustering
# For outlier removal we will calculate TPM values and work in log10 space

## Calculate TPM ----------------------------------------------------------------

# Get more gene_info
genes <- read.csv(file.path("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/resources", "gene_length_info.csv"))
# filter duplicates
genes <- genes %>% distinct(ensembl_gene_id, .keep_all = TRUE)
# Just keep genes with gene length info from biomart
filtered_genes <- inner_join(dge2$genes, genes, by = "ensembl_gene_id")
# 5357 matched rows
# Filter dge2 for matched rows then convert to TPM 
dge2 <- dge2[dge2$genes$ensembl_gene_id %in% filtered_genes$ensembl_gene_id, ]

# Extract count data
counts <- dge2$counts
row.names(counts) <- dge2$genes$ensembl_gene_id

# Convert matrix to Transcripts Per Million (TPM)
counts_tpm <- as.data.frame(tpm(counts, filtered_genes$size))

# Log10 plus pseudo count
pseudo_count = 0.01
logtpm <- log10(counts_tpm + pseudo_count)

### RLE ------------------------------------------------------------------------

# Filter outliers using Relative log expression:
RLEFilterPercent = 0.01
RLEFilterLength <- RLEFilterPercent*ncol(counts_tpm)

rle <- logtpm-apply(logtpm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
iqr <-  apply(rle,2,IQR)
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample", value.name ="TPM", id="ID")
names(rle) <- c("feature","Sample","TPM")
rle_IQR <- rle %>% group_by(Sample) %>% summarise(IQR = IQR(TPM))
rle_IQR_range <- rle_IQR$IQR %>% range %>% abs() %>% max()
rle_IQR_range <- 2*rle_IQR_range %>% ceiling()
bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
par(mar=c(3,3,3,3))
boxplot(TPM ~ bymedian, data=rle, outline=F, ylim = c(-rle_IQR_range, rle_IQR_range), las=2, boxwex=1, col='gray', cex.axis=0.3, main="RLE plot before QC", xlab="", ylab="Residual expression levels", frame=F)
ExpPerSample <- nrow(counts_tpm)
RLEFilterList <- as.character(unique(bymedian[((length(bymedian)-ExpPerSample*RLEFilterLength)+1):length(bymedian)])) #filtered
print(paste0("The right most ", RLEFilterPercent*100, "% samples (N = ", length(RLEFilterList), ") are marked as candidate outliers in this step:") )
RLEFilterList

# check covars of filter list
check <- covar[covar$lls_torid %in% RLEFilterList, ]

### Heirachecal Clustering -------------------------------------------------------

# hcluster thresholds
topk_genes = 100
cluster_percent = 0.6
pvalues.cut = 0.05
treesNum = 5 # Cluster Level
sampleDists <- 1 - cor(logtpm, method='spearman')
hc <- hclust(as.dist(sampleDists), method = "complete")
hcphy <- as.phylo(hc)

plot(hcphy, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Sample clustering before QC (Spearman - Cor.)")

ntop <- topk_genes
Pvars <- apply(logtpm, 1, var)
select <- order(Pvars, decreasing =TRUE)[seq_len(min(ntop, length(Pvars)))]
MD_matrix <- logtpm[select, ]
MahalanobisDistance = mahalanobis(t(MD_matrix), colMeans(t(MD_matrix)), cov(t(MD_matrix))) 
# Note: t(MD_matrix)=sample_row*gene_column, Manhalanobis() returns one vector with length=row number
pvalues = pchisq(MahalanobisDistance, df=nrow(MD_matrix), lower.tail=F)
pvalues.adjust = p.adjust(pvalues, method ="bonferroni") # adjusted pvalues for each sample
pvalues.low <- pvalues.adjust[pvalues.adjust<pvalues.cut]

HCoutliers <- character()
for(x in c(1:treesNum)){
  trees <- cutree(hc,k=x)
  idx <- c(1:x)#which tree is checking
  for(i in idx)
  {
    group <- hc$labels[which(trees == i)]
    if(sum(group %in% names(pvalues.low))/length(group) >= cluster_percent)
    {
      HCoutliers <- union(HCoutliers, group)
    }
  }
}

print(paste(length(HCoutliers), "samples are marked as candidate outlier(s) in this step.", sep = " "))
if(length(HCoutliers)>0){
  print("Sample outliers are marked in red as follows:")
  print(HCoutliers)
  co1 = hc$labels%in%HCoutliers
  co1[which(co1 == "FALSE")]<-"gray0"
  co1[which(co1 == "TRUE")]<-"red"
  par(mar=c(3,3,3,3))
  plot(hcphy, tip.col = co1, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Hierarchical clustering: Outlier Detection")
  Xcol = c("gray0", "red")
  Xtext = c("Normal Sample", "Outliers")
  legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=1)
}else{
  print("No outlier detected.")
}

### Filter outliers -----------------------------------------------------------
# Remake DGEList object with log10 TPM data and filter outliers

dge_tpm <- DGEList(counts = counts_tpm, samples = dge2$samples, genes = dge2$genes)
dge_tpm <- dge_tpm[ ,dge_tpm$samples$lls_torid %!in% HCoutliers] # 1325 samples 

dge2 <- dge2[ ,dge2$samples$lls_torid %!in% HCoutliers]
cpms <- cpm(dge2, log=TRUE)

# 4) PCA -------------------------------------------------------------

# Should use normalised log transformed data, either can be thru edgeR method (cpm)
# Or DSeq2 vst - compare both

pca_result <- prcomp(t(cpms), center=TRUE, scale=FALSE)

plot(pca_result)
biplot(pca_result)
pca_df <- as.data.frame(pca_result$x) %>%
  mutate(ethnic = as.factor(dge2$samples$ethnic.y))

# Update factors
table(pca_df$ethnic)
pca_df$ethnic <- recode(pca_df$ethnic,
                        '3' = 'A.A',
                        '4' = 'His',
                        '5' = 'E.A')

pal <- c("#CC4678D9", "#900DA4D9", "#FCCE25D9")

pca_grid_plot <- create_pca_grid(pca_df)
print(pca_grid_plot)

pc1 <- pca_df %>%
  ggplot(aes(PC1, PC2, color = ethnic)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = pal) +
  theme_bw() +
  labs(color = "Ethnicity") +
  ggtitle("LLS RNA-Seq PCA")

pc1


#### 4.1) PVE ---------------------------------------------------------------

# Step 1: Calculate PVE from the prcomp object
pve_values <- (pca_result$sdev^2) / sum(pca_result$sdev^2)

# Step 2: Create a data frame for PVE values
pve_df <- data.frame(
  PC = factor(seq_along(pve_values)),
  V1 = pve_values * 100  # Multiply by 100 to represent as a percentage
)

num_pcs = 25
pve_df <- pve_df[1:num_pcs,]

# PVE Plot
n_colors = num_pcs
pve_pal <- viridis(
  n = n_colors, 
  option = "plasma",    
  begin = 0.1,          
  end = 0.9,            
  direction = -1,
  alpha = 1
)

# Plot proportion of variance
pve_plot <- ggplot(pve_df, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = pve_pal) +
  theme_bw() + theme(legend.position = "none") 

pve_plot
(pc1/pve_plot)

# Check Elbow:

PCs <-pca_result$x
calc_elbow <-PCAForQTL::runElbow(prcompResult=pca_result)
print(calc_elbow)

# Scree Plots 

PCAForQTL::makeScreePlot(pca_result,
                         labels=c("Elbow"),
                         values=c(calc_elbow),
                         maxNumOfPCsToPlot = 40,
                         titleText="LLS RNA-Seq PCA")


