### Date: 10-10-23
### Aims: Plot LLS PCA data, also plot LLS data projected on 1KG PCs


# Set-up -----------------------------------------------------------------------

# Load packages
library(dplyr)
library(tidylog)
library(patchwork)
library(ggplot2)
library(colorspace)
library(readr)
library(data.table)
library(forcats)
library(viridis)

## Initialise command line arguments
args <- commandArgs(TRUE)
root <- args[1]

# Directories
in_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap")
in_dir2 <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap/qc_ancestry")
plot_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/plots")

# Input data 
# Study PCA data
pca <- fread(file.path(in_dir, "whi_lls.LD.eigenvec"))
eigenval <- fread(file.path(in_dir, "whi_lls.LD.eigenval")) # For PVE plot
covar <- read.csv(file.path(in_dir, file = "LLS_ID_pheno_08-23.csv"))

# Ancestry merged PCA data 
PCA <- read_delim(file.path(in_dir2,"1KG.merged.eigenvec"), 
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# NO EAS data
PCA <- read_delim(file.path(in_dir2,"1KG.noEAS.eigenvec"), 
                  delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# 1KG population phenotypes
pop <- read_delim(file.path(in_dir2, "1kG_ID2Pop.txt"), delim = " ", escape_double = FALSE, trim_ws = TRUE)
pop <- pop %>% rename(IID = `#IID`)

# Just study WGS data ---------------------------------------------------------

### Format Data ----------------------------------------------------------------

covar$IID <- covar$topmed_nwdid
pca$IID <- pca$`#IID`
pca_covar <- left_join(pca, covar, by = c("IID"))

# Recode ethnicity:
pca_covar$ethnicity <- as.factor(recode(pca_covar$ethnic.x, '3' = 'A.A', '4' = 'His','5' = 'E.A'))

# PVE
pve <- data.frame(PC = 1:10, pve = (eigenval/ sum(eigenval) * 100))
pve$PC <- as.factor(pve$PC)
cumsum(pve$V1)

### Palettes -------------------------------------------------------------------

PC_Palette <- heat_hcl(10, h = c(30, -160), c = c(80, NA, 45), l = c(38, 79),
    power = c(0.85, 1.0), fixup = TRUE, gamma = NULL, alpha = 1)

n_colors <- length(levels(pve$PC))

pve_pal <- viridis(
  n = n_colors, 
  option = "plasma",    
  begin = 0.1,          
  end = 0.9,            
  direction = -1,
  alpha = 1
)

# Pal to match 1KG colours # "#FCCE25D9" = Eur, #CC4678D9 = AFR, #900DA4D9 = His
levels(pca_covar$ethnicity)
pal <- c("#CC4678D9", "#FCCE25D9", "#900DA4D9")


### Plots -----------------------------------------------------------------------

# Plot proportion of variance
pve_plot <- ggplot(pve, aes(PC, V1, colour = PC)) +
  geom_line(group = "PC", colour = "#454b51") +
  geom_point(size = 2.5)  +
  ylab("\nProportion of Variance (%)\n") +
  scale_colour_manual(values = pve_pal) +
  theme_bw() + theme(legend.position = "none") 

pve_plot

ggsave("pve_study.png", path = plot_dir)

pca_plot <- pca_covar %>%
  ggplot(aes(PC1, PC2, colour = ethnicity)) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values=pal)+ theme_bw() + labs (color="Ethnicity")
pca_plot


# PCA and PVE combined
pdf(file.path(plot_dir, "LLS_PCA_PVE.pdf"))
(pca_plot/pve_plot) + plot_annotation(title = "LLS WGS PCA")
dev.off()

ggsave(filename = file.path(plot_dir, "LLS_PCA_PVE.png"))


# 1 KG Ancestry ----------------------------------------------------------------

### Format Data ----------------------------------------------------------------

# Rename columns
colnames(PCA) <- c("FID","IID", "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

# Merge with 1KG Population data
merge <- left_join(PCA, pop, by = "IID")

# Rename NAs in Pop data as Study
study_name <- "LLS Study"
merge$SuperPop <- merge$SuperPop %>% replace_na(study_name) 
merge$Population <- merge$Population %>% replace_na(study_name)

# Recode to factors
merge$SuperPop <- as.factor(merge$SuperPop)
levels(merge$SuperPop)

merge$Population <- as.factor(merge$Population)
levels(merge$Population)

# Order labels for plotting
merge$SuperPop <- fct_relevel(merge$SuperPop, "EUR", study_name)
merge$Population <- fct_relevel(merge$Population, study_name)
merge <- merge[order(merge$SuperPop), ]
merge <- merge[order(merge$Population), ]


### Palettes -------------------------------------------------------------------

KG_Palette<-heat_hcl(length(unique(merge$Population)),
                     h = c(300, 75), c. = c(35, 95),
                     l = c(15, 90), power = c(0.8, 1.2),
                     fixup = TRUE, gamma = NULL, alpha = 1)

n_colors <- length(levels(merge$SuperPop))
pal1 <- viridis(
  n = n_colors, 
  option = "plasma",    
  begin = 0.1,          
  end = 0.9,            
  direction = -1,
  alpha = 0.85
)

pal1
pal1 <- c("#F1844BD9","#FCCE25D9","#CC4678D9", "#900DA4D9", "#42049ED9") # reorder colours

n_colors2 <- length(levels(merge$Population))
pal2 <- viridis(
  n = n_colors2, 
  option = "plasma",    
  begin = 0.0,          
  end = 0.9,            
  direction = -1,
  alpha = 0.9
)

### Plots ----------------------------------------------------------------------

a <- ggplot(merge, aes(PC1, PC2, colour = SuperPop)) +
  geom_point(alpha = 0.7) +
  scale_colour_manual(values = (pal1)) + theme_bw() + theme(legend.position = "right") + 
  labs(color = "Population (1KG)")
a 

# , title = paste0("Population Stratification")

pdf(file = file.path(plot_dir, "1KG_superpop.pdf"))
a
dev.off()
ggsave(a, filename = "1KG_superpop.png", path = plot_dir)


b <- ggplot(merge, aes(PC1, PC2, colour = Population)) +
  geom_point() +
  scale_colour_manual(values = (pal2)) +
  theme_bw() + theme(legend.position = "bottom") + 
  labs(color = "")

b

pdf(file = file.path(plot_dir,"1KG_all_pop_pca.pdf"))
b
dev.off()
ggsave(b, filename = "1KG_all_pop.png", path = plot_dir)

# all 1KG plots

pdf(file.path(plot_dir, "1KG_all_pcas.pdf"))
(a|b)
dev.off()
ggsave(filename = file.path(plot_dir, "1KG_all_pcas.png"))

# PCA and PVE combined

pdf(file.path(plot_dir, "1KG_LLS_PCA_PVE.pdf"))
(a/pve_plot)
dev.off()

ggsave(filename = file.path(plot_dir, "1KG_PCA_PVE.png"))

# "#FCCE25D9" = Eur, #CC4678D9 = AFR, #900DA4D9 = His


# EAS check/removal -----------------------------------------------------------

merge %>% filter(SuperPop == "EAS"| SuperPop =="LLS") %>% ggplot(aes(PC1, PC2, colour = SuperPop)) +
  geom_point() +
  theme_bw() +
  ggtitle(paste0(study_name, " Population Stratification (1KG)"))
ggsave(filename = "LLS_vs_EAS.png")

# no participants with EAS genetic ancestry, some near SAS, exclude 

# get list of participants to filter
eas <- filter(merge, SuperPop == "EAS") %>% select(IID)
write.table(eas, file = "EAS_rm.txt", quote = F, sep = "\t", row.names = F, col.names = F)

#./plink2 --bfile qc_ancestry/1KG.merged --remove qc_ancestry/EAS_rm.txt --make-bed --out qc_ancestry/1KG_noEAS
#./plink2 --bfile qc_ancestry/1KG_noEAS --pca --out qc_ancestry/1KG.noEAS

# Use noEAS data for final pop pca plots
