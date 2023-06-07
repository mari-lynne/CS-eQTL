# Cibersort

# Aims:

# Prepare WHI RNAseq data for CibersortX deconvolution
# Plot results

# Read data --------------------------------------------------------------------

source("~/scripts/r-packages.R")
source("~/scripts/functions.R")

wd <- c("~/Documents/CSeQTL/data") # Where main data is saved
plot_dir <- c("~/Documents/CSeQTL/data/plots") # Where to save output plots
id_dir <- ("~/Documents/whi_sca/rna/ids")
meta_dir <- c("~/Documents/whi_sca/rna/meta")
results_dir <- c("~/Documents/CSeQTL/data/results")

data_name <- c("whi_topmed_to6_rnaseq_gene_reads.gct.gz") 
pheno_name <- c("sct_rnaseq_pheno.txt")

setwd(wd)

data <- read_omic(name = data_name, wd = wd)

pheno <- fread(paste(meta_dir, pheno_name, sep = "/"))

pheno <- filter(pheno, flagged_seq_qc != "Yes")
pheno <- pheno %>%
  mutate(plate = coalesce(pheno$pick_plate1, pheno$pick_plate_2_repeat)) %>%
  select(-c("pick_plate1", "pick_plate_2_repeat"))

pheno$plate <- as.factor(pheno$plate)

# Filter data to just include sct genotyped samples 
sub <- data[,colnames(data)%in% pheno$rnaseq_ids]

# Make DGE list object
counts <- as.matrix(sub)

# Genes
gene_info <- data[, 1:2]

# Form DGE list ---------------------------------------------------------------
dge <- DGEList(counts=counts, samples=pheno, genes=gene_info)

# Filter zero counts across all samples
zero_counts <- rowSums(dge$counts==0)==ncol(dge$counts)
dge <- dge[!zero_counts, ,keep.lib.sizes=FALSE]
dim(dge$counts)

# Filter low counts
keep.exprs <- filterByExpr(dge)
dge <- dge[keep.exprs, ,keep.lib.sizes=FALSE]
dim(dge)

# Combat batch correct ---------------------------------------------------------

# Combat (adjust for batch effects) -----------------
exprs <- as.matrix(dge$counts)
pData <- dge$samples
pData <- select(pData, plate, sct) # TODO test with and without sct covar
pData$sct <- as.factor(pData$sct)

combat <- as.data.frame(sva::ComBat_seq(counts=exprs, batch=pData$plate, group=pData$sct, full_mod=TRUE))

row.names(combat) <- dge$genes$Name

dge_bc <- DGEList(counts=combat, samples=pheno, genes=dge$genes)

# TPM --------------------------------------------------------------------------

# get transcript length
# add gene info

dge_bc$genes <- rename(dge_bc$genes, hgnc_symbol = Description)

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

genes <- getBM(
  attributes=c('ensembl_gene_id', 'hgnc_symbol','start_position','end_position'),
  mart = ensembl)

genes$size = genes$end_position - genes$start_position

# 92 duplicates
# TODO update to filter duplicates by lowest expressing transcript
filtered_genes <- 
inner_join(dge_bc$genes %>% group_by(hgnc_symbol) %>% mutate(id = row_number()),
          genes %>% group_by(hgnc_symbol) %>% mutate(id = row_number()), 
          by = c("hgnc_symbol", "id"))

# 2014 transcripts missing ensembl info
# Filter counts table by gene_info
# update DGEList
dge_bc <- dge_bc[row.names(dge_bc$counts) %in% filtered_genes$Name, ]

dge_bc$genes <- filtered_genes

tpm_bc <- as.data.frame(tpm(dge_bc$counts, dge_bc$genes$size))
head(tpm_bc)

tpm_bc$hgnc_symbol <- dge_bc$genes$hgnc_symbol
tpm_bc$dup <- stri_duplicated(tpm_bc$hgnc_symbol)
table(tpm_bc$dup)

tpm_bc <- tidylog::filter(tpm_bc, stri_duplicated(hgnc_symbol)==FALSE) 
tpm_bc <- tidylog::filter(tpm_bc, str_detect(hgnc_symbol, "-")==FALSE)
tpm_bc <- select(tpm_bc,-dup)

tpm_bc <- tpm_bc[, c(806,1:805)]
row.names(tpm_bc) <- tpm_bc$hgnc_symbol
tpm_bc <- tpm_bc %>% remove_rownames

tpm_bc <- select(tpm_bc,-dup)

write.table(tpm_bc, file = paste0(results_dir, "/tpm_rnaseq.txt"), row.names = F, sep = "\t", quote = F)

# Results ----------------------------------------------------------------------

fractions <- clean_names(fread(paste0(results_dir, "/CIBERSORTxGEP_Job1_Fractions-Adjusted.txt")))
colnames(fractions) <- str_c(colnames(fractions), rep("_ciber", ncol(fractions)))

pheno <- clean_names(fread(paste0(meta_dir, "/rnaseq_pheno_all.txt")))
cells <- select(pheno, c("rnaseq_ids", "commonid", "subject_id","age", "bmi_t0","smoking",
                         "neutrophil", "monocytes", "eosinophils", "basophils", "lymphocytes", "rbc_count", "wbc"))

fractions <- rename(fractions, rnaseq_ids = mixture_ciber)
ciber <- left_join(fractions, cells, by = "rnaseq_ids")


flow <- clean_names(read.csv(file = "flow_data__subjectid.csv"))
flow_ciber <- inner_join(flow, ciber, by = "subject_id")
# 72 matches


merged <- clean_names(fread(paste0(results_dir, "/CIBERSORTxGEP_Job3_Weights.txt")))
colnames(merged) <- str_c(colnames(merged), rep("_ciber2", ncol(merged)))
merged <- rename(merged, rnaseq_ids = mixture_ciber2)
flow_ciber <- left_join(flow_ciber, merged, by = "rnaseq_ids")



# Plot correlations

# Monocytes

mono <- 
  ciber %>% ggplot(aes(x=monocytes.y, y =monocytes.x))+
  geom_point() +
  theme_light() +
  labs(x="% cells (clinic)", y = "% cells (ciber)", title = "Monocytes") +
  stat_cor(method = "spearman") +
  geom_smooth(method = lm)

mono

neut <- 
  ciber %>% ggplot(aes(x=neutrophil, y =neutrophils))+
  geom_point() +
  theme_light() +
  labs(x="% cells (clinic)", y = "% cells (ciber)", title = "Neutrophils") +
  stat_cor(method = "spearman") +
  geom_smooth(method = lm)
neut

wbc <-
  ciber %>% ggplot(aes(x=lymphocytes, y =lymphoctes_ciber))+
  geom_point() +
  theme_light() +
  labs(x="% cells (clinic)", y = "% cells (ciber)", title = "Lymphocytes") +
  stat_cor(method = "spearman") +
  geom_smooth(method = lm) + xlim(0,8.5)

(mono|neut|wbc) + plot_annotation(tag_levels = "a")


eo <- 
  ciber %>% ggplot(aes(x=eosinophils.x, y =eosinophils.y))+
  geom_point() +
  theme_light() +
  labs(x="% cells (clinic)", y = "eosinophil % (ciber)") +
  stat_cor(method = "spearman") +
  geom_smooth(method = lm) + xlim(0,0.00001)




# Flow correlations ------------------------------------------------------------
colnames(flow_ciber)

flow_ciber %>% ggplot(aes(x=t_cells_cd8_ciber2, y =cd8_t_cells))+
  geom_point() +
  theme_light() +
  labs(x="% cells (clinic)", y = "% cells (ciber)", title = "Lymphocytes") +
  stat_cor(method = "spearman") +
  geom_smooth(method = lm)

  
  flow_ciber %>% ggplot(aes(x=t_cells, y =lymphoctes_ciber))+
    geom_point() +
    theme_light() +
    labs(x="% cells (clinic)", y = "% cells (ciber)", title = "Lymphocytes") +
    stat_cor(method = "spearman") +
    geom_smooth(method = lm) + xlim(0,8.5) 
  
