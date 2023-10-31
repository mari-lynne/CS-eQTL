# 2/10/23
# Aims
# Make gene coordinate file for splitting gene/snp matricies

# 1) Get gene coords ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(biomaRt)
library(dplyr)

# set up mart
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# get coordinates
genes <- getBM(
  attributes=c('chromosome_name','start_position','end_position','ensembl_gene_id', 'hgnc_symbol'),
  mart = ensembl)
# filter for autsomes and for now only those annotated by hgnc
auto <- genes[genes$chromosome_name %in% c(1:22) & genes$hgnc_symbol != "", ]

# 2) Convert to eQTL coordinates (+-100KB) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
auto$start_position <- auto$start_position - 100000
auto$end_position <- auto$end_position + 100000
# For now replace neg start positions (TODO check boundries)
auto$start_position[auto$start_position < 0] <- 1
# Order DF by chromosome (1:22) and then start position
auto$chromosome_name <- as.numeric(auto$chromosome_name)
auto <- auto[order(auto$chromosome_name, auto$start_position), ]
auto <- select(auto, ensembl_gene_id, chromosome_name, start_position, end_position)

# Split into chromosome files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
chromosomes <- split(auto, auto$chromosome_name) # split into lists

for (chr in names(chromosomes)) {
  chr_df <- chromosomes[[chr]]
  file_name <- paste("~/Documents/CSeQTL/autosome_eqtl_", chr, ".txt", sep="")
  write.table(chr_df, file=file_name, col.names=F, quote=F, sep="\t", row.names=F)
}

# Save file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
write.table(auto, file = "~/Documents/CSeQTL/autosome_eqtl_coord.txt",
            col.names = F, quote = F, sep = "\t", row.names = F)
write.table(auto[5:25,], file = "~/Documents/CSeQTL/test_coord.txt",
            col.names = F, quote = F, sep = "\t", row.names = F)

