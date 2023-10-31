#!/usr/bin/env Rscript

# In terminal of wdir
# wget ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v43.annotation.gtf.gz
# Run script as ml fhR/4.2.2.1-foss-2021b - Rscript /home/mjohnso5/CSeQTL/scripts/1_get-exon-info.R

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(asSeq)
library(CSeQTL)

# Input
gene_info_dir = c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl")
gtf_fn = file.path(gene_info_dir,"gencode.v43.annotation.gtf.gz")
# Output
gtf_rds_fn =  ("exon_by_genes_gencode_v43.rds")

# Prep gene info ---------------------------------------------------------------

# exons by genes
prep_gene_info(gene_info_dir, gtf_fn) # Saves rds object exon_by_genes.rds and gene_info.tsv.gz



