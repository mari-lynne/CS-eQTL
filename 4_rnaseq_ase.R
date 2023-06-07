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

