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

##  Input output files -----------------------------------------------------------

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

gtf_rds_fn = c("exon_by_genes.rds") # Made in 1_get-exon-info
genes = readRDS(gtf_rds_fn)
geno_file = paste0(gen_dir, thisRun$topmed_nwdid, "_hap.txt") 


# 1) Filter RNA bam files for dups and low qual reads ---------------------------------

PE = TRUE

flag1 = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                               isSecondaryAlignment = FALSE, isDuplicate = FALSE,
                               isNotPassingQualityControls = FALSE,
                               isSupplementaryAlignment = FALSE, isProperPair = PE)

param1 = Rsamtools::ScanBamParam(flag = flag1, what = "seq", mapqFilter = 255)

# Flag an integer(2) vector used to filter reads based on their 'flag' entry. This is most easily created with the scanBamFlag() helper function.
# mapqFilter = minimum mapping quality to include. BAM records with mapping qualities less than mapqFilter are discarded.
# param1 is therefore a combination of flag filters and mapqFilter

Rsamtools::filterBam(file = bam_file, destination = bam_filt_fn, param = param1)

print("bam filtering DONE")

# 2) Get total read count per gene ---------------------------------------------

bamfile = Rsamtools::BamFileList(bam_filt_fn, yieldSize = 1000000)

se = GenomicAlignments::summarizeOverlaps(features = genes,     
                                          reads = bamfile, mode = "Union",singleEnd = !PE,
                                          ignore.strand = TRUE, fragments = PE)

ct = as.data.frame(SummarizedExperiment::assay(se)) # reads

# Filter reads by Qname
sort_bam = paste0(rna_dir, bam_id, "-output.filtered.asSeq.sortQ")
sortBam(file = bam_filt_fn, destination = sort_bam, byQname=TRUE)

print("Total Read Count DONE")

