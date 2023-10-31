#!/usr/bin/Rscript

# UPDATES: 31/07/23
# Running script where files are all in the same location...
# Copied .bam, .bai and _hap.txt files to all_files sub dir # Mapq filter set to 20
# TREC step done before sortQ - still hapN

# Modified hap_n files in wgs_prep_redo.sh, didn't filter 0 alleles - THIS FIXED IT!
# Aug 28: Switched 01 notation to bases - UPDATED dirs
# Oct 9th: Updated input files in wgs_data_formatting, can run on all LLS samples, keeping in 0 alleles

###  Directories ---------------------------------------------------------------

# Inherit bash args
args <- commandArgs(trailingOnly = TRUE)
cat("Inherited arguments:", paste(args, collapse = " "), "\n")
cat("array ", Sys.getenv(' SLURM_ARRAY_TASK_ID'))
cat("array ", args[7])
cat("job ", Sys.getenv('SLURM_JOB_NAME'))

# Input Directories
gen_dir <- args[2]
bam_dir <- args[3]
ref_dir <- args[4]

# Output Directories
rna_dir <- args[5]
ase_dir <- args[6]

# Check names
cat(gen_dir, bam_dir, ref_dir, rna_dir, ase_dir)
 
# # Read in the specific manifest you indicated, then select the row as arguments.
manifest <- read.csv(args[1], stringsAsFactors = F)
str(manifest)
thisRun <- manifest[Sys.getenv('SLURM_ARRAY_TASK_ID'),]
str(thisRun)

### Packages ------------------------------------------------------------------
# Installed extra packages to rhino home dir https://sciwiki.fredhutch.org/rModules
library(stringr)
library(data.table)
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(CSeQTL)
library(asSeq)

### File names -----------------------------------------------------------------

# IDs
bam_id <- as.character(thisRun$lls_torid)
geno_id <- as.character(thisRun$topmed_nwdid)

# RNAseq Input Files
bam_file <- paste0(bam_dir, bam_id, ".Aligned.sortedByCoord.out.md.bam")

# Output
bam_filtered <- paste0(rna_dir, bam_id, "-output.filtered.asSeq.bam")
bam_sort <- paste0(rna_dir, bam_id, "-output.filtered.asSeq.sortQ")

# WGS Input Files
gtf_rds_fn <- paste0(ref_dir, "exon_by_genes.rds") # Made in 1_get-exon-info
genes <- readRDS(gtf_rds_fn)
geno_file <- paste0(gen_dir, geno_id, "_hap.txt") # Made in wgs_prep

# Output
ase_out_file <- paste0(ase_dir, geno_id, "-output.trecase.txt")

# Check names
cat("rna files", bam_file, bam_id, bam_filtered, bam_sort)
cat("wgs files", geno_id, geno_file, ase_out_file)

# 1) Filter RNA bam files for dups/low qual reads ------------------------------
start_time <- Sys.time()
start <- format(start_time, format = "%H:%M:%S")
print(paste("Starting ASE Pipeline:", start))

## Set up filter parameters
PE <- TRUE

flag1 <- Rsamtools::scanBamFlag(
isUnmappedQuery = FALSE,
isSecondaryAlignment = FALSE,
isDuplicate = FALSE,
isNotPassingQualityControls = FALSE,
isSupplementaryAlignment = FALSE,
isProperPair = PE
)

param1 <- Rsamtools::ScanBamParam(flag = flag1, what = "seq", mapqFilter = 20)

##  Filter bam file
Rsamtools::filterBam(bam_file,
		     destination = bam_filtered,
		     param = param1,
		     maxMemory=7500,
		     nThreads=2)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("bam filtering - DONE!", "Time taken:", time_taken))


# 2) Get total read count per gene ---------------------------------------------
start_time <- Sys.time()

bamfile <- Rsamtools::BamFileList(bam_filtered, yieldSize = 1000000)
se <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = bamfile,
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)

TReC <- as.data.frame(SummarizedExperiment::assay(se)) 
write.table(TReC, file = paste0(ase_dir, geno_id, "-TReC.txt"), quote = FALSE, sep = "\t")

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("Total Read Count - DONE!", "Time taken:", time_taken))

## 3) Sort reads by Qname ------------------------------------------------------
start_time <- Sys.time()

sortBam(file = bam_filtered,
        destination = bam_sort,
        byQname=TRUE,
        maxMemory=7500,
	nThreads=2)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("bam sorting - DONE!","Time taken:", time_taken))

# 4) Get allele-specific read counts -------------------------------------------
start_time <- Sys.time()

asSeq::extractAsReads(
  input = paste0(bam_sort, ".bam"),
  snpList = geno_file,
  min.avgQ = 20,
  min.snpQ = 20
)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("ASE counting - DONE!","Time taken:", time_taken))

## 4b) Count allele specific read counts ASRC ----------------------------------
start_time <- Sys.time()

se1 <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = sprintf("%s_hap1.bam", bam_sort),
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)
se2 <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = sprintf("%s_hap2.bam", bam_sort),
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)
seN <- GenomicAlignments::summarizeOverlaps(
  features = genes,
  reads = sprintf("%s_hapN.bam", bam_sort),
  mode = "Union",
  singleEnd = !PE,
  ignore.strand = TRUE,
  fragments = PE
)

hap1 <- as.data.frame(SummarizedExperiment::assay(se1))
hap2 <- as.data.frame(SummarizedExperiment::assay(se2))
hapN <- as.data.frame(SummarizedExperiment::assay(seN))
cts <- cbind(TReC,hap1,hap2,hapN)
dim(cts); cts[1:2,]

write.table(cts, file = ase_out_file, quote = FALSE, sep = "\t", eol = "\n")

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("ASE mapping - DONE!","Time taken:", time_taken))
