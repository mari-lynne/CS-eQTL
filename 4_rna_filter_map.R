#!/usr/bin/env Rscript

## Args and manifest -----------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

# Test if there is not one argument: if not, return an error
if (length(args) != 1) {
  stop("One argument must be supplied, the manifest file name", call.=FALSE)
}

# Read in the specific manifest you indicated, then select the row as arguments.
manifest <- read.csv(args[1] , stringsAsFactors = F) 
thisRun <- manifest[Sys.getenv('SLURM_ARRAY_TASK_ID'),] 

print("mainfest str")
str(manifest)
print("array job structure")
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

###  Directories ---------------------------------------------------------------
# Input Directories
main_dir <- c("/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/")
gen_dir <-  paste0(main_dir, "results/genotype/split/sub/") # From 2_wgs prep
bam_dir <- c("/fh/scratch/delete90/kooperberg_c/sct_rnaseq/release_files/")
  

# Output Directories
rna_dir <- paste0(main_dir, "test_july/rnaseq/")
ase_dir <- paste0(main_dir, "test_july/ASE/")

# Check names
cat(main_dir, gen_dir, bam_dir, rna_dir,ase_dir)

### File names -----------------------------------------------------------------
# RNAseq Files
# Input
bam_file <- paste0(bam_dir, thisRun$bam_file)
# Output
bam_id <- as.character(thisRun$nwgc_sample_id)
bam_filtered <- paste0(rna_dir, bam_id, "-output.filtered.asSeq.bam")
bam_sort <- paste0(rna_dir, bam_id, "-output.filtered.asSeq.sortQ")

# Check names
cat(bam_file, bam_id, bam_filtered, bam_sort)

# WGS Input Files
gtf_rds_fn <- c("exon_by_genes.rds") # Made in 1_get-exon-info
genes <- readRDS(gtf_rds_fn)
geno_id <- as.character(thisRun$topmed_nwdid)
geno_file <- paste0(gen_dir, geno_id, "_hap.txt") # Made in 2_wgs_prep

# Check names
cat(geno_id, geno_file)

# 1) Filter RNA bam files for dups/low qual reads ------------------------------
start_time <- Sys.time()
start <- format(start_time, format = "%H:%M:%S")
print(paste("Starting ASE Pipeline:", start))

# Set up filter parameters
PE <- TRUE

flag1 <- Rsamtools::scanBamFlag(
  isUnmappedQuery = FALSE,
  isSecondaryAlignment = FALSE,
  isDuplicate = FALSE,
  isNotPassingQualityControls = FALSE,
  isSupplementaryAlignment = FALSE,
  isProperPair = PE
)

param1 <- Rsamtools::ScanBamParam(flag = flag1,
                                 what = "seq",
                                 mapqFilter = 20?a)
# Filter bam file
Rsamtools::filterBam(bam_file, destination = bam_filtered, param = param1)


end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("bam filtering DONE!", "Time taken:", time_taken))


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

ct <- as.data.frame(SummarizedExperiment::assay(se)) # TReC

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("Total Read Count DONE!", "Time taken:", time_taken))

## 3) Sort reads by Qname ------------------------------------------------------
start_time <- Sys.time()

sortBam(file = bam_filtered, destination = bam_sort, byQname=TRUE, maxMemory=3500)

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("bam sorting DONE!","Time taken:", time_taken))

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
print(paste("ASE counting DONE!","Time taken:", time_taken))

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

ct1 <- as.data.frame(SummarizedExperiment::assay(se1))
ct2 <- as.data.frame(SummarizedExperiment::assay(se2))
ctN <- as.data.frame(SummarizedExperiment::assay(seN))
cts <- cbind(ct,ct1,ct2,ctN) # trec, hap1, hap2, hapN
dim(cts); cts[1:2,]

out_fn <- paste0(bam_id, "-output.trecase.txt")
write.table(cts, file = paste0(ase_dir, out_fn), quote = FALSE,
            sep = "\t", eol = "\n")

end_time <- Sys.time()
time_taken <- format(end_time - start_time, format = "%H:%M:%S")
print(paste("ASE mapping DONE!","Time taken:", time_taken))

sessionInfo()
q(save="no")
