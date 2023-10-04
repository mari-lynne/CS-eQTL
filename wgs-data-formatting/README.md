# CSeQTL WGS Formatting Workflow

## Overview

This repository contains scripts for formatting eQTL data and performing related analyses. The scripts are designed to run on a SLURM cluster. They perform tasks such as transforming genotype data, preparing data for ASeq and CSeQTL analyses, and managing job dependencies. Ultimately we need two file types for aseq and cseqtl analysis from our bcf files:  

#### 1) ASeq (per sample)
| CHR  | POS     | HAP1 | HAP2 |
|------|---------|------|------|
| chr21| 5131973 | C    | C    |
| chr21| 5131987 | G    | G    |

#### 2) CSeQTL (per SNP)
| SAMPLE_ID | GENO |
|-----------|------|
| NWD446825 | 0    |
| NWD446777 | 2    |

---

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [How to Run](#how-to-run)
3. [Scripts](#scripts)
    - [2_eqtl_format_submit.sh](#2_eqtl_format_submitsh)
    - [2_eqtl-format-pt1.sh](#2_eqtl-format-pt1sh)
    - [2_eqtl-format-cseqtl-pt2.sh](#2_eqtl-format-cseqtl-pt2sh)
    - [2_eqtl-format-aseq-pt3.sh](#2_eqtl-format-aseq-pt3sh)

---

## Prerequisites

- Linux/Unix operating system
- SLURM job scheduler
- BCFtools module
- Input directory should contain array file `lls_sample_array_input.csv`
  
---

## How to Run

1. Navigate to the directory containing the scripts.
2. Ensure submission script is executable (chmod +x submit.sh)
3. Run `2_eqtl_format_submit.sh` to initiate the workflow.

---

## Scripts

### 2_eqtl_format_submit.sh

#### Description
The main script to submit jobs for eQTL data formatting. This script loops through specified chromosomes and submits parallel jobs to perform formatting tasks on bcf files for aseq/cseqtl analysis.

#### Variables
- `IN_DIR`: Input directory path.
- `SCRIPT_DIR`: Directory where the other scripts are stored.
- `chromosomes`: Array of chromosome numbers to be processed.

#### Usage
Run this script to initiate the entire workflow.

---

### 2_eqtl-format-pt1.sh

#### Description
This script takes bcf files that have been split by chromosome and creates two intermediate files:
1) ASeq counting (sample_hap_chr.txt)
2) CSeQTL genotype files (chr_temp_snp_data.txt).

#### Variables
- `chr`: Chromosome number passed from the main script.
- `GENO_FILE`: The base name of the genotype file.
- `IN_DIR`: Input directory.
- `OUT_DIR`: Output directory for haplotype files.

#### Usage
This script is called by `2_eqtl_format_submit.sh`.

#### Example output: (dummy data shown)  

ASeq: (per sample)  
chr21  5131973  C    C  
chr21  5131987  G    G  
chr21  5131989  A    A  


CSeQTL: (all samples per chromosome)      
NWD446825 chr21 5131973 0     
NWD446777 chr21 5131973 1

---

### 2_eqtl-format-cseqtl-pt2.sh

#### Description
This script creates per-SNP files for CSeQTL analysis. It aggregates genotype data stored in chr_temp_snp_data.txt by SNP site.

#### Variables
- `chr`: Chromosome number.
- `input_file`: Temporary SNP data file to read.

#### Usage
This script is called by `2_eqtl_format_submit.sh` and is dependent on the successful completion of `2_eqtl-format-pt1.sh`.

#### Output

CSeQTL: (per SNP file e.g chr21_5031973.txt)
SAMPLE_ID GENO
NWD446825 0
NWD446777 2

---

### 2_eqtl-format-aseq-pt3.sh

#### Description
This script concatenates haplotype files which have been split per chromosome/individual, into files grouped per individual (with all chromosome data).

#### Variables
- `SAMPLE_NAME`: Identifier for the sample.
- `IN_DIR`: Input directory for haplotype files.
- `OUT_DIR`: Output directory for concatenated haplotype files.

#### Usage
This script is called by `2_eqtl_format_submit.sh` and is dependent on the successful completion of `2_eqtl-format-pt1.sh`.

#### Output
Same as previous ASeq output, except will be concatenated by chromsome

---
### Contact

For any questions or issues, please contact mjohnso5@fredhutch.org.
