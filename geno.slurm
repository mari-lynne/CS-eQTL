#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --output=Rout/par-%j-%a.out#
#SBATCH --error=Rout/par-%j-%a.out
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

ml BCFtools/1.14-GCC-11.2.0

bash geno_prep.sh
