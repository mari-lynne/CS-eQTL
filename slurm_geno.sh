#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=8000
#SBATCH --output=Rout/par-%j-%a.out   # File to which STDOUT will be written, including job ID (%j)
#SBATCH --error=Rout/par-%j-%a.out    # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

bash geno_prep.sh
