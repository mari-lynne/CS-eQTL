#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=ancestry
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=14
#SBATCH --output=ancestry_%j.out
#SBATCH --error=ancestry_%j.err
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=mjohnso5@fredhutch.org

# Set up directories -------------------------------------------------------

file=whi_lls #name of study PLINK files
refname=all_hg38

studydir=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/genotype/LLS/Sept/LLS_Hap
refdir=/fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/resources

mkdir -p $studydir/qc_ancestry 
qcdir=$studydir/qc_ancestry # qcdir will contain the cleaned study and refernce data
log=$studydir/log

# Download reference data -------------------------------------------------------

cd $refdir

pgen=https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/vx09262b4k1kszy/all_hg38.pvar.zst?dl=1
sample=https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1
king=https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id?dl=1

# wget $pgen
mv 'all_hg38.pgen.zst?dl=1' all_hg38.pgen.zst
$studydir/plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen

wget $pvar
mv 'all_hg38.pvar.zst?dl=1' all_hg38.pvar.zst
wget $sample
mv 'hg38_corrected.psam?dl=1' all_hg38.psam
wget $king
mv 'deg2_hg38.king.cutoff.out.id?dl=1' deg2king.id


echo "1000G data download- DONE"

# Convert ref and study data to Bed
cd ${studydir}

./plink2 \
--pfile $refdir/$refname vzs --max-alleles 2 \
--allow-extra-chr \
--autosome \
--remove $refdir/deg2king.id \
--make-bed \
--out $refdir/$refname

mv $refdir/$refname.log $log

./plink2 \
--pfile $studydir/$file.LD --max-alleles 2 \
--allow-extra-chr \
--autosome \
--make-bed \
--out $qcdir/$file.LD

mv *.log $log

echo "data conversion - DONE!"

# 1) Tidy Study and Reference data---------------------------------------------------

# 1a) Prune Study Data --------------------------------------------------
# Study data already pruned

./plink2 \
--bfile $refdir/$refname \
--indep-pairwise 50 5 0.2 \
--out $refdir/$refname.LD

./plink2 \
--bfile $refdir/$refname \
--allow-extra-chr \
--extract $refdir/$refname.LD.prune.in \
--make-bed \
--out $refdir/$refname.LD

echo "Reference pruning - DONE!"

# 2) Remove AC-GT SNPs -----------------------------------------------------
#(these snps are difficult to merge)

cd ${refdir}
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
$refdir/$refname.LD.bim  > \
$refdir/$refname.acgt #save ref_genome without ac/gt snps list to qc directory

cd ${qcdir}
awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
$qcdir/$file.LD.bim  > \
$qcdir/$file.acgt

echo "AC-GT SNP list done"

cd ${studydir}
# Remove acgt snps
./plink2 \
--bfile $refdir/$refname.LD \
--allow-extra-chr \
--exclude $refdir/$refname.acgt \
--make-bed \
--out $refdir/$refname.no_acgt # clean ref data

./plink2 \
--bfile $qcdir/$file.LD \
--exclude $qcdir/$file.acgt \
--make-bed \
--out $qcdir/$file.no_acgt  # clean study data


# 2b) Remove duplicated snps and keep just atcg snps -------------------------------------------------

./plink2 \
--bfile $refdir/$refname.no_acgt \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-bed \
--out $qcdir/$refname.cleaned

./plink2 \
--bfile $qcdir/$file.no_acgt \
--snps-only just-acgt \
--rm-dup exclude-all \
--make-bed \
--out $qcdir/$file.cleaned

echo "ATCG SNPs removed"

3) Match annotations! --------------------------------------------------------

Update ref chromsome annotation in bim file

awk '{if ($1 != 0) print $2,"chr"$2}' $qcdir/$refname.cleaned.bim > $qcdir/updateFormat.txt

./plink2 \
--bfile $qcdir/$refname.cleaned \
--update-name $qcdir/updateFormat.txt \
--make-bed \
--out $qcdir/$refname.cleaned

echo "Update ref chr format - DONE"

Update study rsids to chr:bp:alt:ref format (as it looks like in 1KG)
./plink2 \
--bfile $qcdir/$file.cleaned \
--set-all-var-ids chr@:#:\$r:\$a \
--make-bed \
--out $qcdir/$file.cleaned

4) Filter reference data for study SNPs --------------------------------------

awk '{print $2}' $qcdir/$file.cleaned.bim > $qcdir/study_keep_list.txt

./plink2 \
--bfile $qcdir/$refname.cleaned \
--extract $qcdir/study_keep_list.txt \
--make-bed \
--out $qcdir/$refname.forMerge

echo "Keep List - Done!"

# 5) Test merge ----------------------------------------------------------------

# Merge study data with reference.forMerge data - will fail first time
# This gives us a list of snps which we can exclude
# You could also flip the snps and try remerging, but that's not necessary here

./plink \
--bfile $qcdir/$file.cleaned \
--bmerge $qcdir/$refname.forMerge --merge-mode 6 \
--out $qcdir/$refname.merge_failures

# Exclude SNPs that fail merge -------------------------------------------------

# extract SNPs that failed to merge from log file
# awk prints fields containing SNP IDs in log file
# Sed removes single quotes from output

grep "chr" $qcdir/$refname.merge_failures.log |\
awk 'BEGIN {OFS="\t"} {
if ($2 == "Multiple")
	print $7;
else
	print $3;
}'| sed -e "s/'//g" > $qcdir/exclude_list.txt

echo "Exclude List - Done!"

# Exclude SNPs
./plink2 \
--bfile $qcdir/$refname.forMerge \
--exclude $qcdir/exclude_list.txt \
--make-bed \
--out $qcdir/$refname.cleanMerge

# 6) Remerge -------------------------------------------------------------------

./plink \
--bfile $qcdir/$file.cleaned \
--bmerge $qcdir/$refname.cleanMerge \
--make-bed \
--out $qcdir/1KG.merged

# QC data

./plink \
--bfile $qcdir/1KG.merged \
--geno 0.01 \
--maf 0.01 \
--hwe 0.0001 \
--make-bed \
--out $qcdir/1KG.merged

# PCA --------------------------------------------------------------------------

./plink2 \
--bfile $qcdir/1KG.merged \
--pca \
--out $qcdir/1KG.merged

echo "PCA - Done!"

# move log files
cd ${qc_dir}

mv *.log $log
rm *~
rm *.nosex

# TODO Manually!! --------------------------------------------------------------

#Eigenvec file has FID (0s) IID PC1-10
#get ID and pop from original psam file
#FID added in ref final fam, if not just use psam

awk '{print $1, $5, $6}' $refdir/$refname.psam > $qcdir/1kG_ID2Pop.txt
#need to add FIDs back in script instead of manual
# awk '{print $1, $2, $6, $7}' $refdir/Ref_final.fam > $qcdir/1kG_ID2Pop.txt

awk '{print $1, $5, $6}' /fh/scratch/delete90/kooperberg_c/mjohnson/cseqtl/results/resources/all_hg38.psam > qc_ancestry/1kG_ID2Pop.txt
