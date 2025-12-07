#!/bin/sh
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=20GB
#SBATCH --job-name="run_admixture_unique"
#SBATCH --output=run_admixture_unique.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/admixture/unique/
cd $WORKDIR

module load admixture/1.3.0

# convert VCF to plink .bed, format required for ADMIXTURE
/projects/academic/mbaiz/software/plink2 --vcf /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/biallelic_DP10_maf_thin_hwe_50.recode.vcf \
  --make-bed --chr-set 30 --allow-extra-chr 0 --allow-extra-chr 0 --out in_file_admixture

# Ks we will test
for K in 1 2 3 4 5 6; do admixture --cv in_file_admixture.bed $K -j8 | tee log${K}.out; done

grep -h CV log*.out > cv_k1-6.txt


