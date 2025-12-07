#!/bin/sh
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=20GB
#SBATCH --job-name="run_admixture_myrtles"
#SBATCH --output=run_admixture_myrtles.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/admixture/unique/substructure_myrtles/
cd $WORKDIR

module load  gcc/11.2.0
module load bcftools/1.14
module load admixture/1.3.0

# filter VCF to include only allopatric myrtles
bcftools view -S allo_ids.txt /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/biallelic_DP10_maf_thin_hwe_50.recode.vcf > substructure_myrtles.vcf

# convert VCF to plink .bed, format required for ADMIXTURE
/projects/academic/mbaiz/software/plink2 --vcf substructure_myrtles.vcf \
  --make-bed --chr-set 30 --allow-extra-chr 0 --out myrtles_admixture

# Ks t0 test
for K in 1 2 3 4; do admixture --cv myrtles_admixture.bed $K -j8 | tee log${K}.out; done

grep -h CV log*.out > cv_k1-4.txt
