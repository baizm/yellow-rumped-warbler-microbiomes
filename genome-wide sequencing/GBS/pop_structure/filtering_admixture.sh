#!/bin/sh
#
#SBATCH --time=2:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=1000M
#SBATCH --job-name="filtering_admixture"
#SBATCH --output=filtering_admixture.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/
cd $WORKDIR

module load gcc/11.2.0
module load bcftools/1.14
module load vcftools/0.1.16

#filter raw vcf to keep only bi-allelic markers
bcftools view -m2 -M2 /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/16s.genotypes.unique.vcf > biallelic.vcf

#filter vcf to keep only sites with minimum mean depth of 10 across individuals
vcftools --vcf biallelic.vcf --min-meanDP 10 --recode --out biallelic_DP10

#filter vcf to get rid of monomorphic sites by setting minor allele frequency low
vcftools --vcf biallelic_DP10.recode.vcf --maf 0.01 --recode --out biallelic_DP10_maf

#filter vcf to reduce effects of linkage by thinning snps within 200bp
vcftools --vcf biallelic_DP10_maf.recode.vcf --thin 200 --recode --out biallelic_DP10_maf_thin

#filter vcf to exclude sites out of HWE at p<0.001
vcftools --vcf biallelic_DP10_maf_thin.recode.vcf --hwe 0.001 --recode --out biallelic_DP10_maf_thin_hwe

#filter vcf to exclude sites not present in 50% of individuals
vcftools --vcf biallelic_DP10_maf_thin_hwe.recode.vcf --max-missing 0.50 --recode --out biallelic_DP10_maf_thin_hwe_50

