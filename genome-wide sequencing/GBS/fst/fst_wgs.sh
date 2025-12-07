#!/bin/bash -l
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="fst_wgs"
#SBATCH --output=fst_wgs.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gemma2/fst_genotypes/
cd $WORKDIR

#load the necessary modules
module load gcc/11.2.0 vcftools/0.1.16 bcftools/1.14

#--extract genotypes genotypes

#get genotypes from the WGS genotype file
#first compress using bgzip for bcftools
bcftools view /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/gemma_idsWGS.vcf -Oz -o gemma_idsWGS_bgzip.vcf.gz
bcftools index gemma_idsWGS_bgzip.vcf.gz

bcftools query -f '%CHROM %POS  %REF  %ALT [ %TGT]\n' -R regions_gemmaWGS.txt /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/gemma_idsWGS_bgzip.vcf.gz > gemmaWGS_genotypes_WGSdataset.txt
#write sample names to file for R
bcftools query -l /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/gemma_idsWGS_bgzip.vcf.gz > gemma_idsWGS_bgzip_ids.txt

