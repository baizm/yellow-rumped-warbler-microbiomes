#!/bin/bash -l
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="fst_gemma4"
#SBATCH --output=fst_gemma4.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/gemma4/no_thin/fst/
cd $WORKDIR

#load the necessary modules
module load gcc/11.2.0 vcftools/0.1.16 bcftools/1.14

#-- calculate FST ---
vcftools --vcf /vscratch/grp-mbaiz/gemma4/no_thin/gemma_ids_no_thin.vcf \
	--weir-fst-pop aud_ids2.txt \
	--weir-fst-pop myr_ids2.txt \
	--out gemma_ids_no_thin_fst.txt

#--extract genotype information
#first compress using bgzip for bcftools
bgzip -c /vscratch/grp-mbaiz/gemma4/no_thin/gemma_ids_no_thin.vcf > gemma_ids_no_thin.vcf.gz
bcftools view gemma_ids_no_thin.vcf.gz -Oz -o gemma_ids_no_thin_bgzip.vcf.gz
bcftools index gemma_ids_no_thin_bgzip.vcf.gz

#now extract genotypes for significant SNPs from gemma4
bcftools query -f '%CHROM %POS  %REF  %ALT [ %TGT]\n' -R regions_gemma4.txt gemma_ids_no_thin_bgzip.vcf.gz > gemma4_genotypes.txt

#write sample names to file for R
bcftools query -l gemma_ids_no_thin_bgzip.vcf.gz > gemma_ids_no_thin_bgzip_ids.txt
