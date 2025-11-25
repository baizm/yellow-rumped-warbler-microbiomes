#!/bin/bash -l
#
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="impute4"
#SBATCH --output=impute4.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/gemma4/no_thin/
cd $WORKDIR

#load the necessary modules
module load gcccore/11.2.0 gatk/4.3.0.0-Java-11.0.16 
module load gcc/11.2.0 openmpi/4.1.1 plink/2.00a3.7 bcftools/1.14 vcftools/0.1.16
module load java/21.0.2


#keep only biallelic snps, remove mitochondria
gatk SelectVariants \
--variant /vscratch/grp-mbaiz/gemma4/16s.genotypes.unique.vcf \
-O /vscratch/grp-mbaiz/gemma4/no_thin/16s.genotypes.unique.bialleic.removeindels.vcf \
-select-type SNP \
--restrict-alleles-to BIALLELIC \
-XL mito \
--select "vc.isNotFiltered()"

#filter the vcf file to include ids_both (individuals with both microbiome + GBS data)
#exclude sites with >70% missiness and minor allele freq <1%
bcftools view -S ids_keep.txt -e'F_MISSING>0.7 & MAF<0.01' -O v \
  /vscratch/grp-mbaiz/gemma4/no_thin/16s.genotypes.unique.bialleic.removeindels.vcf > gemma_ids.vcf

#rename the VF03T03 sample
bcftools reheader -s vcf_names.txt -o gemma_ids_no_thin.vcf gemma_ids.vcf

#impute missing genotypes with BEAGLE_GEMMA
java -Xmx20g -jar /projects/academic/mbaiz/software/beagle.27May24.118.jar \
  gt=/vscratch/grp-mbaiz/gemma4/no_thin/gemma_ids_no_thin.vcf \
  out=/vscratch/grp-mbaiz/gemma4/no_thin/gemma_ids_no_thin.imputed

#unzip the imputed file
gunzip gemma_ids_no_thin.imputed.vcf.gz

#make binary plink file .bed
plink --vcf /vscratch/grp-mbaiz/gemma4/no_thin/gemma_ids_no_thin.imputed.vcf \
  --make-bed --chr-set 30 --allow-extra-chr --out in_file_gemma

