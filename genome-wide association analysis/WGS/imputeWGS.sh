#!/bin/bash -l
#
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="imputeWGS"
#SBATCH --output=imputeWGS.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/
cd $WORKDIR

#load the necessary modules
module load gcccore/11.2.0
module load gatk/4.3.0.0-Java-11.0.16
module load gcc/11.2.0
module load bcftools/1.14
module load java/21.0.2

###---run this script to prepare files for GEMMA

#filter VCF to include only biallelic SNPs. 
gatk SelectVariants \
--variant /projects/academic/mbaiz/YRWA_GBS_2023/wgs_vcf/marcella.recode.vcf \
-O /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/marcella.recode.bialleic.removeindels.vcf \
-select-type SNP \
--restrict-alleles-to BIALLELIC \
--select "vc.isNotFiltered()"

#filter the vcf file to include ids_both (individuals with both microbiome + WGS data)
#exclude sites with >50% missiness and minor allele freq <1%
bcftools view -S ids_bothWGS.txt -e'F_MISSING>0.5 & MAF<0.01' -O v \
  /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/marcella.recode.bialleic.removeindels.vcf > gemma_idsWGS.vcf

#impute missing genotypes with BEAGLE_GEMMA
java -Xmx20g -jar /projects/academic/mbaiz/software/beagle.27May24.118.jar \
  gt=/vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/gemma_idsWGS.vcf \
  out=/vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/gemma_idsWGS.imputed

#unzip the imputed file
gunzip gemma_idsWGS.imputed.vcf.gz

#make binary plink file .bed
/projects/academic/mbaiz/software/plink2 --vcf /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/gemma_idsWGS.imputed.vcf \
  --make-bed --chr-set 30 --allow-extra-chr --out in_file_gemmaWGS
