#!/bin/sh
#
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=60000000K
#SBATCH --job-name="slurm_joint_genotyping_unique"
#SBATCH --output=slurm_joint_genotyping_unique.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/
cd $WORKDIR

#load the necessary modules
module load gcccore/11.2.0
module load gatk/4.3.0.0-Java-11.0.16

#picking up after CombineGVCFs step
gatk GenotypeGVCFs -R /projects/academic/mbaiz/mywa_genome/mywagenomev2.1.fa \
  -V /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/16s.combined.unique.vcf \
  -O /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/16s.genotypes.unique.vcf >& /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/16s.genotypes.unique.vcf.log


