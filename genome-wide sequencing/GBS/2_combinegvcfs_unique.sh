#!/bin/sh
#
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="slurm_combine_gvcfs_unique"
#SBATCH --output=slurm_combine_gvcfs_unique.out
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

##Combine individual .g.vcfs into one joint .g.vcf##
#list of sample files to join, all GBS files w/good 16s data
samples=$(find /vscratch/grp-mbaiz/yrwa_genotyping/gVCFs_unique/ | sed 's/\/\///' | grep -E 'g.vcf$' | sed 's/^/--variant /')

gatk CombineGVCFs $(echo $samples) -O ./16s.combined.unique.vcf \
  -R /projects/academic/mbaiz/mywa_genome/mywagenomev2.1.fa >& ./16s.combined.unique.vcf.log


