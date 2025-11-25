#!/bin/bash -l
#
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="run_gemmaWGS"
#SBATCH --output=run_gemmaWGS.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/
cd $WORKDIR

#load the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 plink/2.00a3.7 gemma/0.98.5 bcftools/1.14


gemma -bfile /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/in_file_gemmaWGS -gk 1 -o relatedness

# linear mixed models, using covariates file
gemma -bfile /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/in_file_gemmaWGS \
  -k /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/output/relatedness.cXX.txt \
  -n 1 \
  -c gemmaWGS_covars_in.txt \
  -lmm 4 -o gemma_outWGS_a

gemma -bfile /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/in_file_gemmaWGS \
  -k /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/output/relatedness.cXX.txt \
  -n 14 \
  -c gemmaWGS_covars_in.txt \
  -lmm 4 -o gemma_outWGS_b1

gemma -bfile /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/in_file_gemmaWGS \
  -k /vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/output/relatedness.cXX.txt \
  -n 15 \
  -c gemmaWGS_covars_in.txt \
  -lmm 4 -o gemma_outWGS_b2

