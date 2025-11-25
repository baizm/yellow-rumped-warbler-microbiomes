#!/bin/bash -l
#
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="run_gemma4"
#SBATCH --output=run_gemma4.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/gemma4/no_thin
cd $WORKDIR

#load the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 gemma/0.98.5


#generate the relatedness matrix
gemma -bfile /vscratch/grp-mbaiz/gemma4/no_thin/in_file_gemma -gk 1 -o relatedness

# linear mixed models, using covariates file
gemma -bfile /vscratch/grp-mbaiz/gemma4/no_thin/in_file_gemma \
  -k /vscratch/grp-mbaiz/gemma4/no_thin/output/relatedness.cXX.txt \
  -n 1 \
  -c gemma_covars_in.txt \
  -lmm 4 -o gemma_out4_a
  
gemma -bfile /vscratch/grp-mbaiz/gemma4/no_thin/in_file_gemma \
 -k /vscratch/grp-mbaiz/gemma4/no_thin/output/relatedness.cXX.txt \
 -n 14 \
 -c gemma_covars_in.txt \
 -lmm 4 -o gemma_out4_b1

gemma -bfile /vscratch/grp-mbaiz/gemma4/no_thin/in_file_gemma \
 -k /vscratch/grp-mbaiz/gemma4/no_thin/output/relatedness.cXX.txt \
 -n 15 \
 -c gemma_covars_in.txt \
 -lmm 4 -o gemma_out4_b2
 
gemma -bfile /vscratch/grp-mbaiz/gemma4/no_thin/in_file_gemma \
 -k /vscratch/grp-mbaiz/gemma4/no_thin/output/relatedness.cXX.txt \
 -n 16 \
 -c gemma_covars_in.txt \
 -lmm 4 -o gemma_out4_beij

gemma -bfile /vscratch/grp-mbaiz/gemma4/no_thin/in_file_gemma \
 -k /vscratch/grp-mbaiz/gemma4/no_thin/output/relatedness.cXX.txt \
 -n 17 \
 -c gemma_covars_in.txt \
 -lmm 4 -o gemma_out4_micro
 
