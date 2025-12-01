#!/bin/bash -l
#
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="gemma.zzzzz"
#SBATCH --output=gemma.zzzzz.stdout
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#load the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 gemma/0.98.5

# Print current date
date

# Change directory to where you submitted the job from, so that relative paths resolve properly
cd /vscratch/grp-mbaiz/perm_GBS

#zzzzz is a placeholder that will be swapped for numbers from 1 to 100 by the wrapper loop
#create iteration-specific directories
mkdir zzzzz
mkdir zzzzz/output

#copy gemma input files in plink bed format to the iteration-specific directories
cp in_file_gemma* zzzzz/
cp gemma_covars_in.txt zzzzz/
cp output/relatedness.cXX.txt zzzzz/output/

echo start
date

#randomly permute phenotype values across individuals, and store in iteration-specific directory
shuf pheno | paste nopheno - > zzzzz/in_file_gemma.fam

cd zzzzz

echo shuffle
date

#run gemma gwas on permuted phenotype values
gemma -bfile in_file_gemma -c gemma_covars_in.txt -k ./output/relatedness.cXX.txt -lmm 4 -o gem.perm

echo assoc
date

#extract p_lrt p-values from gemma results
cd output
cut -f 14 gem.perm.assoc.txt > plrtval.txt

#extract the lowest p-value for this iteration; append to an output text file
cat plrtval.txt | awk  'BEGIN{min=1}{if(($1)<min)  min=($1)}END {print min}' >> ../../smallestpval.2.txt

echo pvals
date

#clean up temporary iteration-specific directory
cd ../..
rm -r zzzzz/

hostname
