#!/bin/bash -l
#
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=256000M
#SBATCH --job-name="get_fst_slidingwindow"
#SBATCH --output=get_fst_slidingwindow.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/angsd_fst/
cd $WORKDIR

#set up path to angsd commands and genome assembly
angsd=/projects/academic/mbaiz/software/angsd/2023.01/software/avx512/Compiler/gcc/11.2.0/angsd/0.940/bin/angsd
realSFS=/projects/academic/mbaiz/software/angsd/2023.01/software/avx512/Compiler/gcc/11.2.0/angsd/0.940/bin/realSFS
mywa_genome=/projects/academic/mbaiz/mywa_genome/mywagenomev2.1.fa


#first calculate per pop saf for each populatoin
$angsd -b myr_ids.txt -anc $mywa_genome -out pop_myr -dosaf 1 -P 16 -gl 1
$angsd -b aud_ids.txt -anc $mywa_genome -out pop_aud -dosaf 1 -P 16 -gl 1

#calculate the 2dsfs prior
$realSFS pop_myr.saf.idx pop_aud.saf.idx -P 38 > pop_myr.pop_aud.ml

#prepare the fst for easy window analysis etc
$realSFS fst index pop_myr.saf.idx pop_aud.saf.idx -sfs pop_myr.pop_aud.ml -P 38 -fstout here

#get the global estimate
$realSFS fst stats here.fst.idx

$realSFS fst stats2 here.fst.idx -win 10000 -step 10000 > sliding_window_fst.txt

