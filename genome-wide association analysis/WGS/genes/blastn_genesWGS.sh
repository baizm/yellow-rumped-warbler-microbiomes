#!/bin/bash -l
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="blastn_genesWGS"
#SBATCH --output=blastn_genesWGS.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gemmaWGS/genes/blastn/
cd $WORKDIR

#load the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 blast+/2.12.0

#Using the zebra finch blast database I created of coding sequences downloaded from ensemble

#alpha genes
blastn -query ../gene_seqs_20kb_a.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_a.txt
 
#beta pcoa1 genes
blastn -query ../gene_seqs_20kb_b1.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_b1.txt

#beta pcoa2 genes 
blastn -query ../gene_seqs_20kb_b2.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_b2.txt