#!/bin/bash -l
#
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=256000M
#SBATCH --job-name="blastn_genes_gemma4"
#SBATCH --output=blastn_genes_gemma4.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#Set up the path to the working directory
WORKDIR=/vscratch/grp-mbaiz/yrwa_genotyping/gemma4/genes_gemma4/blastn/
cd $WORKDIR

#load the necessary modules
module load gcc/11.2.0 openmpi/4.1.1 blast+/2.12.0

#Using the zebra finch blast database I created of coding sequences downloaded from ensemble

#beta1 genes
blastn -query ../gene_seqs_20kb_b1.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_b1.txt

#beij genes
blastn -query ../gene_seqs_20kb_beij.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_beij.txt

#micro genes
blastn -query ../gene_seqs_20kb_micro.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_micro.txt
 
#run 10/5/25 on login node
#all genes
blastn -query ../gene_seqs_20kb_all.fasta -db /projects/academic/mbaiz/taegut_genome/Taeniopygia_guttata.bTaeGut1_v1.p.cds.all.fa \
 -num_threads 4 -outfmt 7 -out taegut_genes_all.txt
