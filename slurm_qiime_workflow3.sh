#!/bin/sh
#
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=50000000K
#SBATCH --job-name="slurm_qiime_workflow3"
#SBATCH --output=slurm_qiime_workflow3.out
#SBATCH --mail-user=mbaiz@buffalo.edu
#SBATCH --mail-type=all
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --reservation=ubhpc-future

#Set up the path to the working directory
WORKDIR=/projects/academic/mbaiz/16s_yrwa/
cd $WORKDIR

#load the necessary modules
module load anaconda3
#activate the qiime environment
source activate qiime2-2023.7

#pick up here after using R decontam to ID contaminant ASVs to remove
#filter out mithochondrial, chloroplast and unassigned ASVs, AND eukaryota
qiime taxa filter-table \
  --i-table feature_table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --p-exclude mitochondria,chloroplast,unassigned,eukaryota \
  --o-filtered-table feature_table_no_MtCpUnEuk.qza

#filter out contaminant ASVs from decontam
qiime feature-table filter-features \
  --i-table feature_table_no_MtCpUnEuk.qza \
  --p-exclude-ids \
  --m-metadata-file asv_to_remove.txt \
  --o-filtered-table feature_table_no_contam.qza

#summarize new feature tabulate
qiime feature-table summarize \
  --i-table feature_table_no_contam.qza \
  --o-visualization feature_table_no_contam.qzv

#filter seqs to remove mt,cp,un,euk,contam for phylogeny
qiime feature-table filter-seqs \
  --i-data rep-seqs.qza \
  --i-table feature_table_no_contam.qza \
  --o-filtered-data rep-seqs_no_contam.qza

#visualize, should have  sequences (subtract cp,mt,un,euk,contam) - !
qiime feature-table tabulate-seqs \
  --i-data rep-seqs_no_contam.qza \
  --o-visualization rep-seqs_no_contam.qzv

#generate a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs_no_contam.qza \
  --o-alignment aligned-rep-seqs_no_contam.qza \
  --o-masked-alignment masked-aligned-rep-seqs_no_contam.qza \
  --o-tree unrooted-tree_final.qza \
  --o-rooted-tree rooted-tree_final.qza

#produce alpha-rarefaction curve
qiime diversity alpha-rarefaction \
  --i-table feature_table_no_contam.qza \
  --i-phylogeny rooted-tree_final.qza \
  --p-max-depth 21000 \
  --m-metadata-file metadata_yrwa_16s.tsv \
  --o-visualization alpha-rarefaction.qzv

#export feature table to otu table (puts it in a folder called out-table.biom)
qiime tools export \
  --input-path feature_table_no_contam.qza \
  --output-path otu_table_final

#convert to tsv file for phyloseq
biom convert -i  ./otu_table_final/feature-table.biom -o ./otu_table_final/feature-table.tsv --to-tsv

#export tree for phyloseq
qiime tools export \
  --input-path rooted-tree_final.qza \
  --output-path exported-tree_final

