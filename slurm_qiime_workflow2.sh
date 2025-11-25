#!/bin/sh
#
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=70000000K
#SBATCH --job-name="slurm_qiime_workflow2"
#SBATCH --output=slurm_qiime_workflow2.out
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

#downloaded silva qza files from: https://docs.qiime2.org/2023.7/data-resources/
#in /projects/academic/mbaiz/16s_yrwa/silva

#train the classifier - onyl need to do once, save in /projects/academic/mbaiz/16s_yrwa/silva/silva_classifier.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva/silva-138-99-seqs-515-806.qza \
  --i-reference-taxonomy silva/silva-138-99-tax-515-806.qza \
  --o-classifier silva/silva_classifier.qza

#use the classifier
qiime feature-classifier classify-sklearn \
  --i-classifier silva/silva_classifier.qza \
  --i-reads rep-seqs.qza \
  --p-n-jobs 4 \
  --o-classification silva_taxonomy.qza

qiime metadata tabulate \
  --m-input-file silva_taxonomy.qza \
  --o-visualization silva_taxonomy.qzv

#visualize classified sequences
qiime taxa barplot \
  --i-table feature_table.qza \
  --i-taxonomy silva_taxonomy.qza \
  --m-metadata-file metadata_yrwa_16s.tsv \
  --o-visualization taxa-bar-plots.qzv

#for decontam:
#export feature table to otu table for decontam (puts it in a folder called otu-table)
qiime tools export \
  --input-path feature_table.qza \
  --output-path otu_table_forDecontam

#convert to tsv file for decontam
biom convert -i  ./otu_table_forDecontam/feature-table.biom -o ./otu_table_forDecontam/feature-table.tsv --to-tsv

#generate a tree for decontam
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree_forDecontam.qza \
  --o-rooted-tree rooted-tree_forDecontam.qza

#export tree for decontam
qiime tools export \
  --input-path rooted-tree_forDecontam.qza \
  --output-path exported-tree_forDecontam

#don't forget to download tsv of silva_taxonomy.qzv from qiime2 view for decontam!
