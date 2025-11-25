#!/bin/sh
#
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=60000000K
#SBATCH --job-name="slurm_qiime_workflow1"
#SBATCH --output=slurm_qiime_workflow1.out
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

#input_path contains all samples sequences(raw_reads directory)
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-path /projects/academic/mbaiz/metabarcoding_data/230110_M07914_0091_000000000-KPV9H/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path /projects/academic/mbaiz/16s_yrwa/demux-seqs.qza

qiime dada2 denoise-paired \
  --p-n-threads 0 \
  --i-demultiplexed-seqs /projects/academic/mbaiz/16s_yrwa/demux-seqs.qza \
  --p-trim-left-f 19 \
  --p-trim-left-r 20 \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 245 \
  --output-dir /projects/academic/mbaiz/16s_yrwa/denoised/ \
  --o-table feature_table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime feature-table summarize \
  --i-table feature_table.qza \
  --o-visualization feature_table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv


