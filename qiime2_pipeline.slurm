#!/bin/sh -l

#SBATCH --nodes=1                    ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20                        ##number of requested CPUs
#SBATCH --time=04:00:00              ##time requested
#SBATCH --job-name ProjectQiime2       ##name of job
#SBATCH -A microbiome                ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist
cd $RCAC_SCRATCH
cd ANSCFinalProject/data_download/qiimeProject

pwd


#De-noise. Error corrections.

echo "dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-min-overlap 8 \
  --p-n-threads 20 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

#Create a .qzv from the output of denoise
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file projectmetadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

echo "Create a phylogenetic tree. (Phylogenetic method)"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

#Alpha rarefaction
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 30000 \
  --m-metadata-file projectmetadata.tsv \
  --o-visualization alpha-rarefaction.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"
