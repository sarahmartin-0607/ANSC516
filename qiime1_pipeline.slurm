#!/bin/sh -l

#SBATCH --nodes=1                    ##number of requested nodes (a node is a group of CPUs)
#SBATCH -n 20                        ##number of requested CPUs
#SBATCH --time=04:00:00              ##time requested
#SBATCH --job-name ProjectQiime1       ##name of job
#SBATCH -A microbiome                ##name of the queue you are using. Could be scholar or microbiome if you are on snyder

#########################################################################

echo "Start time"
date +"%d %B %Y %H:%M:%S"

#step 1, cd into the proper directory. This directory must already exist
cd $RCAC_SCRATCH
cd ANSCFinalProject/data_download/

pwd

echo "loading bioinfo"
module load bioinfo

echo "loading Qiime/2-2023.9"
module load Qiime/2-2023.9

#step 3, Run the code you want
echo "Step 1 of qiime pipeline"

echo "Data Import"
#Make the sequences into a qiime recognizable format

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path ManifestProject.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path demux-untrimmed.qza 

echo "Trimming Primers"
#cut off primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences demux-untrimmed.qza \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r GGACTACNVGGGTWTCTAAT \
--p-match-adapter-wildcards \
--p-match-read-wildcards \
--p-cores 12 \
--o-trimmed-sequences demux.qza \
--verbose

#demux summary
qiime demux summarize \
  --i-data demux.qza \
  --o-visualization demux.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"

