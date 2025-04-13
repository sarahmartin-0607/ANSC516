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


qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-where "[donor_ind]='H-1'" \
  --o-filtered-table DonorTables/H1-table.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-where "[donor_ind]='H-2'" \
  --o-filtered-table DonorTables/H2-table.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-where "[donor_ind]='H-3'" \
  --o-filtered-table DonorTables/H3-table.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-where "[donor_ind]='IBD-1'" \
  --o-filtered-table DonorTables/IBD1-table.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-where "[donor_ind]='IBD-2'" \
  --o-filtered-table DonorTables/IBD2-table.qza

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-where "[donor_ind]='IBD-3'" \
  --o-filtered-table DonorTables/IBD3-table.qza

echo "End time"
date +"%d %B %Y %H:%M:%S"
