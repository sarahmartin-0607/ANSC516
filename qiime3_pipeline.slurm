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


#Calculate diversity metrics. You will need to designate a sampling depth. You might want to consult the rarifaction output prior to this step.
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 20765 \
  --m-metadata-file projectmetadata.tsv \
  --output-dir core-metrics-results

#Alpha Diversity Significance
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file projectmetadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file projectmetadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/shannon_vector.qza \
  --m-metadata-file projectmetadata.tsv \
  --o-visualization core-metrics-results/shannon-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/observed_features_vector.qza \
  --m-metadata-file projectmetadata.tsv \
  --o-visualization core-metrics-results/observed_features-group-significance.qzv

#Beta Diversity Significance
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file projectmetadata.tsv \
  --m-metadata-column diet \
  --o-visualization core-metrics-results/unweighted-unifrac-diet-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file projectmetadata.tsv \
  --m-metadata-column donor_type \
  --o-visualization core-metrics-results/unweighted-unifrac-donor-type-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file projectmetadata.tsv \
  --m-metadata-column donor_ind \
  --o-visualization core-metrics-results/unweighted-unifrac-donor-ind-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file projectmetadata.tsv \
  --m-metadata-column diet \
  --o-visualization core-metrics-results/bray-curtis-diet-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file projectmetadata.tsv \
  --m-metadata-column donor_type \
  --o-visualization core-metrics-results/bray-curtis-donor-type-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
  --m-metadata-file projectmetadata.tsv \
  --m-metadata-column donor_ind \
  --o-visualization core-metrics-results/bray-curtis-donor-ind-significance.qzv \
  --p-pairwise


#Make ordination
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-custom-axes donor_type \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-donor-type.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-custom-axes diet \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-donor-type.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-custom-axes donor_type \
  --o-visualization core-metrics-results/bray-curtis-emperor-donor-type.qzv


qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-custom-axes diet \
  --o-visualization core-metrics-results/bray-curtis-emperor-diet.qzv



echo "End time"
date +"%d %B %Y %H:%M:%S"
