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

#Get the databases to classify the sequences
#wget \
 # -O "gg-13-8-99-515-806-nb-classifier.qza" \
 # "https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes/gg-13-8-99-515-806-nb-classifier.qza"

#echo "Classify and get .qzv"
#qiime feature-classifier classify-sklearn \
 # --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
 # --i-reads rep-seqs.qza \
 # --o-classification taxonomy.qza

#qiime metadata tabulate \
 # --m-input-file taxonomy.qza \
 # --o-visualization taxonomy.qzv

#echo "Make a stacked barplot"
#qiime taxa barplot \
  #--i-table table.qza \
 # --i-taxonomy taxonomy.qza \
 # --m-metadata-file projectmetadata.tsv \
 # --o-visualization taxa-bar-plots.qzv


######
#Diet ancom
echo "Filter away some sample types and do ancom"
#qiime feature-table filter-samples \
 # --i-table table.qza \
 # --m-metadata-file projectmetadata.tsv \
 # --p-where donor_type::IBD \
 # --o-filtered-table IBD-table.qza

#qiime composition add-pseudocount \
 # --i-table IBD-table.qza \
 # --o-composition-table comp-IBD-table.qza

qiime composition ancombc \
  --i-table table.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-formula 'donor_type' \
  --p-reference-levels donor_type::Healthy \
  --o-differentials signficance/ancombc-donor.qza

qiime composition tabulate \
  --i-data signficance/ancombc-donor.qza \
  --o-visualization signficance/ancombc-donor.qzv

qiime composition da-barplot \
  --i-data signficance/ancombc-donor.qza \
  --p-significance-threshold 0.01 \
  --o-visualization signficance/da-barplot-donor.qzv

echo "Collapse levels and ancom-bc again."

qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table signficance/table-l6.qza

qiime composition ancombc \
  --i-table signficance/table-l6.qza \
  --m-metadata-file projectmetadata.tsv \
  --p-formula 'donor_type' \
  --p-reference-levels donor_type::Healthy \
  --o-differentials signficance/l6-ancombc-donor.qza

qiime composition da-barplot \
  --i-data signficance/l6-ancombc-donor.qza \
  --p-significance-threshold 0.01 \
  --p-level-delimiter ';' \
  --o-visualization signficance/l6-da-barplot-donor.qzv

#######
#Healthy samples vs diet

echo "Filter away some sample types and do ancom"
#qiime feature-table filter-samples \
 # --i-table table.qza \
 # --m-metadata-file projectmetadata.tsv \
 # --p-where "[donor_type]='Healthy'" \
 # --o-filtered-table Healthy-table.qza

qiime composition add-pseudocount \
 # --i-table Healthy-table.qza \
 # --o-composition-table comp-Healthy-table.qza
 
#qiime composition ancombc \
 # --i-table Healthy-table.qza \
 # --m-metadata-file projectmetadata.tsv \
 # --p-formula 'diet' \
 # --p-reference-levels diet::P \
 # --o-differentials ancombc-Healthydiet.qza

#qiime composition tabulate \
 # --i-data ancombc-Healthydiet.qza \
 # --o-visualization ancombc-Healthydiet.qzv

#qiime composition da-barplot \
 # --i-data ancombc-Healthydiet.qza \
 # --p-significance-threshold 0.001 \
 # --o-visualization da-barplot-Healthydiet.qzv

echo "Collapse levels and ancom-bc again."

#qiime taxa collapse \
 # --i-table Healthy-table.qza \
 # --i-taxonomy taxonomy.qza \
 # --p-level 6 \
 # --o-collapsed-table Healthy-table-l6.qza

#qiime composition ancombc \
 # --i-table Healthy-table-l6.qza \
  #--m-metadata-file projectmetadata.tsv \
  #--p-formula 'diet' \
  #--p-reference-levels diet::P \
  #--o-differentials l6-ancombc-Healthydiet.qza

#qiime composition da-barplot \
 # --i-data l6-ancombc-Healthydiet.qza \
 # --p-significance-threshold 0.001 \
 # --p-level-delimiter ';' \
 # --o-visualization l6-da-barplot-Healthydiet.qzv

###########
#Donor type vs diet

echo "Filter away some sample types and do ancom"
#qiime feature-table filter-samples \
 # --i-table table.qza \
 # --m-metadata-file projectmetadata.tsv \
 # --p-where "[donor_ind] IN ('IBD-1','IBD-2','IBD-3','H-1','H-2','H-3')" \
 # --o-filtered-table Donor-table.qza

#qiime composition add-pseudocount \
 # --i-table Donor-table.qza \
 # --o-composition-table comp-Donor-table.qza
  
#qiime composition ancombc \
#  --i-table table.qza \
#  --m-metadata-file projectmetadata.tsv \
#  --p-formula 'donor_ind' \
#  --p-reference-levels donor_ind::H1 \
#  --o-differentials ancombc-donor-ind.qza

#qiime composition tabulate \
#  --i-data ancombc-donor-ind.qza \
#  --o-visualization ancombc-donor-ind.qzv

#qiime composition da-barplot \
#  --i-data ancombc-donor-ind.qza \
#  --p-significance-threshold 0.001 \
#  --o-visualization da-barplot-donor-ind.qzv

#echo "Collapse levels and ancom-bc again."

#qiime taxa collapse \
 # --i-table table.qza \
 # --i-taxonomy taxonomy.qza \
 # --p-level 6 \
 # --o-collapsed-table Donor-table-l6.qza

#qiime composition ancombc \
 # --i-table Donor-table-l6.qza \
 # --m-metadata-file projectmetadata.tsv \
 # --p-formula 'donor_ind' \
 # --p-reference-levels donor_ind::H1 \
 # --o-differentials l6-ancombc-donor.qza

#qiime composition da-barplot \
 # --i-data l6-ancombc-donor.qza \
 # --p-significance-threshold 0.001 \
 # --p-level-delimiter ';' \
 # --o-visualization l6-da-barplot-donor.qzv

echo "End time"
date +"%d %B %Y %H:%M:%S"
