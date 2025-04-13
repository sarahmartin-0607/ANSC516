#alpha and beta diversity

#The first step is very important. You need to set your working 
#directory. Just like in unix we have to `cd` to where our data is, the 
#same thing is true for R.

#move to the core-metrics file and set tht as the working directory

list.files()

###Files
# We will use the following files created using the qiime2 moving pictures tutorial but for our own data.

# core-metrics-results/evenness_vector.qza (alpha diversity)
# core-metrics-results/faith_pd_vector.qza (alpha diversity)
# core-metrics-results/observed_features_vector.qza (alpha diversity)
# core-metrics-results/shannon_vector.qza (alpha diversity)
# sample-metadata.tsv (#my file for the project, projectmetadata.tsv)

#download packages if needed
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
install.packages("ggpubr")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

## Load Packages
library(tidyverse)
library(qiime2R)
library(ggpubr)

#reading in metadata
meta<-read_q2metadata("projectmetadata.tsv")
str(meta)

#reading in alpha diversity 
evenness = read_qza("core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\

#Merging all alpha diversities in meta
alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(x=alpha_diversity, y=observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(x=alpha_diversity, y=shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID


# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$pielou_e)
shapiro.test(meta$observed_features)

#Plots the qq-plot for residuals
ggqqplot(meta$shannon_entropy, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$pielou_e, title = "Evenness")
ggqqplot(meta$observed_features, title = "Observed Features")

## **Non-normally distributed metrics**
#diet
kruskal.test(faith_pd ~ diet, data=meta)
kruskal.test(shannon_entropy ~ diet, data=meta)
kruskal.test(pielou_evenness ~ diet, data=meta)
kruskal.test(observed_features ~ diet, data=meta)

pairwise.wilcox.test(meta$faith_pd, meta$diet, p.adjust.method="BH")
pairwise.wilcox.test(meta$shannon_entropy, meta$diet, p.adjust.method="BH")
pairwise.wilcox.test(meta$pielou_evenness, meta$diet, p.adjust.method="BH")
pairwise.wilcox.test(meta$observed_features, meta$diet, p.adjust.method="BH")

#IBD vs Healthy
kruskal.test(observed_features ~ donor_type, data=meta)
kruskal.test(shannon_entropy ~ donor_type, data=meta)
kruskal.test(pielou_evenness ~ donor_type, data=meta)
kruskal.test(observed_features ~ donor_type, data=meta)

pairwise.wilcox.test(meta$faith_pd, meta$donor_type, p.adjust.method="BH")
pairwise.wilcox.test(meta$shannon_entropy, meta$donor_type, p.adjust.method="BH")
pairwise.wilcox.test(meta$pielou_evenness, meta$donor_type, p.adjust.method="BH")
pairwise.wilcox.test(meta$observed_features, meta$donor_type, p.adjust.method="BH")

#Donor Idnetifier
kruskal.test(observed_features ~ donor_ind, data=meta)
kruskal.test(shannon_entropy ~ donor_ind, data=meta)
kruskal.test(pielou_evenness ~ donor_ind, data=meta)
kruskal.test(observed_features ~ donor_ind, data=meta)

pairwise.wilcox.test(meta$faith_pd, meta$donor_ind, p.adjust.method="BH")
pairwise.wilcox.test(meta$shannon_entropy, meta$donor_ind, p.adjust.method="BH")
pairwise.wilcox.test(meta$pielou_evenness, meta$donor_ind, p.adjust.method="BH")
pairwise.wilcox.test(meta$observed_features, meta$donor_ind, p.adjust.method="BH")

##ggplots
#diets
faith_pd_boxplot <- ggplot(meta, aes(diet, faith_pd)) + 
  geom_boxplot(aes(color = diet)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/pd_diet.png", faith_pd_boxplot, height = 3, width = 3)

shannon_boxplot <- ggplot(meta, aes(diet, shannon_entropy)) + 
  geom_boxplot(aes(color = diet)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Diversity", x = "") 
ggsave("output/shannon_diet.png", shannon_boxplot, height = 3, width = 3)

evenness_boxplot <- ggplot(meta, aes(diet, pielou_evenness)) + 
  geom_boxplot(aes(color = diet)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's Evenness", x = "") 
ggsave("output/evenness_diet.png", evenness_boxplot, height = 3, width = 3)

observed_boxplot <- ggplot(meta, aes(diet, observed_features)) + 
  geom_boxplot(aes(color = diet)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("output/observed_diet.png", observed_boxplot, height = 3, width = 3)

#Healthy v IBD
faith_pd_boxplot <- ggplot(meta, aes(donor_type, faith_pd)) + 
  geom_boxplot(aes(color = donor_type)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/pd_donortype.png", faith_pd_boxplot, height = 3, width = 3)

shannon_boxplot <- ggplot(meta, aes(donor_type, shannon_entropy)) + 
  geom_boxplot(aes(color = donor_type)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Diversity", x = "") 
ggsave("output/shannon_donortype.png", shannon_boxplot, height = 3, width = 3)

evenness_boxplot <- ggplot(meta, aes(donor_type, pielou_evenness)) + 
  geom_boxplot(aes(color = donor_type)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's Evenness", x = "") 
ggsave("output/evenness_donor_type.png", evenness_boxplot, height = 3, width = 3)

observed_boxplot <- ggplot(meta, aes(donor_type, observed_features)) + 
  geom_boxplot(aes(color = donor_type)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("output/observed_donor_type.png", observed_boxplot, height = 3, width = 3)

#Donor Identifier
faith_pd_boxplot <- ggplot(meta, aes(donor_ind, faith_pd)) + 
  geom_boxplot(aes(color = donor_ind)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/pd_donorind.png", faith_pd_boxplot, height = 3, width = 3)

shannon_boxplot <- ggplot(meta, aes(donor_ind, shannon_entropy)) + 
  geom_boxplot(aes(color = donor_ind)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon Diversity", x = "") 
ggsave("output/shannon_donorind.png", shannon_boxplot, height = 3, width = 3)

evenness_boxplot <- ggplot(meta, aes(donor_ind, pielou_evenness)) + 
  geom_boxplot(aes(color = donor_ind)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's Evenness", x = "") 
ggsave("output/evenness_donor_ind.png", evenness_boxplot, height = 3, width = 3)

observed_boxplot <- ggplot(meta, aes(donor_ind, observed_features)) + 
  geom_boxplot(aes(color = donor_ind)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed Features", x = "") 
ggsave("output/observed_donor_ind.png", observed_boxplot, height = 3, width = 3)
