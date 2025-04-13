#Install the packages, IF YOU NEED TO :)
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)

#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza

getwd()
#make working directory where core-metrix is at
#How to load a file into R
#sep is the sperator (tab separated), header =T (true),
metadata2 <- read.delim("projectmetadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
#square bracket- indexing [row,column]
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column, (-) means remove row 1
metadata2 <- metadata2[-1,]
metadata <- metadata2

#row.names(assigns a new row name)
row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

#read_qza is a function
bc_PCoA<-read_qza("core-metrics-results/bray_curtis_pcoa_results.qza")
wUF <- read_qza("core-metrics-results/weighted_unifrac_pcoa_results.qza")
uUF <- read_qza("core-metrics-results/unweighted_unifrac_pcoa_results.qza")
jc_PCoA <- read_qza("core-metrics-results/jaccard_pcoa_results.qza")


bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

wUF_meta <- wUF$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

uUF_meta <- uUF$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

jc_PCoA_meta <- jc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

# Now we are going to make an ordination plot
#diet
my_column <- "Diet"
diet_colors <- c("Blue", "Green", "Gray", "Pink", "Purple")
ggplot(bc_meta, aes(x=PC1, y=PC2, color=diet)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=diet_colors, name = "Diet")
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(wUF_meta, aes(x=PC1, y=PC2, color=diet)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=diet_colors, name = "Diet")
ggsave(paste0("output/wUF-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(uUF_meta, aes(x=PC1, y=PC2, color=diet)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=diet_colors, name = "Diet")
ggsave(paste0("output/uUF-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(jc_PCoA_meta, aes(x=PC1, y=PC2, color=diet)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=diet_colors, name = "Diet")
ggsave(paste0("output/jc-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

#donor_type
donor_colors <- c("Blue","Pink")
my_column <- "Donor_Type"
ggplot(bc_meta, aes(x=PC1, y=PC2, color=donor_type)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Disease Status")
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(wUF_meta, aes(x=PC1, y=PC2, color=donor_type)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Disease Status")
ggsave(paste0("output/wUF-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(uUF_meta, aes(x=PC1, y=PC2, color=donor_type)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Disease Statue")
ggsave(paste0("output/uUF-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(jc_PCoA_meta, aes(x=PC1, y=PC2, color=donor_type)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Disease Status")
ggsave(paste0("output/jc-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches


#donor_ind
donor_colors <- c("Blue","Pink", "Purple","Red","Yellow","Orange")
my_column <- "Donor_Ind"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=donor_ind)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Donor")
ggsave(paste0("output/BC-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(wUF_meta, aes(x=PC1, y=PC2, color=donor_ind)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Donor")
ggsave(paste0("output/wUF-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(uUF_meta, aes(x=PC1, y=PC2, color=donor_ind)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Donor")
ggsave(paste0("output/uUF-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches

ggplot(jc_PCoA_meta, aes(x=PC1, y=PC2, color=donor_ind)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=donor_colors, name = "Donor")
ggsave(paste0("output/jc-basic_", my_column,".tiff"), height=2, width=3, device="tiff") # save a PDF 3 inches by 4 inches


#Run some PERMANOVAs
bc_dist_mat<-read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(bc_dm ~ diet, data = metadata_sub)

write.table(PERMANOVA_out,"output/Diet_Adonis_overall.csv",sep=",", row.names = TRUE) 

bc_dist_mat<-read_qza("core-metrics-results/bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$SampleID),]
rownames(bc_dm) == metadata_sub$SampleID ## all these values need to be "TRUE"

PERMANOVA_out <- adonis2(bc_dm ~ donor_ind, data = metadata_sub)

write.table(PERMANOVA_out,"output/Donorind_Adonis_overall.csv",sep=",", row.names = TRUE) 
