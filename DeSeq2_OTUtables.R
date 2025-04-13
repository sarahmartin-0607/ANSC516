##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#sample-metadata.tsv
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza

if(!dir.exists("output/taxa"))
  dir.create("output/taxa")

# to install phyloseq:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#to install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(qiime2R)
library(phyloseq)
#library(zoo)
library(tidyverse)
library(DESeq2)

##Qiime2R method of creating a phyloseq object
physeq <- qza_to_phyloseq(
  features="core-metrics-results/rarefied_table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "projectmetadata.tsv"
)

asv_table <- data.frame(otu_table(physeq), check.names = F)
metadata <- data.frame(sample_data(physeq), check.names = F)
taxonomy <- data.frame(tax_table(physeq), check.names = F)

#Clean up metadata, slightly
levels(metadata$diet)
metadata$diet = factor(metadata$diet, c("P", "A", "B", "C", "D"))
levels(metadata$diet)

#Clean up taxonomy
head(taxonomy)
tax.clean <- taxonomy

#removing NA's
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("uncl_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("uncl_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("uncl_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("uncl_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("uncl_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("uncl_",tax.clean$Genus[i], sep = "_")
  }
}

#################################################################
##Taxa barplot
#################################################################

#taxasum function will set the bar plot to 100%

#Assign our edited and formatted tables as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(asv_table), taxa_are_rows=TRUE)
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)

# Set colors for plotting
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
)

####Diet
#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "diet"  #this is the metadata column that we will use in the taxa barplot
my_column_ordered <- c("P", "A", "B", "C", "D")

abund_filter <- 0.02  # Our abundance threshold

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  physeq.taxa.mean <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.mean=mean(Abundance.average))
  
  physeq.taxa.mean <- as.data.frame(physeq.taxa.mean)
  colnames(physeq.taxa.mean)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  physeq_meta <- merge (physeq_meta, physeq.taxa.mean)
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$my_column_ordered = factor(physeq_meta_filtered[[my_column]], my_column_ordered)
  
  physeq_meta_filtered[[ml]] <- factor(physeq_meta_filtered[[ml]])
  y = tapply(physeq_meta_filtered$overall.mean, physeq_meta_filtered[[ml]], function(y) max(y))
  y = sort(y, TRUE)
  physeq_meta_filtered[[ml]] = factor(as.character(physeq_meta_filtered[[ml]]), levels=names(y))
  levels(physeq_meta_filtered[[ml]])
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = my_column_ordered, y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~.) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/taxa/", ml, "BarPlot_", my_column, ".png"), height = 5, width = 4)
}



### Donor type 
my_level <- c("Phylum", "Family", "Genus")
my_column <- "donor_type"  #this is the metadata column that we will use in the taxa barplot
my_column_ordered <- c("Healthy", "IBD")

abund_filter <- 0.02  # Our abundance threshold

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  physeq.taxa.mean <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.mean=mean(Abundance.average))
  
  physeq.taxa.mean <- as.data.frame(physeq.taxa.mean)
  colnames(physeq.taxa.mean)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  physeq_meta <- merge (physeq_meta, physeq.taxa.mean)
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$my_column_ordered = factor(physeq_meta_filtered[[my_column]], my_column_ordered)
  
  physeq_meta_filtered[[ml]] <- factor(physeq_meta_filtered[[ml]])
  y = tapply(physeq_meta_filtered$overall.mean, physeq_meta_filtered[[ml]], function(y) max(y))
  y = sort(y, TRUE)
  physeq_meta_filtered[[ml]] = factor(as.character(physeq_meta_filtered[[ml]]), levels=names(y))
  levels(physeq_meta_filtered[[ml]])
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = my_column_ordered, y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~.) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/taxa/", ml, "BarPlot_", my_column, ".png"), height = 5, width = 4)
}

### Donor
#If you want different taxonomic level, find and replace the taxonomic level listed here
my_level <- c("Phylum", "Family", "Genus")
my_column <- "donor_ind"  #this is the metadata column that we will use in the taxa barplot
my_column_ordered <- c("H-1", "H-2", "H-3", "IBD-1", "IBD-2", "IBD-3")

abund_filter <- 0.02  # Our abundance threshold

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  physeq.taxa.mean <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.mean=mean(Abundance.average))
  
  physeq.taxa.mean <- as.data.frame(physeq.taxa.mean)
  colnames(physeq.taxa.mean)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  physeq_meta <- merge (physeq_meta, physeq.taxa.mean)
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
  #str(physeq_meta_filtered)
  
  physeq_meta_filtered$my_column_ordered = factor(physeq_meta_filtered[[my_column]], my_column_ordered)
  
  physeq_meta_filtered[[ml]] <- factor(physeq_meta_filtered[[ml]])
  y = tapply(physeq_meta_filtered$overall.mean, physeq_meta_filtered[[ml]], function(y) max(y))
  y = sort(y, TRUE)
  physeq_meta_filtered[[ml]] = factor(as.character(physeq_meta_filtered[[ml]]), levels=names(y))
  levels(physeq_meta_filtered[[ml]])
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = my_column_ordered, y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~.) +
    geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = TRUE, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    xlab(my_column) +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) in at least 1 sample")) 
  ggsave(paste0("output/taxa/", ml, "BarPlot_", my_column, ".png"), height = 5, width = 4)
}



#################################################################
###Differential Abundance with DESeq2
#################################################################


#Adapted from https://joey711.github.io/phyloseq-extensions/DESeq2.html

#First load DESeq2.
#If you need help  with DESeq2 install, see this website
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html

#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1
#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq), check.names = FALSE)

OTU.clean2 <- physeq_otu_table + 1


#Now make the phyloseq object:


OTU.physeq = otu_table(as.matrix(OTU.clean2), taxa_are_rows=TRUE)

#We then merge these into an object of class phyloseq.


physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)

#The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~body.site term). The DESeq function does the rest of the testing, in this case with default testing framework, but you can actually use alternatives.

#########change the column to be analyzed
diagdds = phyloseq_to_deseq2(physeq_deseq, ~ donor_ind)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

alpha = 0.05

run_deseq2 <- function(my_factor, x, y){
  
  my_contrast <- c(my_factor, x, y)
  res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)
  
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))
  
  
  #  ###Volcano Plot
  
  #  with(res, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", xlim=c(-15,15)))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  #  with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  #  with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  
  #Let's look at the OTUs that were significantly different between the two treatment groups. The following makes a nice ggplot2 summary of the results.
  
  
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Genus order
  x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  DESeq_fig = ggplot(sigtab, aes(x=Genus, y = log2FoldChange, color=Phylum)) + 
    geom_point(size=3) + 
    ylab(paste0("(", my_contrast[2], "/", my_contrast[3], ")\n", "log2FoldChange")) +
    scale_color_manual(values = my_colors[c(4,6,8,10,12,14,16,18,20)]) +
    #ylim(0,8) +
    geom_text(color="black", x=length(unique(sigtab$Genus))-1, y=max(sigtab$log2FoldChange)-1, label=my_contrast[2], show_guide = F) +
    geom_text(color="black", x=3, y=min(sigtab$log2FoldChange)+1, label=my_contrast[3], show_guide = F) +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
  
  ggsave(paste0("output/taxa/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 5, width = 10)
}


run_deseq2("diet", "P", "A")
run_deseq2("diet", "P", "B")
run_deseq2("diet", "P", "C")
run_deseq2("diet", "P", "D")
run_deseq2("donor_type", "Healthy", "IBD")
run_deseq2("donor_ind", "H-1", "H-2")
run_deseq2("donor_ind", "H-1", "H-3")
run_deseq2("donor_ind", "H-1", "IBD-1")
run_deseq2("donor_ind", "H-1", "IBD-2")
run_deseq2("donor_ind", "H-1", "IBD-3")

