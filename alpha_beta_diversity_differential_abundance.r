title: "Phyloseq-Microbiome"
output: html_document
date: "2025-10-25"

#' Install required packages for analyses!! 
cran_packages <- c("tidyverse", "phyloseq", "microbiome", "vegan")

for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
} 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("Biostrings", "ANCOMBC")

for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

#' require libraries
require(tidyverse)
require(phyloseq)
require(microbiome)
require(Biostrings)
require(vegan)
require(ANCOMBC)
require(phyloseq)
require(DESeq2)

#' read phyloseq obj
phyloseq<- readRDS("2025.060.phyloseq.rds")
phyloseq

#View tax, otu, and sameple data 
view(phyloseq@sam_data)
view(phyloseq@otu_table)
view(phyloseq@tax_table)


#'load metadata 
metadata <- read.csv("metadata_second_batch.csv", 
                     stringsAsFactors = TRUE)
rownames(metadata) <- metadata$sample_id
metadata$sample_id <- NULL
all(rownames(metadata) %in% sample_names(phyloseq))

#add metadata to phyloseq obj
sample_data(phyloseq) <- sample_data(metadata)


#split Arthritis and thyroid 
Phyloseq_Thyroid <- subset_samples(phyloseq, Group == "Thyroid")
Phyloseq_Thyroid <- prune_taxa(taxa_sums(Phyloseq_Thyroid) > 0, Phyloseq_Thyroid)

#view the component of object
view(Phyloseq_Thyroid@sam_data)
view(Phyloseq_Thyroid@otu_table)
view(Phyloseq_Thyroid@tax_table)


#' Arthritis project
Phyloseq_arthritis <- subset_samples(phyloseq, Group == "Arthritis")
Phyloseq_arthritis <- prune_taxa(taxa_sums(Phyloseq_arthritis) > 0, Phyloseq_arthritis)

#view the component of object
view(Phyloseq_arthritis@sam_data)
view(Phyloseq_arthritis@otu_table)
view(Phyloseq_arthritis@tax_table)


#'Alpha diversity (Arthritis)
require(tidyverse)
alpha_div <- read.table("2025.060.rarefied..no_of_classified_reads_and_alpha_diversity.txt", 
                        header = TRUE, 
                        sep = "\t", stringsAsFactors = FALSE)
head(alpha_div)


#'Change column name "X" to "sample_id"

colnames(alpha_div)[colnames(alpha_div) == "X"] <- "sample_id"

#' metadata for alpha diversity 
metadata <- read.csv("metadata_second_batch.csv", 
                     stringsAsFactors = FALSE)
alpha_meta <- merge(alpha_div, metadata, by = "sample_id")
view(alpha_meta)

# write csv file or excel if you like to do it in prism 
write.csv(alpha_meta, file = "Alpha_Diversity.csv")

# Select only arthritis to plot 
alpha_diversity_arthritis<- alpha_meta %>% filter(Group == "Arthritis")
alpha_diversity_arthritis<- alpha_diversity_arthritis %>% filter(Treatment %in% c("HC", 
                                                                                  "Pre-ICI",
                                                                                  "ICI X 2", 
                                                                                  "ICI X 3", 
                                                                                  "ICI X 4"))

#'Order the group
alpha_diversity_arthritis$Treatment<- factor(alpha_diversity_arthritis$Treatment, 
                                             levels = c("HC", 
                                                        "Pre-ICI",
                                                        "ICI X 2", 
                                                        "ICI X 3", 
                                                        "ICI X 4"))

write.csv(alpha_diversity_arthritis, file = "alpha_diversity_arthritis.csv")

#' ggplot
alpha.diversity.plot<- ggplot(alpha_diversity_arthritis, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2) +
  theme_classic() +
  labs(title = "Alpha Diversity") +
 theme(plot.title = element_text(hjust = 0.5, face= "bold", size= 9), 
        axis.title.y = element_text(size= 14, face = "bold"),
        axis.text.x = element_text(size= 10,
                                   face = "bold",
                                   angle = 45, 
                                   vjust=0.5), 
        axis.title.x = element_blank(), 
       legend.position = "none") 

alpha.diversity.plot
 
anova<- aov(Shannon ~ Treatment,alpha_diversity_arthritis)
summary(anova)

#' save high quality plot 
ggsave("alpha.diversity_second.tiff", plot = alpha.diversity.plot, dpi = 300, width = 3, height = 2, units = "in")


#' Beta diversity (split the two projects)
Beta_bray <- read.table("2025.060.rarefied.bray_relative_dist.txt", 
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

#'  Read metadata
metadata <- read.csv("metadata_second_batch.csv", stringsAsFactors = FALSE)

#'  Match sample IDs between matrix and metadata
common_samples <- intersect(rownames(Beta_bray), metadata$sample_id)

# Subset matrix and metadata to matching samples only
beta_mat <- Beta_bray[common_samples, common_samples]
metadata <- metadata[metadata$sample_id %in% common_samples, ]

#' Split by Group
thyroid_samples <- metadata$sample_id[metadata$Group == "Thyroid"]
arthritis_samples <- metadata$sample_id[metadata$Group == "Arthritis"]

beta_thyroid <- beta_mat[thyroid_samples, thyroid_samples]
beta_arthritis <- beta_mat[arthritis_samples, arthritis_samples]

#' Save each subset
write.csv(beta_thyroid, "beta_Thyroid.csv", quote = FALSE)
write.csv(beta_arthritis, "beta_Arthritis.csv", quote = FALSE)

# Optional: check dimensions
cat("Thyroid matrix dimensions: ", dim(beta_thyroid), "\n")
cat("Arthritis matrix dimensions: ", dim(beta_arthritis), "\n")

#' Beta Diversity for arthritis 

#' Read bata-diversity file for arthritis back 
beta_arthritis <- read.csv("beta_Arthritis.csv",
                           row.names = 1, 
                           check.names = FALSE)

metadata<- read.csv("metadata_arthritis_second_batch.csv", 
                    stringsAsFactors = TRUE)
rownames(metadata) <- metadata$Sample

#' filter some groups

metadata<- metadata %>% filter(split %in% c("HC", "Pre_ICI", "ICI X 2", "ICI X 3", "ICI X 4"))

#' Match metadat and the distance matrix 
common_samples <- intersect(rownames(beta_arthritis), metadata$Sample)
beta_arthritis <- beta_arthritis[common_samples, common_samples]
metadata <- metadata[common_samples, ]

#' Run PCOA

#' Convert to distance object
dist_arthritis <- as.dist(beta_arthritis)

#' Perform Principal Coordinates Analysis (PCoA)
pcoa_arthritis <- ape::pcoa(dist_arthritis)

#' extract coordinate 
pcoa_points <- as.data.frame(pcoa_arthritis$vectors[, 1:2])
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$Sample <- rownames(pcoa_points)

# Merge with metadata
pcoa_points <- left_join(pcoa_points, metadata, by = "Sample")

# Calculate % variance explained
variance_explained <- round(100 * pcoa_arthritis$values$Relative_eig[1:2], 1)
library(dplyr)
library(ggplot2)
library(ape)

#' Calculate centroids (mean position for each group)
centroids <- pcoa_points %>%
  group_by(split) %>%
  summarize(
    centroid_x = mean(PCoA1),
    centroid_y = mean(PCoA2)
  )

#' Merge centroid info back to each sample
pcoa_points <- left_join(pcoa_points, centroids, by = "split")

#' Plot: draw lines from each point to its group centroid
pcoa_points$split<- factor(pcoa_points$split, levels = c("HC", "Pre_ICI", "ICI X 2", "ICI X 3", "ICI X 4"))
beta.diversity <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = split)) +
  geom_segment(aes(xend = centroid_x, yend = centroid_y), alpha = 0.7, linewidth = 0.7) +  # lines to centroid
  geom_point(size = 3) +  # sample points
  theme_bw(base_size = 14) +
  theme_classic() +
  xlab(paste0("PCoA1 (", variance_explained[1], "%)")) +
  ylab(paste0("PCoA2 (", variance_explained[2], "%)")) +
  ggtitle("Beta Diversity (PCoA)") +
  theme(
    legend.position = "right",
    plot.title = element_blank(), 
    legend.title = element_blank(), 
    axis.text = element_text(size= 12, face = "bold"), 
    axis.title = element_text(size= 12, face = "bold"), 
    legend.text = element_text(size= 12, face = "bold"))

beta.diversity


#save high quality plot 
ggsave("beta.diversity.second.batch.tiff", plot = beta.diversity, dpi = 300, width = 3.9, height = 2, units = "in")

#' PCOA test 
adonis2(as.dist(beta_arthritis) ~ group, data = metadata, permutations = 999)

#' Relative abundance arhritis project 
phyloseq_rel_RA <- transform_sample_counts(Phyloseq_arthritis, function(x) x / sum(x))
view(Phyloseq_arthritis@sam_data)


#' Aggregate to Genus
phyloseq_family <- tax_glom(phyloseq_rel_RA, taxrank = "Family")

# Convert to data frame and summarize by group
# Melt to long format
df_family <- psmelt(phyloseq_family)

#' Summarize mean relative abundance per group
df_group_mean <- df_family %>%
  group_by(Treatment, Family) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  ungroup()

#' Identify top 10
# Identify top 15 genera overall
top_genera <- df_group_mean %>%
  group_by(Family) %>%
  summarise(total = sum(mean_abundance)) %>%
  top_n(10, total) %>%
  pull(Family)

#' Filter to only those top genera
df_group_top <- df_group_mean %>%
  filter(Family %in% top_genera)

#' select relevant group 
relevant<- df_group_top %>% filter(Treatment %in% c("HC", "Pre-ICI", "ICI X 2", "ICI X 3", "ICI X 4"))

#' Replot with top genera
Top10_bacteria<- ggplot(relevant, aes(x = Treatment, y = mean_abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  ggtitle("Top 10 Genera by Average Relative Abundance") +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
Top10_bacteria

#' Define color palette (CIBERSORT-inspired, 10 distinct colors)
cibersort_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#E31A1C"
)


#' Replot with custom colors
relevant$Treatment<- factor(relevant$Treatment, levels = c("HC", "Pre-ICI", "ICI X 2", "ICI X 3", "ICI X 4"))
Top10_bacteria <- ggplot(relevant, aes(x = Treatment, y = mean_abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size= 11, face = "bold"),
        axis.text.y = element_text(
                                   hjust = 1,
                                   size = 11, 
                                   face = "bold"),
        axis.text.x = element_text(angle= 45,
                                   hjust = 1,
                                   size = 11, 
                                   face = "bold"),
        plot.title = element_text(hjust = 0.5), 
        legend.title = element_text(face= "bold", size= 12), 
        legend.text = element_text(face = "bold", size= 7)) +
  scale_fill_manual(values = cibersort_colors)
Top10_bacteria

# ---- Normalize each Day_factor to fractions ----
relevant_norm <- relevant %>%
  group_by(Treatment) %>%
  mutate(fraction = (mean_abundance / sum(mean_abundance, na.rm = TRUE)) * 100)

#' ---- Define color palette ----
cibersort_colors <- c(
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#E31A1C"
)

#' ---- Plot ----
relevant$Treatment<- factor(relevant$Treatment, levels = c("HC", "Pre-ICI", "ICI X 2", "ICI X 3", "ICI X 4"))
Top10_bacteria <- ggplot(relevant_norm, aes(x = Treatment, y = fraction, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_fill_manual(values = cibersort_colors) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(face = "bold", size = 8)
  ) +
  labs(y = "Relative Abundance (%)")
Top10_bacteria

ggsave("Top15_bacteria_family.tiff", plot = Top10_bacteria, dpi = 300, width = 3.7, height = 3.5, units = "in")

#' Order

#' Aggregate to Genus
phyloseq_order <- tax_glom(phyloseq_rel_RA, taxrank = "Order")

#' Convert to data frame and summarize by group
# Melt to long format
df_order <- psmelt(phyloseq_order)

#' Summarize mean relative abundance per group
df_group_mean <- df_order %>%
  group_by(Day_factor, Order) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  ungroup()

#' Identify top 15 genera overall
top_genera <- df_group_mean %>%
  group_by(Order) %>%
  summarise(total = sum(mean_abundance)) %>%
  top_n(10, total) %>%
  pull(Order)

#' Filter to only those top genera
df_group_top <- df_group_mean %>%
  filter(Order %in% top_genera)

#' select relevant group 
relevant<- df_group_top %>% filter(Day_factor %in% c("D19", "D23", "D26", "D29"))

cibersort_colors <- c(
  "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
  "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"
)

#' Replot with custom colors
relevant$Day_factor<- factor(relevant$Day_factor, levels = c("D19", "D23", "D26", "D29"))
Top10_bacteria <- ggplot(relevant, aes(x = Day_factor, y = mean_abundance, fill = Order)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  ggtitle("Top 10 Genera by Average Relative Abundance") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = cibersort_colors)

Top10_bacteria

ggsave("Top10_bacteria_order.tiff", plot = Top10_bacteria, dpi = 300, width = 7, height = 4, units = "in") 

#' Class 
#' Aggregate to Genus
phyloseq_class <- tax_glom(phyloseq_rel_RA, taxrank = "Class")

#' Convert to data frame and summarize by group
#' Melt to long format
df_class <- psmelt(phyloseq_class)

# Summarize mean relative abundance per group
df_group_mean <- df_class %>%
  group_by(Treatment, Class) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  ungroup()

#' Identify top 15 genera overall
top_genera <- df_group_mean %>%
  group_by(Class) %>%
  summarise(total = sum(mean_abundance)) %>%
  top_n(10, total) %>%
  pull(Class)

# Filter to only those top genera
df_group_top <- df_group_mean %>%
  filter(Class %in% top_genera)

#' select relevant group 
relevant<- df_group_top %>% filter(Treatment %in% c("CII-D19", "CII-ICI.D23", "CII-ICI.D26", "CII-ICI.D29"))
cibersort_colors <- c(
  "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
  "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"
)

#' Replot with custom colors
relevant$Treatment<- factor(relevant$Treatment, levels = c("HC", "CII-D19", "CII-ICI.D23", "CII-ICI.D26", "CII-ICI.D29"))
Top10_bacteria <- ggplot(relevant, aes(x = Treatment, y = mean_abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  ggtitle("Top 10 Genera by Average Relative Abundance") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = cibersort_colors)
Top10_bacteria

ggsave("Top10_bacteria_class.tiff", plot = Top10_bacteria, dpi = 300, width = 7, height = 4, units = "in")
 
# Phylum 
# Aggregate to Genus
phyloseq_phylum <- tax_glom(phyloseq_rel_RA, taxrank = "Phylum")

# Convert to data frame and summarize by group
# Melt to long format
df_phylum <- psmelt(phyloseq_phylum)

# Summarize mean relative abundance per group
df_group_mean <- df_phylum %>%
  group_by(Treatment, Phylum) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  ungroup()

#' Identify top 15 genera overall
top_genera <- df_group_mean %>%
  group_by(Phylum) %>%
  summarise(total = sum(mean_abundance)) %>%
  top_n(10, total) %>%
  pull(Phylum)

#' Filter to only those top genera
df_group_top <- df_group_mean %>%
  filter(Phylum %in% top_genera)

#' select relevant group 
relevant<- df_group_top %>% filter(Treatment %in% c("HC", "CII-D19", "CII-ICI.D23", "CII-ICI.D26", "CII-ICI.D29"))
cibersort_colors <- c(
  "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
  "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6"
)

#' Replot with custom colors
relevant$Treatment<- factor(relevant$Treatment, levels = c("HC", "CII-D19", "CII-ICI.D23", "CII-ICI.D26", "CII-ICI.D29"))
Top10_bacteria <- ggplot(relevant, aes(x = Treatment, y = mean_abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw(base_size = 14) +
  ggtitle("Top 10 Genera by Average Relative Abundance") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = cibersort_colors)
Top10_bacteria
ggsave("Top10_bacteria_phylum.tiff", plot = Top10_bacteria, dpi = 300, width = 7, height = 4, units = "in")

#' Relative abundance using DESeq2 
require(DESeq2)
view(Phyloseq_arthritis@tax_table)
Phyloseq_arthritis <- subset_samples(
  Phyloseq_arthritis,
  !(Treatment %in% c("CII-D23", "CII-D26", "CII-D29", "Exclusion"))
)
view(Phyloseq_arthritis@sam_data)
dds <- phyloseq_to_deseq2(Phyloseq_arthritis, ~ Treatment)
dds <- DESeq(dds)
res_d23<- results(dds, contrast = c("Treatment", "CII-D19", "HC"))
res_sigd23 <- res_d23[which(res_d23$padj < 0.05), ]
res_sigd23 <- res_sigd23[order(res_sigd23$log2FoldChange, decreasing = TRUE), ]
head(res_sigd23)

#' Add Taxa info
tax <- as.data.frame(tax_table(Phyloseq_arthritis))
res_sig_tax <- cbind(as(res_sigd23, "data.frame"), tax[rownames(res_sigd23), ])
all(rownames(res_sigd23) %in% rownames(tax))  # return True meaning it annotate correctly!! 
head(res_sig_tax)

#In case you would like to save as excel or csv file
library(openxlsx)
write.xlsx(res_sig_tax, file = "XXXX.xlsx")
write.csv(res_sig_tax, file = "Pre-ICI_vs_HC.csv")


preICIvs_HC<- read.csv("Pre-ICI_vs_HC_.csv", stringsAsFactors = TRUE)
eg <- ggplot(preICIvs_HC, 
             aes(x = log2FoldChange, 
                 y = reorder(Genus, log2FoldChange), 
                 color = padj, 
                 shape = top, 
                 size = abs(log2FoldChange))) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#BB0000", high = "#0000FF", 
                       name = "adj. p-value") +
  theme_classic()+
  labs(x = "log2 Fold Change", 
       y = "Genus",
       color = "adj. p-value",
       shape = "Group",
       size = "log2FC",
       title = "Differential Abundance") +
  theme(axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
eg

# Anther way 
library(ggplot2)
library(dplyr)
library(forcats)

# Ensure Genus is a factor for proper ordering
preICIvs_HC <- preICIvs_HC %>%
  mutate(Genus = fct_reorder(Genus, log2FoldChange))

# Create a 'sign' or 'group' column if needed
preICIvs_HC$sign <- preICIvs_HC$top  # top = "Pre-ICI" or "HC"

# Plot
eg <- ggplot(preICIvs_HC, 
             aes(x = log2FoldChange, 
                 y = Genus, 
                 color = padj, 
                 size = abs(log2FoldChange))) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#BB0000", high = "#0000FF", 
                       name = "adj. p-value") +
  facet_grid(. ~ sign, scales = "free_y", 
             space = "free_y",
          ) + # split by group
  theme_bw() +
  labs(x = "log2 Fold Change", 
       y = "Genus",
       color = "adj. p-value",
       size = "log2FC",
       title = "Differential Abundance by Group") +
  theme(axis.text.y = element_text(face = "bold", size= 13),
        axis.text.x= element_text(face= "bold", size=10),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold")) + 
  xlim(-25, 30)

eg 
ggsave("Pre-IcI_vs_HC.tiff", plot = eg, dpi = 300, width = 7.5, height = 3.3, units = "in")

#' Volcano plot 
res_df <- as.data.frame(res_d23)
res_df$taxa_id <- rownames(res_df)  # Add ASV/OTU ID

#Add taxonomy (e.g., Genus) to results
tax_df <- as.data.frame(tax_table(Phyloseq_arthritis))
tax_df$taxa_id <- rownames(tax_df)

# Merge tax info with DESeq2 results
res_annotated <- merge(res_df, tax_df, by = "taxa_id")

#Add significance label
res_annotated$Significant <- ifelse(res_annotated$padj < 0.05, "Significance", "Not significance")

#volcano plot 
library(ggplot2)
library(ggrepel)  # for nicer text labels

iciarthritisd23_vs_preici<- ggplot(res_annotated, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  geom_text_repel(
    data = subset(res_annotated, padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = Genus),
    size = 3.5,
    max.overlaps = 10
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face= "bold", size= 14), 
        axis.title.y = element_text(size= 14, face = "bold"),
        axis.text = element_text(size= 13, face = "bold"), 
        axis.title.x = element_text(size= 12, face = "bold")) +
  scale_color_manual(values = c("Not significance" = "grey", "Significance" = "red")) +
  labs(title = "CII-ICI-D29 vs. CII-ICI-D26",
       x = "log2 Fold Change",
       y = "-log10 adjusted p-value")
iciarthritisd23_vs_preici

#save high quality plot 
ggsave("ICI_D29_ICI_D26_second.tiff", plot = iciarthritisd23_vs_preici, dpi = 300, width = 7, height = 4, units = "in")

#' Time course analysis 

#' exclude some groups 
view(Phyloseq_arthritis@sam_data)
phyloseq_filtered <- subset_samples(
  Phyloseq_arthritis,
  !(split %in% c("HC", "Non-IA", "Exclusion"))
)

#' Convert phyloseq to DESeq2, including Treatment and Day_nu
dds <- phyloseq_to_deseq2(phyloseq_filtered, ~ split + Day_numeric)

#' Make sure Treatment is a factor
dds$split <- factor(dds$split)

#' Set No_arthritis as the reference level
dds$split <- relevel(dds$split, ref = "CII")

#' Run DESeq2
dds <- DESeq(dds)

#' Check available coefficients
resultsNames(dds)

#' Main effect of Treatment (ICI_IA vs No_arthritis), adjusting for day
res <- results(dds, contrast = c("split", "IA", "CII"))
# Order by adjusted p-value
res <- res[order(res$padj), ]
res_sig <- res[which(res$padj < 0.05), ]
res_sig <- res_sig[order(res_sig$log2FoldChange, decreasing = TRUE), ]
head(res_sig)

tax <- as.data.frame(tax_table(phyloseq_filtered))
res_tax <- cbind(as(res, "data.frame"), tax[rownames(res), ])
all(rownames(res) %in% rownames(tax))  # return True meaning it annotate correctly!! 

#'save data 
res_ordered <- res_tax[order(-res_tax$log2FoldChange, -res_tax$baseMean), ]
write.csv(res_ordered, file = "CII_ICIvs_CII.all.csv")

#' plot 
library(ggplot2)
library(dplyr)
library(forcats)

#' Read data 
ICI_IA_vs_CII<- read.csv("ICI_IA_vs_CII.csv", 
                         stringsAsFactors = TRUE)


#' Ensure Genus is a factor for proper ordering
ICI_IA_vs_CII <- ICI_IA_vs_CII %>%
  mutate(Genus = fct_reorder(Genus, log2FoldChange))

#' Create a 'sign' or 'group' column if needed
ICI_IA_vs_CII$sign <- ICI_IA_vs_CII$top  # top = "Pre-ICI" or "HC"
ICI_IA_vs_CII$sign<- factor(ICI_IA_vs_CII$sign, levels = c("Pre-ICI", "ICI-IA"))

#' Plot
eg <- ggplot(ICI_IA_vs_CII, 
             aes(x = log2FoldChange, 
                 y = Genus, 
                 color = padj, 
                 size = abs(log2FoldChange))) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#BB0000", high = "#0000FF", 
                       name = "adj. p-value") +
  facet_grid(. ~ sign, scales = "free_y", 
             space = "free_y",
          ) + # split by group
  theme_bw() +
  labs(x = "log2 Fold Change", 
       y = "Genus",
       color = "adj. p-value",
       size = "log2FC",
       title = "Differential Abundance by Group") +
  theme(axis.text.y = element_text(face = "bold"),
        axis.text.x= element_text(face= "bold", size=10),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold")) + 
  xlim(-25, 35)

eg

#' new code 
ICI_IA_vs_CII <- read.csv("ICI_IA_vs_CII.csv", stringsAsFactors = TRUE)

#' Ensure Genus is a factor for proper ordering
ICI_IA_vs_CII <- ICI_IA_vs_CII %>%
  mutate(Genus = fct_reorder(Genus, log2FoldChange),
         sign = factor(top, levels = c("ICI-IA", "Pre-ICI")))  # set facet order here

#' Plot
eg <- ggplot(ICI_IA_vs_CII, 
             aes(x = log2FoldChange, 
                 y = Genus, 
                 color = padj, 
                 size = abs(log2FoldChange))) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#BB0000", high = "#0000FF", name = "adj. p-value") +
  facet_grid(. ~ sign, scales = "free_y", space = "free_y") +  # use 'sign' with proper levels
  theme_bw() +
  labs(x = "log2 Fold Change", 
       y = "Genus",
       color = "adj. p-value",
       size = "log2FC",
       title = "Differential Abundance by Group") +
  theme(axis.text.y = element_text(face = "bold", size= 13),
        axis.text.x= element_text(face= "bold", size=10),
        axis.title = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 14, face = "bold")) +
  xlim(-25, 35)
eg 
ggsave("ICI-IA_vs_pre-ici.tiff", plot = eg, dpi = 300, width = 7.5, height = 3.3, units = "in")



#' Add Taxa info
tax <- as.data.frame(tax_table(phyloseq_filtered))
res_sig_tax <- cbind(as(res_sig, "data.frame"), tax[rownames(res_sig), ])
all(rownames(res_sig) %in% rownames(tax))  # return True meaning it annotate correctly!! 
head(res_sig_tax)
write.csv(res_sig_tax, file = "ICI_IA_vs_CII.csv")

#' Visualization 
library(ggplot2)
library(dplyr)

res <- read.csv("ICI_IA_vs_CII.csv")

#' Filter for significant ASVs (padj < 0.05) and remove empty Genus
res_sig <- res %>%
  filter(padj < 0.05 & Genus != "")

#' Horizontal bar plot showing only Genus
ggplot(res_sig, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = padj)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars
  scale_fill_gradient(low = "#bb0c00", high = "#F0F0F0", name = "padj") +
  theme_classic(base_size = 14) +
  labs(x = "Genus", y = "log2FoldChange", title = "Differentially Abundant Bacteria") +
  theme(axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))

#'Lippo plot 
ggplot(res_sig, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange)) +
  geom_segment(aes(x = Genus, xend = Genus, y = 0, yend = log2FoldChange),
               color = "grey") +
  geom_point(aes(color = padj), size = 4) +
  coord_flip() +
  scale_color_gradient(low = "#bb0c00", high = "#F0F0F0") +
  theme_classic(base_size = 14) +
  labs(x = "Genus", y = "log2 Fold Change", title = "Differential Abundance (ICI-IA vs. Pre-ICI)" +
  theme(axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))

#' dot plot 
ggplot(res_sig, aes(x = log2FoldChange, y = reorder(Genus, log2FoldChange), 
                    size = baseMean, color = padj)) +
  geom_point() +
  scale_color_gradient(low = "#bb0c00", high = "#F0F0F0") +
  theme_classic(base_size = 14) +
  labs(x = "log2 Fold Change", y = "Genus", title = "Differential Abundance Dot Plot", size="Mean abundance") +
  theme(plot.title = element_text(hjust = 0.5))

#' LIppo
deg<- ggplot(res_sig, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange)) +
  geom_segment(aes(x = Genus, xend = Genus, y = 0, yend = log2FoldChange),
               color = "grey") +
  geom_point(aes(color = padj), size = 4) +
  coord_flip() +
  scale_color_gradient(low = "#bb0c00", high = "#F0F0F0") +
  theme_classic(base_size = 14) +
  labs(x = "Genus", y = "log2 Fold Change",
       title = "Differential Abundance (ICI-IA vs. Pre-ICI)") +
  theme(axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold", 
                                  size= 12),
        legend.title = element_text(size= 12, face = "bold"),
        legend.text = element_text(size= 12, face = "bold"),
        plot.title = element_blank())
deg

#' save high quality plot 
ggsave("ICI_IA.vs_pre_ICI.tiff", plot = deg, dpi = 300, width = 6, height = 3, units = "in")

#' Volcano plot!! 
res_df <- as.data.frame(res)
res_df$taxa_id <- rownames(res_df)  # Add ASV/OTU ID

#' Add taxonomy (e.g., Genus) to results
tax_df <- as.data.frame(tax_table(phyloseq_filtered))
tax_df$taxa_id <- rownames(tax_df)

#' Merge tax info with DESeq2 results
res_annotated <- merge(res_df, tax_df, by = "taxa_id")

#' Add significance label
res_annotated$Significant <- ifelse(res_annotated$padj < 0.05, "Significance", "Not significance")

#' volcano plot 
library(ggplot2)
library(ggrepel)  # for nicer text labels
iciarthritisd23_vs_preici<- ggplot(res_annotated, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significant)) +
  geom_text_repel(
    data = subset(res_annotated, padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = Genus),
    size = 3.5,
    max.overlaps = 10
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face= "bold", size= 14), 
        axis.title.y = element_text(size= 14, face = "bold"),
        axis.text = element_text(size= 13, face = "bold"), 
        axis.title.x = element_text(size= 12, face = "bold")) +
  scale_color_manual(values = c("Not significance" = "grey", "Significance" = "red")) +
  labs(title = "ICI-IA vs. pre-ICI",
       x = "log2 Fold Change",
       y = "-log10 adjusted p-value")
iciarthritisd23_vs_preici

#' save high quality plot 
ggsave("ICI-IAvs_CII_second.tiff", plot = iciarthritisd23_vs_preici, dpi = 300, width = 7, height = 4, units = "in")

#' ANCOM-BC to validate DEseq2 output 
# Global model !! 
set.seed(123)
output <- ancombc2(
  Phyloseq_arthritis,
  fix_formula = "Treatment",   # variable to adjust for (can be the same)
  group = "Treatment",         # REQUIRED for pairwise or global tests
  p_adj_method = "holm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 2,
  global = FALSE,
  pairwise = TRUE,
  verbose = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100)
)

res_pair <- output$res_pair
colnames(res_pair)

#' More Restricted model
set.seed(123)
Phyloseq_subset <- subset_samples(
  Phyloseq_arthritis,
  Treatment %in% c("CII-D19", "CII-ICI.D23", "CII-ICI.D26", "CII-ICI.D29")
)

#'  Relevel Treatment so CII-D19 is the reference
sample_data(Phyloseq_subset)$Treatment <- factor(
  sample_data(Phyloseq_subset)$Treatment,
  levels = c("CII-D19", "CII-ICI.D23", "CII-ICI.D26", "CII-ICI.D29")
)

#' Run ANCOM-BC2
output <- ancombc2(
  Phyloseq_subset,
  fix_formula = "Treatment",
  group = "Treatment",         # required for pairwise comparisons
  p_adj_method = "holm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 2,
  global = FALSE,              # skip global test
  pairwise = TRUE,             # do pairwise tests
  verbose = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100)
)

#'  Extract pairwise results
res_pair <- output$res_pair

#'  Check available LFC columns
colnames(res_pair)



#model with just two groups == 

#' select only the columns for CII-ICI.D23 vs CII-D19
df_sig_D23 <- res_pair[, c(
  "taxon",
  "lfc_TreatmentCII-ICI.D23",
  "q_TreatmentCII-ICI.D23"
)]

# filter taxa that are significant in this comparison
df_sig_D23_filtered <- df_sig_D23[
  df_sig_D23$`q_TreatmentCII-ICI.D23` < 0.1,
]

#' view top significant taxa
head(df_sig_D23_filtered)

# Map taxa to genus
tax_table_df <- as.data.frame(tax_table(Phyloseq_arthritis))

#' Keep only the columns we need (assume column "Genus" exists)
tax_table_df$taxon <- rownames(tax_table_df)
tax_table_genus <- tax_table_df[, c("taxon", "Genus")]

#' Merge genus info with your ANCOM-BC results
df_sig_genus <- merge(
  df_sig_D23_filtered,
  tax_table_genus,
  by = "taxon",
  all.x = TRUE
)

#' Optional: summarize by genus (average log-fold change if multiple ASVs per genus)
library(dplyr)
df_genus_summary <- df_sig_genus %>%
  group_by(Genus) %>%
  summarize(
    mean_lfc = mean(`lfc_TreatmentCII-ICI.D23`),
    min_q = min(`q_TreatmentCII-ICI.D23`),
    .groups = "drop"
  )

#' ggplot 
ggplot(df_genus_summary, aes(x = reorder(Genus, mean_lfc), y = mean_lfc)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Mean log-fold change (CII-ICI.D23 vs CII-D19)",
    title = "Differential abundance at genus level"
  ) +
  theme_minimal()



#' Start with just two groups 
groups_of_interest <- c("CII-D19", "CII-ICI.D23")
Phyloseq_subset <- subset_samples(Phyloseq_arthritis, Treatment %in% groups_of_interest)

# Prune taxa with zero counts across all samples
Phyloseq_subset <- prune_taxa(taxa_sums(Phyloseq_subset) > 0, Phyloseq_subset)

# Remove taxa with zero variance across samples
taxa_var <- apply(otu_table(Phyloseq_subset), 1, var)
Phyloseq_subset <- prune_taxa(taxa_var > 0, Phyloseq_subset)

# Set Treatment factor with reference level

sample_data(Phyloseq_subset)$Treatment <- factor(
  sample_data(Phyloseq_subset)$Treatment,
  levels = groups_of_interest  # first level = reference
)

#' Check
levels(sample_data(Phyloseq_subset)$Treatment)

# Run ANCOM-BC
output <- ancombc2(
  Phyloseq_subset,
  fix_formula = "Treatment",
  group = "Treatment",
  p_adj_method = "holm",
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 2,
  global = FALSE,
  pairwise = TRUE,
  verbose = TRUE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100)
)

# Extract pairwise results
res_pair <- output$res_pair
colnames(res_pair)

#' Extract significant taxa for the pair
#' =============================================
#' Only need columns of interest
df_sig <- res_pair[, c(
  "taxon",
  "lfc_TreatmentCII-ICI.D23",
  "q_TreatmentCII-ICI.D23"
)]

# Filter significant taxa
df_sig_filtered <- df_sig %>%
  filter(q_TreatmentCII-ICI.D23 < 0.05)

# View top significant taxa
head(df_sig_filtered)

#  Optional: Merge taxonomy to genus level
tax_table_df <- as.data.frame(tax_table(Phyloseq_subset))
tax_table_df$taxon <- rownames(tax_table_df)

df_sig_filtered_genus <- df_sig_filtered %>%
  left_join(tax_table_df[, c("taxon", "Genus")], by = "taxon")

head(df_sig_filtered_genus)

#  Optional: Visualization (bar plot example)

library(ggplot2)

ggplot(df_sig_filtered_genus, aes(x = Genus, y = lfc_TreatmentCII.ICI.D23)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Genus",
    y = "Log fold change (CII-ICI.D23 vs CII-D19)",
    title = "Significant taxa between two groups"
  )

# Session info
sessionInfo()
