
title: "Bacterial metabolic Pathway analysis inferred from 16S rRNAseq"
output: html_document
date: "2025-11-07"

#' Load all required libraries 
BiocManager::install("DESeq2")
require(IRanges)
require(Tax4Fun2)
require(phyloseq)
require(Biostrings)
require(biomformat)
require(DESeq2)
require(tidyverse)
require(pheatmap)
require(data.table)

#' Export ASV sequences (FASTA)
asv_seqs <- refseq(Phyloseq_arthritis)
Biostrings::writeXStringSet(asv_seqs, "asv_sequences.fasta")
otu_mat <- as(otu_table(Phyloseq_arthritis), "matrix")

#' Export abundance table as BIOM
biom_obj <- make_biom(data = otu_mat)
write_biom(biom_obj, "asv_table.biom")

#' Turn to bash terminal 
module load miniforge3 
eval "$(/risapps/rhel8/miniforge3/24.5.0-0/bin/conda shell.bash hook)"

conda create -n picrust2 -c conda-forge -c bioconda picrust2=2.5.1
conda activate picrust2
picrust2_pipeline.py \
    -s asv_sequences.fasta \
    -i asv_table.biom \
    -o picrust2_output \
    -p 8

#' Read PICRUSt2 unstratified pathway table and align with metadata
path_abun <- fread("picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")
path_abun

#' since the first column is pathway !! 
path_abun <- as.data.frame(path_abun)
rownames(path_abun) <- path_abun$pathway
path_abun$pathway <- NULL
path_abun

#' Ensure sample names match your metadata
metadata <- as.data.frame(sample_data(Phyloseq_arthritis))
all(colnames(path_abun) %in% rownames(metadata))

#' Build a phyloseq object for pathway data
OTU <- otu_table(as.matrix(path_abun), taxa_are_rows = TRUE)
phy_path <- phyloseq(OTU, sample_data(metadata))
phy_path@sam_data

#' Differential pathway analysis with DESeq2
phyloseq_filtered <- subset_samples(
  phy_path,
  !(split %in% c("HC", "Non-IA", "Exclusion"))
)

dds <- phyloseq_to_deseq2(phyloseq_filtered, ~ split + Day_numeric)
dds$split <- factor(dds$split)
dds$split <- relevel(dds$split, ref = "CII")

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast = c("split", "IA", "CII"))
res_df <- as.data.frame(res)
res_df$pathway <- rownames(res_df)

res_df <- res_df[order(res_df$padj), ]
res_sig <- subset(res_df, padj < 0.05)
res_sig <- res_sig[order(-res_sig$log2FoldChange), ]

write.csv(res_df, "CII_ICIvs_CII.all.csv")
write.csv(res_sig, "CII_ICIvs_CII.sig.csv")

#' Another option
install.packages("broom.mixed")
library(tidyverse)
library(lme4)
library(broom.mixed)
library(ggplot2)

#' Read PWY abundances per sample
pwy_abun <- read_tsv("picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")

#' First column usually is pathway ID
colnames(pwy_abun)[1] <- "PWY"

#' Pivot longer for easier handling
pwy_long <- pwy_abun %>%
  pivot_longer(-PWY, names_to = "Sample", values_to = "Abundance")

#' Merge with metadata
metadata <- read.csv("metadata_arthritis_second_batch.csv", stringsAsFactors = TRUE)

#' Keep only relevant treatment groups
pwy_long <- pwy_long %>%
  left_join(metadata, by = "Sample")

pwy_long<- pwy_long %>%
  filter(!Treatment %in% c("HC", "Non-IA", "exclusion"))

#' Check factor levels
unique(pwy_long$Treatment)
#' Make sure it matches what you filter in the model
pwy_long$Treatment <- factor(pwy_long$Treatment, levels = c("Pre-ICI", "IA"))

#' Filter pathways present in at least 2 samples per treatment:
pwy_filtered <- pwy_long %>%
  group_by(PWY, Treatment) %>%
  filter(n() >= 2) %>%
  ungroup()

pwy_filtered <- pwy_filtered %>%
  mutate(Abundance_log = log10(Abundance))

results_pwy <- pwy_long %>%
  group_by(PWY) %>%
  group_modify(~{
    tryCatch({
      mod <- lm(Abundance ~ Treatment + Day_numeric, data = .x)
      tidy_mod <- broom::tidy(mod)
      treat_row <- tidy_mod %>% filter(term == "TreatmentIA")
      tibble(
        Estimate = treat_row$estimate,
        StdError = treat_row$std.error,
        p_value = treat_row$p.value
      )
    }, error = function(e){
      tibble(Estimate = NA, StdError = NA, p_value = NA)
    })
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

results_pwy<- results_pwy %>% arrange(desc(abs(Estimate)))
head(results_pwy)

#' save as csv file 
write.csv(results_pwy, file= "linearmodel.csv")

#================================================================================================
                    # map pathway id #
library(stringr)

#' read the CSV
map_pathway <- read.csv("Pathway_id_map.csv", stringsAsFactors = FALSE)

#' fix column names (remove trailing space in "Name ")
colnames(map_pathway) <- trimws(colnames(map_pathway))

#' remove HTML tags like <i>S</i>
map_pathway$Name <- str_replace_all(map_pathway$Name, "<.*?>", "")

#' inspect
head(map_pathway)

#' read linear model output 
linearmodel<- read.csv("linearmodel.csv", stringsAsFactors = TRUE)

#' merge them 
#' Clean column names first
colnames(map_pathway) <- trimws(colnames(map_pathway))
colnames(linearmodel) <- trimws(colnames(linearmodel))

#' Clean IDs (important!)
map_pathway$ID <- trimws(map_pathway$ID)
linearmodel$ID       <- trimws(linearmodel$ID)

# Merge: keep the pathway name from map_pathway
merged_df <- merge(
  linearmodel, 
  map_pathway[, c("ID", "Name")],
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)
merged_df

library(dplyr)
merged_df<- as.data.frame(merged_df)
sig_paths <- merged_df[merged_df$p_value < 0.05, c("ID", "Name", "Estimate", "p_adj")]


#' Write csv file 
write.csv(merged_df, file= "metabolic.pathway.IAvspreici.csv")

top_positive <- merged_df %>%
  arrange(desc(Estimate)) %>%
  slice(1:5)

top_negative <- merged_df %>%
  arrange(Estimate) %>%   # ascending = most negative first
  slice(1:5)

top10 <- rbind(top_positive, top_negative)

#' Mean standardize difference 
library(tidyverse)

#' Read pathway abundance table (PWY pathways)
path_abun <- read_tsv("picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")

#' Suppose you have a metadata table with group info
metadata <- read.csv("metadata_arthritis_second_batch.csv", stringsAsFactors = TRUE)

#' Join abundance table with metadata
path_long <- path_abun %>%
  pivot_longer(-pathway, names_to="Sample", values_to="Abundance") %>%
  left_join(metadata, by="Sample")

#'  Filter only the two groups of interest
path_long <- path_long %>%
  filter(Treatment %in% c("Pre-ICI", "IA"))

#' Calculate SMD (IA vs Pre-ICI)
smd_table <- path_long %>%
  group_by(pathway) %>%
  summarize(
    mean_IA     = mean(Abundance[Treatment == "IA"]),
    mean_PreICI = mean(Abundance[Treatment == "Pre-ICI"]),
    sd_IA       = sd(Abundance[Treatment == "IA"]),
    sd_PreICI   = sd(Abundance[Treatment == "Pre-ICI"]),
    n_IA        = sum(Treatment == "IA"),
    n_PreICI    = sum(Treatment == "Pre-ICI")
  ) %>%
  mutate(
    pooled_sd = sqrt(((n_IA - 1)*sd_IA^2 + (n_PreICI - 1)*sd_PreICI^2) /
                     (n_IA + n_PreICI - 2)),
    SMD = (mean_IA - mean_PreICI) / pooled_sd
  )

pval_table <- path_long %>%
  group_by(pathway) %>%
  summarize(
    p_value = wilcox.test(Abundance ~ Treatment)$p.value
  )

#' Combine
final_table <- smd_table %>%
  left_join(pval_table, by = "pathway")

#' next is to map with pathway annotation !! 

                    ## map pathway id 
library(stringr)

#' read the CSV
map_pathway <- read.csv("Pathway_id_map.csv", stringsAsFactors = FALSE)

#' fix column names (remove trailing space in "Name ")
colnames(map_pathway) <- trimws(colnames(map_pathway))

#' remove HTML tags like <i>S</i>
map_pathway$Name <- str_replace_all(map_pathway$Name, "<.*?>", "")

#' inspect
head(map_pathway)

#' merge them 
#' Clean column names first
colnames(map_pathway) <- trimws(colnames(map_pathway))
colnames(final_table) <- trimws(colnames(final_table))

#' Clean IDs (important!)
map_pathway$pathway <- trimws(map_pathway$pathway)
final_table$pathway     <- trimws(final_table$pathway)

#' Merge: keep the pathway name from map_pathway
merged_df <- merge(
  final_table, 
  map_pathway[, c("pathway", "Name")],
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)
merged_df

sig_df <- merged_df %>%
  filter(p_value < 0.05)

merged_df$padj <- p.adjust(merged_df$p_value, method = "fdr")

#' Filter top-5 
top5_positive <- sig_df %>%
  arrange(desc(SMD)) %>%      # highest SMD first
  slice(1:6)

top5_negative <- sig_df %>%
  arrange(SMD) %>%            # most negative first
  slice(1:6)

top12 <- bind_rows(
  top5_positive %>% mutate(direction = "Higher in IA"),
  top5_negative %>% mutate(direction = "Higher in Pre-ICI")
)

#' Add p_sig to top12
top12_with_p <- top12 %>%
  mutate(p_sig = case_when(
    p_value < 0.01 ~ "< 0.01",
    p_value < 0.05 ~ "< 0.05",
    TRUE ~ ">= 0.05"
  ))

#' Plot
library(tidyverse)

metabolic_pathway<- ggplot(top12, aes(x = SMD, y = reorder(Name, SMD), fill = p_value)) +
  geom_col() +
  theme_void() +
  xlab("Standardized Mean Difference (IA vs Pre-ICI)") +
  ylab("Pathway") +
  scale_fill_gradient(
    name = "p-value",
    low = "red",  # red for low p-values
    high = "blue"   # grey for high p-values
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(), 
    legend.position = "bottom")  
metabolic_pathway

ggsave("metabolic_pathway.tiff", plot = metabolic_pathway, dpi = 300, width = 4, height = 2.6, units = "in")

#' Predict metabolite from corresponding pathways 
library(tidyverse)
library(ggalluvial)

df_alluvial <- metabolite %>%
  mutate(
    Order = ifelse(Order == "Top", 
                   "Top (ICI-IA)", 
                   "Bottom (Pre-ICI)")
  )

df_alluvial <- df_alluvial %>%
  mutate(
    is_alluvia_form = TRUE
  )

ggplot(df_alluvial,
       aes(axis1 = Group, axis2 = Order, axis3 = Metabolite,
           y = 1)) +
  geom_alluvium(aes(fill = Group), width = 0.2, alpha = 0.8) +
  geom_stratum(width = 0.2, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("ICI-IA" = "#1f78b4",
                               "Pre-ICI" = "#e31a1c")) +
  labs(title = "Metabolite Flow: Group → Order → Metabolite",
       x = "", y = "") +
  theme_minimal(base_size = 14)


#' Alluvial 
library(tidyverse)
library(ggalluvial)

 Load your data
df <- read.csv("Metabolite_ICI_IAvsPre_ICI.csv", stringsAsFactors = FALSE)

#' Prepare data for alluvial
df_alluvial <- df %>%
  mutate(is_alluvia_form = TRUE)

#' Plot
pathway_metabolite<- library(tidyverse)
library(ggalluvial)

df <- read.csv("Metabolite_ICI_IAvsPre_ICI.csv", stringsAsFactors = FALSE)

df_alluvial <- df %>%
  mutate(is_alluvia_form = TRUE)

Alluvial.metabolite<- ggplot(df_alluvial,
       aes(axis1 = Pathway,
           axis2 = Metabolite,
           axis3 = Group,
           y = 0.5)) +

  # ---- STRAIGHT LINES ----
  geom_flow(aes(fill = Group),
            stat = "alluvium",
            lwd = 0.5,
            alpha = 0.9,
            width = 0.25) +

  geom_stratum(width = 0.45,
               fill = "grey95",
               color = "black") +

  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3, 
            fontface= "bold") +

  scale_fill_manual(values = c(
    "ICI-IA" = "#1f78b4",
    "Pre-ICI" = "#e31a1c"
  )) +

  scale_x_discrete(
    limits = c("Pathway", "Metabolite", "Group"),
    expand = c(.05, .05)
  ) +

  labs(
    title = "Pathways → Predicted Metabolites",
    x = "",
    y = ""
  ) +

  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold", hjust=0.5)
  )
Alluvial.metabolite

#' save high quality plot 
ggsave("metabolite_ICIIAvspreICI.tiff", plot = Alluvial.metabolite, dpi = 300, width = 6.3, height = 3, units = "in")

#' Mean difference adjusting for time!! 
library(tidyverse)
library(tidyverse)
library(emmeans)
library(purrr)
library(broom)

#' Read pathway abundance table (PWY pathways)
path_abun <- read_tsv("picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")

#' Suppose you have a metadata table with group info
metadata <- read.csv("metadata_arthritis_second_batch.csv", stringsAsFactors = TRUE)

#' Join abundance table with metadata
path_long <- path_abun %>%
  pivot_longer(-pathway, names_to="Sample", values_to="Abundance") %>%
  left_join(metadata, by="Sample")

#' 4. Filter only the two groups of interest
path_long <- path_long %>%
  filter(Treatment %in% c("Pre-ICI", "IA"))

#' Make sure Treatment levels are consistent
#' path_long <- path_long %>%
 # mutate(Treatment = factor(Treatment, levels = c("Pre-ICI", "IA"))) # change levels


#' Function to fit model and extract metrics for one pathway
analyze_pathway <- function(df_path){
  tryCatch({
    if (is.numeric(df_path$Day_factor)) {
      mod <- lm(Abundance ~ Treatment + Day_factor, data = df_path)
    } else {
      mod <- lm(Abundance ~ Treatment + factor(Day_factor), data = df_path)
    }

    emm <- emmeans(mod, ~Treatment)
    contr <- contrast(emm, method = "pairwise")
    contr_df <- as.data.frame(contr)
    emm_df <- as.data.frame(emm)

    # Adjusted mean difference (IA - Pre-ICI)
    est <- contr_df$estimate[1]
    pval <- contr_df$p.value[1]

    resid_sd <- sigma(mod)

    mu_pre <- emm_df$emmean[emm_df$Treatment == "Pre-ICI"]
    mu_ia  <- emm_df$emmean[emm_df$Treatment == "IA"]

    SMD_adj <- (mu_ia - mu_pre) / resid_sd

    df_resid <- df.residual(mod)
    J <- ifelse(df_resid > 2, 1 - (3 / (4 * df_resid - 1)), NA_real_)
    hedges_g <- SMD_adj * J

    tibble(
      estimate = est,
      SMD_adj = SMD_adj,
      hedges_g = hedges_g,
      p_value = pval,
      n = nrow(df_path)
    )
  }, error = function(e){
    tibble(
      estimate = NA_real_,
      SMD_adj = NA_real_,
      hedges_g = NA_real_,
      p_value = NA_real_,
      n = nrow(df_path)
    )
  })
}

results <- path_long %>%
  group_by(pathway) %>%
  group_split() %>%
  purrr::map_df(analyze_pathway)


#' apply per pathway
results <- path_long %>%
  group_by(pathway) %>%
  group_map(~ analyze_pathway(.x), .keep = TRUE) %>%
  bind_rows(.id = "pathway") %>%
  # if group_map returned integer names, ensure Pathway column is proper
  mutate(pathway = unique(path_long$pathway)[as.integer(pathway)]) %>%
  select(pathway, everything())

results %>% arrange(p_value) %>% slice(1:50)

#' Map the pathways!! 
library(stringr)

#' read the CSV
map_pathway <- read.csv("Pathway_id_map.csv", stringsAsFactors = FALSE)

#' fix column names (remove trailing space in "Name ")
colnames(map_pathway) <- trimws(colnames(map_pathway))

#' remove HTML tags like <i>S</i>
map_pathway$Name <- str_replace_all(map_pathway$Name, "<.*?>", "")

#' inspect
head(map_pathway)

#' merge them 
#' Clean column names first
colnames(map_pathway) <- trimws(colnames(map_pathway))
colnames(results) <- trimws(colnames(results))

#' Clean IDs (important!)
map_pathway$pathway <- trimws(map_pathway$pathway)
results$pathway     <- trimws(results$pathway)

#' Merge: keep the pathway name from map_pathway
merged_df <- merge(
  results, 
  map_pathway[, c("pathway", "Name")],
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)
merged_df
#' Merge: keep the pathway name from map_pathway
merged_df <- merge(
  results, 
  map_pathway[, c("pathway", "Name")],
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)
merged_df

sig_df <- merged_df %>%
  filter(p_value < 0.05)

merged_df$padj <- p.adjust(merged_df$p_value, method = "fdr")

#' Filter top-5 
top5_positive <- sig_df %>%
  arrange(desc(SMD_adj)) %>%      # highest SMD first
  slice(1:6)

top5_negative <- sig_df %>%
  arrange(SMD_adj) %>%            # most negative first
  slice(1:6)

top12 <- bind_rows(
  top5_positive %>% mutate(direction = "Higher in IA"),
  top5_negative %>% mutate(direction = "Higher in Pre-ICI")
)

library(tidyverse)

#' Add p_sig to top12
top12_with_p <- top12 %>%
  mutate(p_sig = case_when(
    p_value < 0.01 ~ "< 0.01",
    p_value < 0.05 ~ "< 0.05",
    TRUE ~ ">= 0.05"
  ))

#' Plot
library(tidyverse)

metabolic_pathway<- ggplot(top12, aes(x = SMD, y = reorder(Name, SMD), fill = p_value)) +
  geom_col() +
  theme_void() +
  xlab("Standardized Mean Difference (IA vs Pre-ICI)") +
  ylab("Pathway") +
  scale_fill_gradient(
    name = "p-value",
    low = "red",  # red for low p-values
    high = "blue"   # grey for high p-values
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(), 
    legend.position = "bottom")  
metabolic_pathway

ggsave("metabolic_pathway.tiff", plot = metabolic_pathway, dpi = 300, width = 4, height = 2.6, units = "in")

#' Pathway to metabolites 

#' Plot
pathway_metabolite<- library(tidyverse)
library(ggalluvial)

df <- read.csv("Metabolite_ICI_IAvsPre_ICI.csv", stringsAsFactors = FALSE)

df_alluvial <- df %>%
  mutate(is_alluvia_form = TRUE)

Alluvial.metabolite<- ggplot(df_alluvial,
       aes(axis1 = Pathway,
           axis2 = Metabolite,
           axis3 = Group,
           y = 0.5)) +

  # ---- STRAIGHT LINES ----
  geom_flow(aes(fill = Group),
            stat = "alluvium",
            lwd = 0.5,
            alpha = 0.9,
            width = 0.25) +

  geom_stratum(width = 0.45,
               fill = "grey95",
               color = "black") +

  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3, 
            fontface= "bold") +

  scale_fill_manual(values = c(
    "ICI-IA" = "#1f78b4",
    "Pre-ICI" = "#e31a1c"
  )) +

  scale_x_discrete(
    limits = c("Pathway", "Metabolite", "Group"),
    expand = c(.05, .05)
  ) +

  labs(
    title = "Metabolites",
    x = "",
    y = ""
  ) +

  theme_void() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 13, face = "bold", hjust=0.5)
  )
Alluvial.metabolite

#' Mix model !! 
#' Run linear mixed-effects models per pathway  #did not work

# Formula: Abundance ~ Treatment + Day + (1|Sample)
results_pwy <- pwy_long %>%
  group_by(PWY) %>%
  group_modify(~{
    tryCatch({
      mod <- lmer(Abundance ~ Treatment + Day_numeric + (1|Sample), data = .x)
      tidy_mod <- broom.mixed::tidy(mod, effects = "fixed")
      
      # Extract Treatment coefficient
      treat_row <- tidy_mod %>% filter(term == "TreatmentIA")  # adjust if your factor levels differ
      
      tibble(
        Estimate = treat_row$estimate,
        StdError = treat_row$std.error,
        p_value = treat_row$p.value
      )
    }, error = function(e){
      tibble(Estimate = NA, StdError = NA, p_value = NA)
    })
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"))

#' Optional: filter significant pathways

sig_pwy <- results_pwy %>%
  filter(!is.na(p_adj)) %>%
  arrange(p_adj)

#' 6. Plot bar figure (Figure 1 style)
#' Merge coefficient with PWY names for plotting
plot_df <- results_pwy %>% filter(!is.na(Estimate))

ggplot(plot_df, aes(x = Estimate, y = reorder(PWY, Estimate))) +
  geom_col(aes(fill = p_adj < 0.05)) +
  scale_fill_manual(values = c("grey", "red")) +
  theme_minimal() +
  xlab("Treatment effect (Coefficient, IA vs Pre-ICI)") +
  ylab("Pathway") +
  theme(legend.position = "none") +
  ggtitle("Differential PWY pathway analysis") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

#' Check for NA in Treatment
sum(is.na(pwy_long$Treatment))

# Check how many samples per pathway
pwy_long %>%
  group_by(PWY, Treatment) %>%
  summarise(n = n())

#' Run KEGG 
# Load libraries
library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)  # placeholder, not used for KEGG KO

#' Read KO abundances per sample
ko_abun <- read_tsv("picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
colnames(ko_abun)[1] <- "KO"

#' 2. Pivot longer

ko_long <- ko_abun %>%
  pivot_longer(-KO, names_to = "Sample", values_to = "Abundance")

#' Merge with metadata

metadata <- read.csv("metadata_arthritis_second_batch.csv", stringsAsFactors = TRUE) %>%
  filter(!Treatment %in% c("HC", "exclusion"))

ko_long <- ko_long %>%
  filter(Sample %in% metadata$Sample) %>%
  left_join(metadata, by = "Sample")

#' Aggregate per pathway (optional, can keep KO-level)
# Here we can keep KO-level abundances for DE analysis
library(dplyr)
library(tidyr)
library(tibble)

ko_path <- ko_long %>%
  dplyr::select(Sample, KO, Abundance) %>%   # explicitly use dplyr::select
  pivot_wider(names_from = Sample, values_from = Abundance) %>%
  column_to_rownames("KO")

#' DESeq2 differential KOs (Treatment vs Pre-ICI, adjusting for Day)
metadata_deseq <- metadata %>%
  column_to_rownames("Sample")

dds <- DESeqDataSetFromMatrix(
  countData = round(ko_path),  # DESeq2 expects counts
  colData = metadata_deseq,
  design = ~ split
)

dds <- DESeq(dds)
res_deseq <- results(dds, contrast = c("split", "ICI_IA_D29", "Pre_ICI"))
head(res_deseq)


#' KEGG enrichment with clusterProfiler

all_kos <- rownames(res_deseq)[!is.na(res_deseq$log2FoldChange)]
# Ensure format is clean
all_kos <- toupper(str_trim(all_kos))

# Run KEGG enrichment without filtering
kk <- enrichKEGG(
  gene = res_deseq,
  organism = "ko",
  keyType = "kegg",
  universe = res_deseq,
  pvalueCutoff = 1,   # set high so no pathways are filtered out
  qvalueCutoff = 1    # same for FDR
)

# Inspect results
write.csv(kk, file = "EnrichKEGG.csv")

#' Visualize enrichment

library(enrichplot)

barplot(kk_up, showCategory = 20, title = "KEGG Pathways Up in Treatment")
dotplot(kk_down, showCategory = 20, title = "KEGG Pathways Down in Treatment")


#' Inspect results tables
head(as.data.frame(kk_up))
head(as.data.frame(kk_down))

#' Test KO different 
library(data.table)
ko_per_sample <- fread("picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", data.table=FALSE)
rownames(ko_per_sample) <- ko_per_sample$`function`
ko_per_sample$`function` <- NULL
ko_per_sample$

#check table 
dim(ko_per_sample)         # number of KOs × number of samples
colnames(ko_per_sample)    # sample IDs (use to link to metadata)
head(ko_per_sample[,1:5])  # preview first 5 samples

#' metadata 
metadata <- read.csv("metadata_arthritis_second_batch.csv", stringsAsFactors = TRUE) 
ko_per_sample <- ko_per_sample[, metadata$Sample]


#' Convert KO table to long format for mixed modeling
library(tidyr)
library(dplyr)

#' Transpose KO table: samples as rows, KOs as columns
ko_long <- as.data.frame(t(ko_per_sample))
ko_long$Sample <- rownames(ko_long)

#' Merge metadata
ko_long <- left_join(ko_long, metadata, by = c("Sample" = "Sample"))

# Pivot longer: columns = KOs, rows = Sample x KO
ko_long <- pivot_longer(ko_long, 
                        cols = starts_with("K"), 
                        names_to = "KO", 
                        values_to = "abundance")


ko_long <- ko_long %>% filter(!Treatment %in% c("exclusion", "HC"))
write.csv(ko_long, file= "ko_long.csv")

#' fit linear mixed model 
library(lme4)
library(lme4)
library(dplyr)

#' Ensure 'group' is a factor with reference
ko_long$Treatment <- factor(ko_long$Treatment, levels = c("Pre-ICI", "IA"))  # adjust as needed

#' Fit mixed model per KO
results <- ko_long %>%
  group_by(KO) %>%
  do({
    model <- lm(abundance ~ Treatment + Day_numeric, data = .)
    coef_summary <- summary(model)$coefficients
    if("TreatmentIA" %in% rownames(coef_summary)) {
      data.frame(
        estimate = coef_summary["TreatmentIA","Estimate"],
        pvalue  = coef_summary["TreatmentIA","Pr(>|t|)"]
      )
    } else {
      data.frame(
        estimate = NA,
        pvalue = NA
      )
    }
  }) %>% ungroup()

results$FDR <- p.adjust(results$pvalue, method = "BH")
sig_ko <- results %>% filter(!is.na(FDR) & FDR < 0.05)

ko_per_mouse <- ko_long %>%
  group_by(real_id, Treatment, KO) %>%
  summarise(mean_abundance = mean(abundance), .groups="drop")

#' Then run paired t-test per KO
results <- ko_per_mouse %>%
  group_by(KO) %>%
  summarise(
    t_test = t.test(mean_abundance ~ Treatment)$p.value
  )
results$FDR <- p.adjust(results$t_test, method="BH")

#' Keep only significant hits (FDR < 0.05)
sig_results <- results[results$FDR < 0.05, ]

# Session info
sessionInfo()
