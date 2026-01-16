##############################################################
# 16S rRNA Microbial Community Analysis using Phyloseq & Qiime2R
# Author: Dr. Sobia Ajaz
# Updated: October 2025
##############################################################

# Load required packages
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(rstatix)
library(ggpubr)
library(RColorBrewer)

##############################################################
# 1️⃣ Data Import from QIIME2
##############################################################

# Import feature table (ASVs/OTUs)
otu_table <- read_qza("table.qza")$data

# Import taxonomy
taxonomy <- read_qza("taxonomy.qza")$data

# Import rooted phylogenetic tree (if available)
phy_tree <- read_qza("rooted-tree.qza")$data

# Import sample metadata (edit to your filename)
metadata <- read_tsv("sample-metadata.tsv")

# Convert data to phyloseq components
OTU <- otu_table(otu_table, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(read_qza("taxonomy.qza")$data))
META <- sample_data(metadata)
TREE <- phy_tree(phy_tree)

# Combine into a single phyloseq object
physeq <- phyloseq(OTU, TAX, META, TREE)

##############################################################
# 2️⃣ Basic Filtering and Cleanup
##############################################################

# Remove unassigned or unwanted taxa
physeq <- subset_taxa(physeq, !is.na(Family) & Family != "")

# Remove chloroplasts and mitochondria
physeq <- subset_taxa(physeq, !Family %in% c("Chloroplast", "Mitochondria"))

# Optionally remove low-abundance taxa
physeq <- prune_taxa(taxa_sums(physeq) > 10, physeq)

##############################################################
# 3️⃣ Rarefaction / Normalization
##############################################################

# Rarefy to even sequencing depth (for diversity analysis)
set.seed(123)
min_depth <- min(sample_sums(physeq))
physeq_rarefied <- rarefy_even_depth(physeq, sample.size = min_depth, rngseed = 123, replace = FALSE)

##############################################################
# 4️⃣ Taxonomic Agglomeration (Family Level)
##############################################################

family_physeq <- tax_glom(physeq_rarefied, taxrank = "Family")

# Transform to relative abundance
family_rel <- transform_sample_counts(family_physeq, function(x) x / sum(x))

##############################################################
# 5️⃣ Stacked Bar Plot (Relative Abundance)
##############################################################

family_df <- psmelt(family_rel)
family_df$Family <- as.character(family_df$Family)

# Optional: select top 15 families by mean abundance
top_families <- family_df %>%
  group_by(Family) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  top_n(15, mean_abundance) %>%
  pull(Family)

family_df$Family_grouped <- ifelse(family_df$Family %in% top_families, family_df$Family, "Other")

# Plot
p_bar <- ggplot(family_df, aes(x = Sample, y = Abundance, fill = Family_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(length(unique(family_df$Family_grouped)))) +
  theme_bw(base_size = 14) +
  labs(title = "Relative Abundance of Top Bacterial Families",
       x = "Sample", y = "Relative Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

print(p_bar)

##############################################################
# 6️⃣ Alpha Diversity (Shannon Index)
##############################################################

alpha_div <- estimate_richness(physeq_rarefied, measures = "Shannon")
alpha_div$SampleID <- rownames(alpha_div)
alpha_div <- merge(alpha_div, metadata, by.x = "SampleID", by.y = "SampleID")

# Statistical comparison
shannon_test <- alpha_div %>% wilcox_test(Shannon ~ Treatment)
print(shannon_test)

# Boxplot
p_alpha <- ggboxplot(alpha_div, x = "Treatment", y = "Shannon",
                     fill = "Treatment", palette = "Set2") +
  stat_compare_means(comparisons = list(c("Control", "Treatment")),
                     method = "wilcox.test") +
  theme_bw(base_size = 14) +
  labs(title = "Alpha Diversity (Shannon Index)",
       x = "Treatment", y = "Shannon Index")

print(p_alpha)

##############################################################
# 7️⃣ Beta Diversity (Bray–Curtis + PCoA)
##############################################################

# Compute Bray–Curtis dissimilarity
beta_dist <- phyloseq::distance(family_physeq, method = "bray")

# Ordination (PCoA)
ordination <- ordinate(family_physeq, method = "PCoA", distance = beta_dist)

# Prepare ordination dataframe
ordination_scores <- as.data.frame(ordination$vectors)
ordination_scores$SampleID <- rownames(ordination_scores)
ordination_scores <- left_join(ordination_scores, metadata, by = "SampleID")

# Plot ordination
p_beta <- plot_ordination(family_physeq, ordination, color = "Plant", shape = "Clamp") +
  geom_point(aes(color = Plant, shape = Clamp), size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = interaction(Plant, Treatment_Sequencing),
                   linetype = Treatment_Sequencing),
               type = "t", size = 1) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
  theme_bw(base_size = 16) +
  labs(
    title = "Bray–Curtis PCoA",
    x = paste0("PCoA1 (", round(ordination$values$Relative_eig[1] * 100, 1), "%)"),
    y = paste0("PCoA2 (", round(ordination$values$Relative_eig[2] * 100, 1), "%)")
  ) +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

print(p_beta)

##############################################################
# 8️⃣ PERMANOVA (Statistical test on beta diversity)
##############################################################

adonis_result <- adonis2(beta_dist ~ Plant * Clamp * Treatment_Sequencing, data = metadata)
print(adonis_result)

##############################################################
# 9️⃣ Save All Objects for Future Use
##############################################################

save(
  physeq, physeq_rarefied, family_physeq, family_rel, family_df,
  alpha_div, beta_dist, ordination, adonis_result,
  file = "16S_phyloseq_analysis.RData"
)

##############################################################
# ✅ End of Workflow
##############################################################
