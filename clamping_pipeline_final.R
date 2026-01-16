##############################
# Microbiome Analysis Pipeline
# QIIME2 → Phyloseq → Diversity → PICRUSt2 → Pathways
##############################

# ==== 1. Load Libraries ====
library(phyloseq)
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggplot2)
library(microbiome)
library(reshape2)
library(phyloseq)
library(qiime2R)
library(readxl)
library(rstudioapi)
library(tidyr)
library(ggplot2)
library(maps)
library(sf)
library(nlme)
library(plyr)
library(microViz)
library(vegan)  # For diversity calculations
library(ggplot2)  # For visualization
library(dplyr)  # For data manipulation
library(rstatix)  # For statistical tests
library(gridExtra)  # For arranging multiple plots
library(biomformat)
library(Biostrings)
library(pheatmap)
library(KEGGREST)
library(reshape2)
library(httr)
library(jsonlite)

ASVs <- read_qza("feature-table-16S-Plant-all.qza")      # ASV table
metadata <- read_q2metadata("metadata_Plant.txt")        # sample metadata
taxonomy <- read_qza("taxonomy_SILVA-16S-Plant-all.qza") # taxonomy
taxonomy <- parse_taxonomy(taxonomy$data)

physeq <- qza_to_phyloseq(
  features = "feature-table-16S-Plant-all.qza",
  tree = "rooted-tree-16S-Plant-all.qza",
  taxonomy = "taxonomy_SILVA-16S-Plant-all.qza",
  metadata = "metadata_Plant.txt"
)


# ==== 2. Import QIIME2 Data ====
asv_tab <- read_qza("feature-table-16S-Plant-all.qza")
tax_tab <- read_qza("taxonomy_SILVA-16S-Plant-all.qza")
tree <- read_qza("rooted-tree-16S-Plant-all.qza")
meta_data <- read_tsv("metadata.tsv") %>% as.data.frame()
rownames(meta_data) <- meta_data$`#SampleID`

# ==== 3. Build Phyloseq Object ====
otu <- otu_table(asv_tab$data, taxa_are_rows=TRUE)
samp <- sample_data(meta_data)
tax <- tax_tab$data %>%
  separate(Taxon, into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep="; ", fill="right") %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()
tree <- phy_tree(tree$data)

ps <- phyloseq(otu, samp, tax_table(tax), tree)

# ==== 4. Rarefy Data ====
ps.rarefied <- rarefy_even_depth(ps, rngseed=123, sample.size=min(sample_sums(ps)), replace=FALSE)

# ==== 5. Taxa Aggregation ====
ps.family <- tax_glom(ps.rarefied, taxrank="Family")
ps.family.rel <- transform_sample_counts(ps.family, function(x) x/sum(x))

# ==== 6. Plotting Functions ====
plot_barplot <- function(ps, topN=20, rel=TRUE, rm_chloro=TRUE) {
  df <- psmelt(ps)
  if (rm_chloro) {
    df <- df %>% filter(!Family %in% c("Mitochondria", "Chloroplast"))
  }
  top_families <- df %>% group_by(Family) %>%
    summarise(Abundance=sum(Abundance)) %>%
    top_n(topN, Abundance) %>% pull(Family)
  df <- df %>% filter(Family %in% top_families)
  
  ggplot(df, aes(x=Sample, y=Abundance, fill=Family)) +
    geom_bar(stat="identity", position=ifelse(rel, "fill", "stack")) +
    facet_wrap(~Treatment, scales="free_x") +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
}

# Example barplots
plot_barplot(ps.family, topN=20, rel=FALSE)
plot_barplot(ps.family.rel, topN=20, rel=TRUE)

# ==== 7. Alpha Diversity ====
alpha_df <- estimate_richness(ps.rarefied, measures=c("Observed","Shannon")) %>%
  rownames_to_column("Sample") %>%
  left_join(meta_data, by=c("Sample"="#SampleID")) %>%
  pivot_longer(cols=c("Observed","Shannon"), names_to="Index", values_to="value")

# Test function (Wilcoxon if 2 groups, Kruskal if >2)
test_alpha <- function(df, index="Observed") {
  d <- df %>% filter(Index==index)
  if (n_distinct(d$Treatment) > 2) {
    kruskal.test(value ~ Treatment, data=d)
  } else {
    wilcox.test(value ~ Treatment, data=d)
  }
}

# Example tests
test_alpha(alpha_df, "Observed")
test_alpha(alpha_df, "Shannon")

# Plot alpha diversity
alpha_df %>% ggplot(aes(x=Treatment, y=value, fill=Treatment)) +
  geom_boxplot() + facet_wrap(~Index, scales="free_y") + theme_bw()

# ==== 8. Beta Diversity ====
dist <- phyloseq::distance(ps.rarefied, method="bray")
ord <- ordinate(ps.rarefied, method="PCoA", distance=dist)

plot_ordination(ps.rarefied, ord, color="Treatment") +
  geom_point(size=3) +
  stat_ellipse(type="t", level=0.95) +
  theme_bw()

# PERMANOVA test
adonis2(dist ~ Treatment, data=meta_data)

# ==== 9. Export for PICRUSt2 ====
write.table(as.data.frame(otu_table(ps)), "ASV_table.tsv", sep="\t", quote=FALSE, col.names=NA)
write.table(as.data.frame(tax_table(ps)), "ASV_taxonomy.tsv", sep="\t", quote=FALSE, col.names=NA)

# ==== 10. Import PICRUSt2 Pathways ====
path_abun <- read.delim("path_abun_unstrat.tsv", skip=1, row.names=1, check.names=FALSE)
path_desc <- read.delim("path_abun_unstrat_descrip.tsv")

# Merge & reshape
top_pathways <- path_abun %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to="Sample", values_to="Abundance") %>%
  left_join(path_desc, by="Pathway") %>%
  left_join(meta_data %>% rownames_to_column("Sample"), by="Sample")

# Select top 30
top30 <- top_pathways %>% group_by(Pathway) %>%
  summarise(Total=sum(Abundance)) %>% top_n(30, Total) %>% pull(Pathway)

heatmap_df <- top_pathways %>% filter(Pathway %in% top30) %>%
  dcast(Pathway ~ Sample, value.var="Abundance") %>%
  column_to_rownames("Pathway")

pheatmap::pheatmap(log1p(heatmap_df), scale="row", clustering_distance_rows="euclidean")
