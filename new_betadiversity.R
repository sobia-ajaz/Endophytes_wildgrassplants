library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(compositions) 
library(rstudioapi)
# For CLR transformation

# Set working directory
setwd(dirname(getActiveDocumentContext()$path))

# Load and prepare data
physeq <- qza_to_phyloseq(
  features = "feature-table-16S-Plant-all.qza",
  tree = "rooted-tree-16S-Plant-all.qza",
  taxonomy = "taxonomy_SILVA-16S-Plant-all.qza",
  metadata = "metadata_Plant.txt"
)

# Rarefy
set.seed(1)
S16_scaled <- rarefy_even_depth(physeq, sample.size = 40000, replace = FALSE, rngseed = 1)

# Subset to plants
plant_subset_physeq <- subset_samples(S16_scaled, Plant %in% c("Buttercup", "Holcus", "White clover"))

# Remove chloroplast and mitochondria
clean_physeq <- subset_taxa(plant_subset_physeq, 
                            !grepl("chloroplast|Chloroplast", Order, ignore.case = TRUE) &
                              !grepl("chloroplast|Chloroplast", Family, ignore.case = TRUE) &
                              !grepl("mitochondria|Mitochondria", Family, ignore.case = TRUE))

print(paste("Final clean data:", ntaxa(clean_physeq), "ASVs after removing chloroplast/mitochondria"))

# Create simplified treatment groups (3 main types)
sample_data(clean_physeq)$TreatmentSimple <- case_when(
  sample_data(clean_physeq)$Treatment == "Sterilized" & sample_data(clean_physeq)$sequencing == "Optimized" ~ "Sterilized Optimized",
  sample_data(clean_physeq)$Treatment == "Sterilized" & sample_data(clean_physeq)$sequencing == "Standard" ~ "Sterilized Standard", 
  sample_data(clean_physeq)$Treatment == "Washed" ~ "Washed",
  TRUE ~ "Other"
)

# Remove any "Other" category
clean_physeq <- subset_samples(clean_physeq, TreatmentSimple != "Other")

# Create separate variables for plotting
sample_data(clean_physeq)$ClampStatus <- ifelse(sample_data(clean_physeq)$Clamp == "clamp", "Clamp", "No Clamp")

# Set factor orders
sample_data(clean_physeq)$TreatmentSimple <- factor(sample_data(clean_physeq)$TreatmentSimple,
                                                    levels = c("Sterilized Standard", "Sterilized Optimized", "Washed"))
sample_data(clean_physeq)$Plant <- factor(sample_data(clean_physeq)$Plant,
                                          levels = c("Buttercup", "Holcus", "White clover"))
sample_data(clean_physeq)$ClampStatus <- factor(sample_data(clean_physeq)$ClampStatus,
                                                levels = c("Clamp", "No Clamp"))

# =============================================================================
# CENTER LOG RATIO (CLR) TRANSFORMATION
# =============================================================================

print("=== CENTER LOG RATIO (CLR) TRANSFORMATION ===")

# Extract OTU table
otu_table <- as(otu_table(clean_physeq), "matrix")
if (!taxa_are_rows(clean_physeq)) {
  otu_table <- t(otu_table)
}

print(paste("OTU table dimensions:", nrow(otu_table), "taxa x", ncol(otu_table), "samples"))

# Add pseudocount to handle zeros (required for CLR)
pseudocount <- min(otu_table[otu_table > 0]) / 2
otu_table_clr <- otu_table + pseudocount

print(paste("Added pseudocount:", pseudocount))

# Perform CLR transformation
clr_transformed <- apply(otu_table_clr, 2, function(x) {
  log(x) - mean(log(x))
})

# Convert back to matrix and transpose if needed
clr_matrix <- t(clr_transformed)

print("CLR transformation completed")
print(paste("CLR matrix dimensions:", nrow(clr_matrix), "samples x", ncol(clr_matrix), "taxa"))

# Create new phyloseq object with CLR-transformed data
original_taxa_are_rows <- taxa_are_rows(clean_physeq)
print(paste("Original taxa_are_rows:", original_taxa_are_rows))

if (original_taxa_are_rows) {
  clr_otu <- otu_table(t(clr_matrix), taxa_are_rows = TRUE)
} else {
  clr_otu <- otu_table(clr_matrix, taxa_are_rows = FALSE)
}

clr_physeq <- phyloseq(
  otu_table = clr_otu,
  tax_table = tax_table(clean_physeq),
  sample_data = sample_data(clean_physeq)
)

print("Created new phyloseq object with CLR-transformed data")

# =============================================================================
# BETA DIVERSITY WITH CLR-TRANSFORMED DATA
# =============================================================================

print("=== BETA DIVERSITY WITH CLR-TRANSFORMED DATA ===")

# Calculate Euclidean distance on CLR-transformed data
clr_euclidean_dist <- phyloseq::distance(clr_physeq, method = "euclidean")
print("Euclidean distance matrix calculated from CLR-transformed data")

# PCoA ordination on CLR data
pcoa_clr <- ordinate(clr_physeq, method = "PCoA", distance = "euclidean")

# Calculate variance explained
variance_explained_clr <- round(100 * pcoa_clr$values$Eigenvalues / sum(pcoa_clr$values$Eigenvalues), 2)

print(paste("Variance explained (CLR) - PCoA1:", variance_explained_clr[1], "%, PCoA2:", variance_explained_clr[2], "%"))

# Extract ordination scores for custom plotting
ordination_scores_clr <- as.data.frame(pcoa_clr$vectors[, 1:2])
colnames(ordination_scores_clr) <- c("PCoA1", "PCoA2")
ordination_scores_clr$Sample <- rownames(ordination_scores_clr)

# Merge with metadata
metadata <- as(sample_data(clr_physeq), "data.frame")
metadata$Sample <- rownames(metadata)
plot_data_clr <- merge(ordination_scores_clr, metadata, by = "Sample")

# Define custom aesthetics - THREE TREATMENT LINE TYPES
treatment_linetypes <- c(
  "Sterilized Standard" = "solid",
  "Sterilized Optimized" = "dashed", 
  "Washed" = "dotted"
)

clamp_shapes <- c("Clamp" = 16, "No Clamp" = 17)
plant_colors <- c("Buttercup" = "#E69F00", "Holcus" = "#56B4E9", "White clover" = "#009E73")

# Plot 1: Main CLR ordination with THREE treatment line types - SIZE 8 POINTS & THICKER LINES
p_clr_ordination <- ggplot(plot_data_clr, aes(x = PCoA1, y = PCoA2)) +
  # Points: Color by Plant, Shape by ClampStatus - SIZE 8
  geom_point(aes(color = Plant, shape = ClampStatus), size = 8, alpha = 0.8) +
  
  # Ellipses: Linetype by TreatmentSimple, Color by Plant - INCREASED LINE SIZE
  stat_ellipse(aes(linetype = TreatmentSimple, color = Plant, group = interaction(Plant, TreatmentSimple)), 
               type = "t", level = 0.7, size = 1.5, alpha = 0.7) +  # Increased from 0.8 to 1.5
  
  # Scales
  scale_color_manual(values = plant_colors, name = "Plant Species") +
  scale_shape_manual(values = clamp_shapes, name = "Clamp Status") +
  scale_linetype_manual(values = treatment_linetypes, name = "Treatment") +
  
  # Theme
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.box = "vertical"
  ) +
  labs(
    title = "Beta Diversity with CLR Transformation",
    subtitle = "PCoA - Euclidean Distance | Plants (Colors), Clamp Status (Shapes), Treatments (Line Types)",
    x = paste0("PCoA 1 (", variance_explained_clr[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained_clr[2], "%)")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 8)),  # Increased to match size 8
    shape = guide_legend(override.aes = list(size = 8)),              # Increased to match size 8
    linetype = guide_legend(override.aes = list(color = "black", size = 1.5))  # Increased from 1 to 1.5
  )

print(p_clr_ordination)

# =============================================================================
# FACETED PLOTS FOR BETTER VISUALIZATION - SIZE 8 POINTS & THICKER LINES
# =============================================================================

# Plot 2: Faceted by Treatment - SIZE 8 POINTS & THICKER LINES
p_faceted_treatment <- ggplot(plot_data_clr, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Plant, shape = ClampStatus), size = 8, alpha = 0.8) +  # SIZE 8
  stat_ellipse(aes(color = Plant, group = Plant), 
               type = "t", level = 0.7, size = 1.2, alpha = 0.7) +  # Increased from 0.6 to 1.2
  facet_wrap(~ TreatmentSimple, nrow = 1) +
  scale_color_manual(values = plant_colors, name = "Plant Species") +
  scale_shape_manual(values = clamp_shapes, name = "Clamp Status") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 13, face = "bold")
  ) +
  labs(
    title = "Beta Diversity by Treatment Type (CLR-transformed)",
    subtitle = "PCoA - Euclidean Distance | Each treatment shown separately",
    x = paste0("PCoA 1 (", variance_explained_clr[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained_clr[2], "%)")
  )

print(p_faceted_treatment)

# Plot 3: Faceted by Plant - SIZE 8 POINTS & THICKER LINES
p_faceted_plant <- ggplot(plot_data_clr, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(shape = ClampStatus, color = TreatmentSimple), size = 8, alpha = 0.8) +  # SIZE 8
  stat_ellipse(aes(linetype = TreatmentSimple, color = TreatmentSimple, group = TreatmentSimple), 
               type = "t", level = 0.7, size = 1.2, alpha = 0.7) +  # Increased from 0.6 to 1.2
  facet_wrap(~ Plant, nrow = 1) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = "Treatment") +
  scale_shape_manual(values = clamp_shapes, name = "Clamp Status") +
  scale_linetype_manual(values = treatment_linetypes, name = "Treatment") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 13, face = "bold")
  ) +
  labs(
    title = "Beta Diversity by Plant Species (CLR-transformed)",
    subtitle = "PCoA - Euclidean Distance | Each plant species shown separately",
    x = paste0("PCoA 1 (", variance_explained_clr[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained_clr[2], "%)")
  )

print(p_faceted_plant)

# =============================================================================
# COMPARISON WITH BRAY-CURTIS - SIZE 8 POINTS & THICKER LINES
# =============================================================================

print("=== COMPARISON WITH BRAY-CURTIS ===")

# Calculate Bray-Curtis on original data for comparison
bray_dist <- phyloseq::distance(clean_physeq, method = "bray")
pcoa_bray <- ordinate(clean_physeq, method = "PCoA", distance = bray_dist)
variance_explained_bray <- round(100 * pcoa_bray$values$Eigenvalues / sum(pcoa_bray$values$Eigenvalues), 2)

# Extract Bray-Curtis ordination scores
ordination_scores_bray <- as.data.frame(pcoa_bray$vectors[, 1:2])
colnames(ordination_scores_bray) <- c("PCoA1", "PCoA2")
ordination_scores_bray$Sample <- rownames(ordination_scores_bray)
plot_data_bray <- merge(ordination_scores_bray, metadata, by = "Sample")

# Plot comparison - SIZE 8 POINTS & THICKER LINES
p_bray <- ggplot(plot_data_bray, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Plant, shape = ClampStatus), size = 8, alpha = 0.8) +  # SIZE 8
  stat_ellipse(aes(linetype = TreatmentSimple, color = Plant, group = interaction(Plant, TreatmentSimple)), 
               type = "t", level = 0.7, size = 1.2, alpha = 0.7) +  # Increased from 0.6 to 1.2
  scale_color_manual(values = plant_colors) +
  scale_shape_manual(values = clamp_shapes) +
  scale_linetype_manual(values = treatment_linetypes) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "Bray-Curtis",
    x = paste0("PCoA 1 (", variance_explained_bray[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained_bray[2], "%)")
  )

p_clr_comp <- ggplot(plot_data_clr, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Plant, shape = ClampStatus), size = 8, alpha = 0.8) +  # SIZE 8
  stat_ellipse(aes(linetype = TreatmentSimple, color = Plant, group = interaction(Plant, TreatmentSimple)), 
               type = "t", level = 0.7, size = 1.2, alpha = 0.7) +  # Increased from 0.6 to 1.2
  scale_color_manual(values = plant_colors) +
  scale_shape_manual(values = clamp_shapes) +
  scale_linetype_manual(values = treatment_linetypes) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "none",
    plot.title = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "CLR + Euclidean",
    x = paste0("PCoA 1 (", variance_explained_clr[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained_clr[2], "%)")
  )

# Combine comparison plots
library(patchwork)
p_comparison <- p_bray + p_clr_comp + 
  plot_annotation(
    title = "Comparison: Bray-Curtis vs CLR + Euclidean Distance",
    subtitle = "Three Treatments: Sterilized Standard, Sterilized Optimized, Washed",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 14))
  ) +
  plot_layout(guides = "collect")

print(p_comparison)

# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

print("=== STATISTICAL ANALYSIS ===")

# PERMANOVA on CLR Euclidean distance
permanova_clr_treatment <- adonis2(clr_euclidean_dist ~ TreatmentSimple, data = metadata, permutations = 999)
permanova_clr_plant <- adonis2(clr_euclidean_dist ~ Plant, data = metadata, permutations = 999)
permanova_clr_clamp <- adonis2(clr_euclidean_dist ~ ClampStatus, data = metadata, permutations = 999)

print("PERMANOVA Results (CLR + Euclidean):")
print("Treatment:")
print(permanova_clr_treatment)
print("Plant Species:")
print(permanova_clr_plant)
print("Clamp Status:")
print(permanova_clr_clamp)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save all plots
ggsave("beta_diversity_clr_3treatments_size8_thicklines.png", p_clr_ordination, width = 12, height = 8, dpi = 300)
ggsave("beta_diversity_faceted_treatment_size8_thicklines.png", p_faceted_treatment, width = 14, height = 6, dpi = 300)
ggsave("beta_diversity_faceted_plant_size8_thicklines.png", p_faceted_plant, width = 14, height = 6, dpi = 300)
ggsave("comparison_bray_vs_clr_3treatments_size8_thicklines.png", p_comparison, width = 16, height = 8, dpi = 300)

# Save statistical results
clr_results <- list(
  PERMANOVA_Treatment = permanova_clr_treatment,
  PERMANOVA_Plant = permanova_clr_plant,
  PERMANOVA_Clamp = permanova_clr_clamp,
  Variance_Explained = variance_explained_clr
)

save(clr_results, file = "clr_analysis_results_3treatments.RData")

print("CLR transformation and analysis complete!")
print("Three treatment groups:")
print("1. Sterilized Standard (Solid lines)")
print("2. Sterilized Optimized (Dashed lines)") 
print("3. Washed (Dotted lines)")
print("All results saved to files with SIZE 8 points and THICKER LINES.")