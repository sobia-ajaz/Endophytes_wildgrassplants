library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(vegan)
library(rstudioapi)

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

# Create treatment groups with separate variables for plotting
sample_data(clean_physeq)$TreatmentGroup <- case_when(
  sample_data(clean_physeq)$Clamp == "clamp" & sample_data(clean_physeq)$Treatment == "Sterilized" & 
    sample_data(clean_physeq)$sequencing == "Standard" ~ "Sterilized Standard Clamp",
  sample_data(clean_physeq)$Clamp == "clamp" & sample_data(clean_physeq)$Treatment == "Sterilized" & 
    sample_data(clean_physeq)$sequencing == "Optimized" ~ "Sterilized Optimized Clamp", 
  sample_data(clean_physeq)$Clamp == "clamp" & sample_data(clean_physeq)$Treatment == "Washed" ~ "Washed Clamp",
  sample_data(clean_physeq)$Clamp == "no_clamp" & sample_data(clean_physeq)$Treatment == "Sterilized" ~ "Sterilized No Clamp",
  sample_data(clean_physeq)$Clamp == "no_clamp" & sample_data(clean_physeq)$Treatment == "Washed" ~ "Washed No Clamp",
  TRUE ~ "Other"
)

# Remove any "Other" category
clean_physeq <- subset_samples(clean_physeq, TreatmentGroup != "Other")

# Create separate variables for plotting
sample_data(clean_physeq)$ClampStatus <- ifelse(sample_data(clean_physeq)$Clamp == "clamp", "Clamp", "No Clamp")
sample_data(clean_physeq)$TreatmentType <- case_when(
  sample_data(clean_physeq)$TreatmentGroup %in% c("Sterilized Standard Clamp", "Sterilized Optimized Clamp", "Sterilized No Clamp") ~ "Sterilized",
  sample_data(clean_physeq)$TreatmentGroup %in% c("Washed Clamp", "Washed No Clamp") ~ "Washed"
)

# Set factor orders
sample_data(clean_physeq)$TreatmentGroup <- factor(sample_data(clean_physeq)$TreatmentGroup,
                                                   levels = c("Sterilized Standard Clamp", "Sterilized Optimized Clamp", 
                                                              "Sterilized No Clamp", "Washed Clamp", "Washed No Clamp"))
sample_data(clean_physeq)$Plant <- factor(sample_data(clean_physeq)$Plant,
                                          levels = c("Buttercup", "Holcus", "White clover"))
sample_data(clean_physeq)$ClampStatus <- factor(sample_data(clean_physeq)$ClampStatus,
                                                levels = c("Clamp", "No Clamp"))
sample_data(clean_physeq)$TreatmentType <- factor(sample_data(clean_physeq)$TreatmentType,
                                                  levels = c("Sterilized", "Washed"))

# =============================================================================
# BETA DIVERSITY ANALYSIS WITH CUSTOM VISUALIZATION
# =============================================================================

print("=== BETA DIVERSITY ANALYSIS WITH CUSTOM VISUALIZATION ===")

# Calculate Bray-Curtis distance
bray_dist <- phyloseq::distance(clean_physeq, method = "bray")
print("Bray-Curtis distance matrix calculated")

# PCoA ordination
pcoa <- ordinate(clean_physeq, method = "PCoA", distance = bray_dist)

# Calculate variance explained
variance_explained <- round(100 * pcoa$values$Eigenvalues / sum(pcoa$values$Eigenvalues), 2)

print(paste("Variance explained - PCoA1:", variance_explained[1], "%, PCoA2:", variance_explained[2], "%"))

# Extract ordination scores for custom plotting
ordination_scores <- as.data.frame(pcoa$vectors[, 1:2])
colnames(ordination_scores) <- c("PCoA1", "PCoA2")
ordination_scores$Sample <- rownames(ordination_scores)

# Merge with metadata
metadata <- as(sample_data(clean_physeq), "data.frame")
metadata$Sample <- rownames(metadata)
plot_data <- merge(ordination_scores, metadata, by = "Sample")

# Define custom linetypes for treatment groups
treatment_linetypes <- c(
  "Sterilized Standard Clamp" = "solid",
  "Sterilized Optimized Clamp" = "dashed", 
  "Sterilized No Clamp" = "dotted",
  "Washed Clamp" = "solid",
  "Washed No Clamp" = "dashed"
)

# Define shapes for clamp status
clamp_shapes <- c("Clamp" = 16, "No Clamp" = 17)  # 16 = circle, 17 = triangle

# Define colors for plants
plant_colors <- c("Buttercup" = "#E69F00", "Holcus" = "#56B4E9", "White clover" = "#009E73")

# Plot 1: Main ordination with custom aesthetics
p_custom_ordination <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  # Points: Color by Plant, Shape by ClampStatus
  geom_point(aes(color = Plant, shape = ClampStatus), size = 4, alpha = 0.8) +
  
  # Ellipses: Linetype by TreatmentGroup, Color by Plant
  stat_ellipse(aes(linetype = TreatmentGroup, color = Plant, group = interaction(Plant, TreatmentGroup)), 
               type = "t", level = 0.7, size = 0.8, alpha = 0.7) +
  
  # Scales
  scale_color_manual(values = plant_colors, name = "Plant Species") +
  scale_shape_manual(values = clamp_shapes, name = "Clamp Status") +
  scale_linetype_manual(values = treatment_linetypes, name = "Treatment Group") +
  
  # Theme
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    legend.box = "vertical",
    legend.spacing.y = unit(0.5, "cm")
  ) +
  labs(
    title = "Beta Diversity: Plants (Colors), Clamp Status (Shapes), Treatments (Line Types)",
    subtitle = "PCoA - Bray-Curtis | Chloroplast and Mitochondria ASVs Removed",
    x = paste0("PCoA 1 (", variance_explained[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained[2], "%)")
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 5)),  # Squares for color legend
    shape = guide_legend(override.aes = list(size = 4)),
    linetype = guide_legend(override.aes = list(color = "black", size = 1))
  )

print(p_custom_ordination)

# Plot 2: Faceted by Treatment Type (Sterilized vs Washed)
p_faceted_treatment <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Plant, shape = ClampStatus), size = 3, alpha = 0.8) +
  stat_ellipse(aes(linetype = TreatmentGroup, color = Plant, group = interaction(Plant, TreatmentGroup)), 
               type = "t", level = 0.7, size = 0.6, alpha = 0.7) +
  facet_wrap(~ TreatmentType, nrow = 1) +
  scale_color_manual(values = plant_colors, name = "Plant Species") +
  scale_shape_manual(values = clamp_shapes, name = "Clamp Status") +
  scale_linetype_manual(values = treatment_linetypes, name = "Treatment Group") +
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
    title = "Beta Diversity by Treatment Type",
    subtitle = "PCoA - Bray-Curtis | Chloroplast and Mitochondria ASVs Removed",
    x = paste0("PCoA 1 (", variance_explained[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained[2], "%)")
  )

print(p_faceted_treatment)

# Plot 3: Simplified version focusing on main effects
p_simplified <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Plant, shape = ClampStatus), size = 4, alpha = 0.8) +
  stat_ellipse(aes(linetype = TreatmentType, color = Plant, group = interaction(Plant, TreatmentType)), 
               type = "t", level = 0.7, size = 0.8, alpha = 0.7) +
  scale_color_manual(values = plant_colors, name = "Plant Species") +
  scale_shape_manual(values = clamp_shapes, name = "Clamp Status") +
  scale_linetype_manual(values = c("Sterilized" = "solid", "Washed" = "dashed"), name = "Treatment Type") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "Beta Diversity: Simplified View",
    subtitle = "Plants (Colors), Clamp Status (Shapes), Treatment Type (Line Types)",
    x = paste0("PCoA 1 (", variance_explained[1], "%)"),
    y = paste0("PCoA 2 (", variance_explained[2], "%)")
  )

print(p_simplified)

# =============================================================================
# PERMANOVA STATISTICAL TESTS
# =============================================================================

print("=== PERMANOVA STATISTICAL ANALYSIS ===")

# PERMANOVA for main effects
permanova_treatment <- adonis2(bray_dist ~ TreatmentGroup, data = metadata, permutations = 999)
permanova_plant <- adonis2(bray_dist ~ Plant, data = metadata, permutations = 999)
permanova_clamp <- adonis2(bray_dist ~ ClampStatus, data = metadata, permutations = 999)
permanova_interaction <- adonis2(bray_dist ~ TreatmentGroup * Plant, data = metadata, permutations = 999)

print("PERMANOVA Results:")
print("Treatment Group:")
print(permanova_treatment)
print("Plant Species:")
print(permanova_plant)
print("Clamp Status:")
print(permanova_clamp)
print("Interaction (Treatment Group * Plant):")
print(permanova_interaction)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save all plots
ggsave("beta_diversity_custom.png", p_custom_ordination, width = 12, height = 8, dpi = 300)
ggsave("beta_diversity_faceted_treatment.png", p_faceted_treatment, width = 14, height = 6, dpi = 300)
ggsave("beta_diversity_simplified.png", p_simplified, width = 10, height = 8, dpi = 300)

# Save PERMANOVA results
permanova_results <- list(
  TreatmentGroup = permanova_treatment,
  Plant = permanova_plant,
  ClampStatus = permanova_clamp,
  Interaction = permanova_interaction
)

save(permanova_results, file = "permanova_results_custom.RData")

print("Custom beta diversity analysis complete!")
print("Visualization features:")
print("- Plants: Different colors (Buttercup, Holcus, White clover)")
print("- Clamp Status: Different shapes (Clamp = circle, No Clamp = triangle)")
print("- Treatment Groups: Different line types")
print("All plots and results saved to files.")