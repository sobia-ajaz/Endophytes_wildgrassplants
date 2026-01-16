# Load required libraries
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggtext)

# Set working directory to current script location
setwd(dirname(getActiveDocumentContext()$path))

# Load QIIME2 data
ASVs <- read_qza("feature-table-16S-Plant-all.qza")
metadata <- read_q2metadata("metadata_Plant.txt")
taxonomy <- read_qza("taxonomy_SILVA-16S-Plant-all.qza")
taxonomy <- parse_taxonomy(taxonomy$data)

# Create phyloseq object
physeq <- qza_to_phyloseq(
  features = "feature-table-16S-Plant-all.qza",
  tree = "rooted-tree-16S-Plant-all.qza",
  taxonomy = "taxonomy_SILVA-16S-Plant-all.qza",
  metadata = "metadata_Plant.txt"
)

# Rarefy to even depth
set.seed(1)
S16_scaled <- rarefy_even_depth(physeq, sample.size = 40000, replace = FALSE, rngseed = 1)

# Subset to specific plants
plant_subset_physeq <- subset_samples(S16_scaled, Plant %in% c("Buttercup", "Holcus", "White clover"))

# Aggregate at Family level
family_physeq <- tax_glom(plant_subset_physeq, taxrank = "Family")

# REMOVE CHLOROPLAST AND MITOCHONDRIA
# Filter out chloroplast and mitochondria taxa
family_physeq_microbes <- subset_taxa(family_physeq, 
                                      !grepl("Chloroplast", Family, ignore.case = TRUE) & 
                                        !grepl("Mitochondria", Family, ignore.case = TRUE))

cat("After removing chloroplast and mitochondria:\n")
cat("Original number of taxa:", ntaxa(family_physeq), "\n")
cat("Microbial taxa remaining:", ntaxa(family_physeq_microbes), "\n")

# Create relative abundance version for top microbial families
family_physeq_microbes_rel <- transform_sample_counts(family_physeq_microbes, function(x) x / sum(x))
top20_microbial_families <- names(sort(taxa_sums(family_physeq_microbes_rel), decreasing = TRUE)[1:20])
family_physeq_top20_microbes <- prune_taxa(top20_microbial_families, family_physeq_microbes_rel)

# Ensure factors are properly ordered
sample_data(family_physeq_top20_microbes)$Clamp <- factor(
  sample_data(family_physeq_top20_microbes)$Clamp,
  levels = c("clamp", "no_clamp")
)

# CORRECTED: Create proper Treatment_agg variable that captures all categories
sample_data(family_physeq_top20_microbes)$Treatment_agg <- NA

for(i in 1:nsamples(family_physeq_top20_microbes)) {
  s_name <- as.character(sample_data(family_physeq_top20_microbes)$S_Name[i])
  clamp_status <- as.character(sample_data(family_physeq_top20_microbes)$Clamp[i])
  
  # Check for Clamp samples first (with optimized/standard distinction)
  if(clamp_status == "clamp") {
    if(grepl("optimized", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Clamp_Optimized_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Clamp_Optimized_Washed"
      }
    } else if(grepl("standard", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Clamp_Standard_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Clamp_Standard_Washed"
      }
    } else {
      # Generic clamp samples (just clamp + washed/sterilized)
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Clamp_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Clamp_Washed"
      }
    }
  } 
  # Check for No Clamp samples
  else if(clamp_status == "no_clamp") {
    if(grepl("sterilized", s_name, ignore.case = TRUE)) {
      sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "NoClamp_Sterilized"
    } else if(grepl("washed", s_name, ignore.case = TRUE)) {
      sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "NoClamp_Washed"
    }
  }
  
  # If still not assigned, use a default
  if(is.na(sample_data(family_physeq_top20_microbes)$Treatment_agg[i])) {
    sample_data(family_physeq_top20_microbes)$Treatment_agg[i] <- "Other"
  }
}

# CORRECTED ORDER: Sterilized first, then Washed
desired_order <- c(
  # Clamp Sterilized treatments
  "Clamp_Optimized_Sterilized",
  "Clamp_Standard_Sterilized", 
  "Clamp_Sterilized",
  # No Clamp Sterilized treatments
  "NoClamp_Sterilized",
  # Clamp Washed treatments
  "Clamp_Optimized_Washed",
  "Clamp_Standard_Washed",
  "Clamp_Washed",
  # No Clamp Washed treatments
  "NoClamp_Washed",
  "Other"
)

# Use only the categories that actually exist in our data
existing_treatments <- unique(sample_data(family_physeq_top20_microbes)$Treatment_agg)
final_order <- desired_order[desired_order %in% existing_treatments]

sample_data(family_physeq_top20_microbes)$Treatment_agg <- factor(
  sample_data(family_physeq_top20_microbes)$Treatment_agg,
  levels = final_order
)

# Melt data for plotting
df_microbes <- psmelt(family_physeq_top20_microbes)

# Colorblind-friendly palette for families
colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73",  "#66CCEE", 
  "#D55E00", "#999999", "#CC79A7", "#0072B2", 
  "#E69F99", "#F0E442", "#B15928", "#FF7F00", "#6A3D9A",
  "#FFFF99", "#1F78B4", "#33A02C", "#B2DF8A", 
  "#5D5D5D", "#882255", "#FB9A99"
)

# Create dynamic color assignment for families
num_families <- length(unique(df_microbes$Family))
if(num_families <= 20) {
  family_colors <- colorblind_palette[1:num_families]
} else {
  family_colors <- colorRampPalette(colorblind_palette)(num_families)
}

# Create function to determine treatment type for coloring
get_treatment_type <- function(treatment) {
  if(grepl("Sterilized", treatment)) {
    return("Sterilized")
  } else if(grepl("Washed", treatment)) {
    return("Washed")
  } else {
    return("Other")
  }
}

# Add treatment type to the data
df_microbes$Treatment_Type <- sapply(df_microbes$Treatment_agg, get_treatment_type)

# PLOT: Top 20 Microbial Families (without chloroplast/mitochondria)
microbial_plot <- ggplot(df_microbes, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Treatment_Type + Plant, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Treatment Groups",
    y = "Relative Abundance",
    title = "Top 20 Microbial Family Composition (Chloroplast & Mitochondria Removed)",
    fill = "Bacterial Family"
  ) +
  scale_fill_manual(values = family_colors) +
  scale_y_continuous(labels = scales::percent_format())

print(microbial_plot)

# Save microbial plot
ggsave("microbial_family_abundance_plot.png", microbial_plot, width = 16, height = 10, dpi = 700)

# Show the top 20 microbial families
cat("\nTop 20 Microbial Families (after removing chloroplast & mitochondria):\n")
top_families_summary <- df_microbes %>%
  group_by(Family) %>%
  summarise(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

print(top_families_summary)

# Create a summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total samples:", nsamples(family_physeq_top20_microbes), "\n")
cat("Total microbial taxa (after filtering):", ntaxa(family_physeq_microbes), "\n")
cat("Top 20 microbial families plotted:", nrow(top_families_summary), "\n")
cat("Treatment order:", paste(final_order, collapse = " -> "), "\n")
cat("Plants analyzed:", paste(unique(df_microbes$Plant), collapse = ", "), "\n")

# Optional: Create a version with individual samples
individual_microbial_plot <- ggplot(df_microbes, aes(x = S_Name, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Treatment_Type + Plant, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Individual Samples",
    y = "Relative Abundance",
    title = "Top 20 Microbial Family Composition by Individual Sample",
    subtitle = "Chloroplast & Mitochondria Removed",
    fill = "Bacterial Family"
  ) +
  scale_fill_manual(values = family_colors) +
  scale_y_continuous(labels = scales::percent_format())

print(individual_microbial_plot)

# Save individual sample plot
ggsave("microbial_family_individual_samples.png", individual_microbial_plot, width = 20, height = 10, dpi = 700)