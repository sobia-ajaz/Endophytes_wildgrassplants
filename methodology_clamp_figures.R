# Load required libraries
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(vegan)
library(gtable)
library(grid)
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

# Create relative abundance version for top families
family_physeq_rel <- transform_sample_counts(family_physeq, function(x) x / sum(x))
top20_families <- names(sort(taxa_sums(family_physeq_rel), decreasing = TRUE)[1:20])
family_physeq_top20 <- prune_taxa(top20_families, family_physeq_rel)

# Ensure factors are properly ordered
sample_data(family_physeq_top20)$Clamp <- factor(
  sample_data(family_physeq_top20)$Clamp,
  levels = c("clamp", "no_clamp")
)

# CORRECTED: Create proper Treatment_agg variable that captures all categories
sample_data(family_physeq_top20)$Treatment_agg <- NA

for(i in 1:nsamples(family_physeq_top20)) {
  s_name <- as.character(sample_data(family_physeq_top20)$S_Name[i])
  clamp_status <- as.character(sample_data(family_physeq_top20)$Clamp[i])
  
  # Check for Clamp samples first (with optimized/standard distinction)
  if(clamp_status == "clamp") {
    if(grepl("optimized", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Optimized_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Optimized_Washed"
      }
    } else if(grepl("standard", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Standard_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Standard_Washed"
      }
    } else {
      # Generic clamp samples (just clamp + washed/sterilized)
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Washed"
      }
    }
  } 
  # Check for No Clamp samples
  else if(clamp_status == "no_clamp") {
    if(grepl("sterilized", s_name, ignore.case = TRUE)) {
      sample_data(family_physeq_top20)$Treatment_agg[i] <- "NoClamp_Sterilized"
    } else if(grepl("washed", s_name, ignore.case = TRUE)) {
      sample_data(family_physeq_top20)$Treatment_agg[i] <- "NoClamp_Washed"
    }
  }
  
  # If still not assigned, use a default
  if(is.na(sample_data(family_physeq_top20)$Treatment_agg[i])) {
    sample_data(family_physeq_top20)$Treatment_agg[i] <- "Other"
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
existing_treatments <- unique(sample_data(family_physeq_top20)$Treatment_agg)
final_order <- desired_order[desired_order %in% existing_treatments]

sample_data(family_physeq_top20)$Treatment_agg <- factor(
  sample_data(family_physeq_top20)$Treatment_agg,
  levels = final_order
)

# Melt data for plotting
df_all <- psmelt(family_physeq_top20)

# Colorblind-friendly palette for families
colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#5D5D5D", "#66CCEE", 
  "#D55E00", "#999999", "#CC79A7", "#0072B2", "#882255",
  "#E69F99", "#F0E442", "#B15928", "#FF7F00", "#6A3D9A",
  "#FFFF99", "#1F78B4", "#33A02C", "#B2DF8A", "#FB9A99"
)

# Colors for treatment type boxes
sterilized_color <- "#E6F3FF"  # Light blue for sterilized
washed_color <- "#FFF0E6"      # Light orange for washed

# Create dynamic color assignment for families
num_families <- length(unique(df_all$Family))
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

# Create custom x-axis labels with colored backgrounds
create_treatment_labels <- function(treatment_levels) {
  labels <- c()
  for(treatment in treatment_levels) {
    treatment_type <- get_treatment_type(treatment)
    
    if(treatment == "Clamp_Optimized_Sterilized") {
      label_text <- "Clamp\nOptimised\n(Sterilised)"
    } else if(treatment == "Clamp_Standard_Sterilized") {
      label_text <- "Clamp\nStandard\n(Sterilised)"
    } else if(treatment == "Clamp_Sterilized") {
      label_text <- "Clamp\n(Sterilised)"
    } else if(treatment == "Clamp_Optimized_Washed") {
      label_text <- "Clamp\nOptimized\n(Washed)"
    } else if(treatment == "Clamp_Standard_Washed") {
      label_text <- "Clamp\nStandard\n(Washed)"
    } else if(treatment == "Clamp_Washed") {
      label_text <- "Clamp\n(Washed)"
    } else if(treatment == "NoClamp_Sterilized") {
      label_text <- "No Clamp\n(Sterilised)"
    } else if(treatment == "NoClamp_Washed") {
      label_text <- "No Clamp\n(Washed)"
    } else {
      label_text <- treatment
    }
    
    # Apply color based on treatment type
    if(treatment_type == "Sterilized") {
      labels <- c(labels, paste0("<span style='background-color:", sterilized_color, "; padding: 5px; border-radius: 3px;'>", label_text, "</span>"))
    } else if(treatment_type == "Washed") {
      labels <- c(labels, paste0("<span style='background-color:", washed_color, "; padding: 5px; border-radius: 3px;'>", label_text, "</span>"))
    } else {
      labels <- c(labels, label_text)
    }
  }
  return(labels)
}

# Create the plot with colored treatment labels
main_plot <- ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme_minimal() +
  theme(
    # Increased x-axis label size with HTML parsing
    axis.text.x = element_markdown(angle = 90, hjust = 1, size = 12),
    # Increased y-axis label size
    axis.text.y = element_text(size = 12),
    # Increased axis title sizes
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    # Increased legend text size
    legend.text = element_text(size = 11),
    # Increased legend title size
    legend.title = element_text(size = 13, face = "bold"),
    # Increased plot title size
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14),
    # Increased facet label size
    strip.text = element_text(face = "bold", size = 13),
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(
    x = "Treatment Groups",
    y = "Relative Abundance",
    title = "",
    subtitle = "",
    fill = "Bacterial Family"
  ) +
  scale_fill_manual(values = family_colors) +
  scale_x_discrete(labels = create_treatment_labels(levels(df_all$Treatment_agg)))

# Print the plot
print(main_plot)

# Alternative method using annotation if HTML doesn't work
# Create a version with manual color coding in the legend instead
cat("\nIf the colored boxes don't appear, creating alternative version...\n")

# Alternative plot with treatment type as fill and faceting
df_all$Treatment_Type <- sapply(df_all$Treatment_agg, get_treatment_type)

alternative_plot <- ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Treatment_Type + Plant, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Treatment Groups",
    y = "Relative Abundance",
    title = "",
    subtitle = "",
    fill = "Bacterial Family"
  ) +
  scale_fill_manual(values = family_colors)

print(alternative_plot)
# Save both plots
ggsave("family_abundance_plot_colored_labels.png", main_plot, width = 18, height = 10, dpi = 700)
ggsave("family_abundance_plot_grouped_treatments.png", alternative_plot, width = 18, height = 10, dpi = 700)

# Create a summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total samples:", nsamples(family_physeq_top20), "\n")
cat("Treatment order:", paste(final_order, collapse = " -> "), "\n")
cat("Plants analyzed:", paste(unique(df_all$Plant), collapse = ", "), "\n")
cat("Sterilized treatments colored with:", sterilized_color, "\n")
cat("Washed treatments colored with:", washed_color, "\n")

####################################################
##################################################
#Individual samples

alternative_plot_individual <- ggplot(df_all, aes(x = S_Name, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Treatment_Type + Plant, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Treatment Groups",
    y = "Relative Abundance",
    title = "",
    subtitle = "",
    fill = "Bacterial Family"
  ) +
  scale_fill_manual(values = family_colors)

print(alternative_plot_individual)
# Save both plots

ggsave("family_abundance_plot_grouped_treatments_individual_samples.png", alternative_plot, width = 18, height = 10, dpi = 700)

##########################################################
#Categorize families into Chloroplast, Mitochondria, and Microbes
########################################################

# Convert phyloseq object to data frame
df_all <- psmelt(family_physeq)

# Categorize families into Chloroplast, Mitochondria, and Microbes
df_all$Category <- ifelse(df_all$Family == "Chloroplast", "Chloroplast",
                          ifelse(df_all$Family == "Mitochondria", "Mitochondria", "Bacteria"))

# Aggregate abundance values for microbes
df_agg <- df_all %>%
  group_by(Sample,S_Name, Treatment_agg, Plant, Category) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")


# Define custom colors for categories
color_palette <- c("Chloroplast" = '#009E73', "Mitochondria" = '#1F78B4', "Bacteria" = '#FF7F00')

# Ensure "Microbes" is always at the bottom of the stacked bars
df_agg$Category <- factor(df_agg$Category, levels = c("Chloroplast", "Mitochondria", "Bacteria"))
df_agg <- df_agg %>%
  filter(Category %in% c("Chloroplast", "Mitochondria", "Bacteria")) %>%
  group_by(Treatment_agg, Plant) %>%
  mutate(Abundance_relative = Abundance / sum(Abundance) * 100)
df_agg_summarized <- df_agg %>%
  filter(Category %in% c("Bacteria", "Chloroplast", "Mitochondria")) %>%
  group_by(Treatment_agg, Plant, Category) %>%
  summarise(Abundance_relative = sum(Abundance_relative), .groups = "drop")  # Sum up within groups

###### plot Final aggregated 

ggplot(df_agg_summarized, aes(x = Treatment_agg, y = Abundance_relative, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    aes(label = sprintf("%.0f%%", Abundance_relative )), 
    position = position_fill(vjust = 0.5), 
    size = 3, color = "black"
  ) +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "Treatments and Clamping ",
    y = "Relative Abundance (%)\n16S rRNA Gene Sequences",
    fill = "Family"
  ) +
  scale_fill_manual(values = color_palette) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))
##




