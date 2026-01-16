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

# Add treatment type to the data
df_all$Treatment_Type <- sapply(df_all$Treatment_agg, get_treatment_type)

# PLOT 1: Family-level composition
family_plot <- ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
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
    title = "Top 20 Family-level Microbial Composition",
    fill = "Bacterial Family"
  ) +
  scale_fill_manual(values = family_colors) +
  scale_y_continuous(labels = scales::percent_format())

print(family_plot)

# Save family plot
ggsave("family_abundance_plot.png", family_plot, width = 16, height = 10, dpi = 700)

##########################################################
# Categorize families into Chloroplast, Mitochondria, and Bacteria
########################################################

# SIMPLER APPROACH: Use psmelt which already includes sample data
df_all_categories <- psmelt(family_physeq)

# Check what columns are available
cat("Columns in df_all_categories:\n")
print(colnames(df_all_categories))

# Categorize families into Chloroplast, Mitochondria, and Bacteria
df_all_categories$Category <- ifelse(grepl("Chloroplast", df_all_categories$Family, ignore.case = TRUE), "Chloroplast",
                                     ifelse(grepl("Mitochondria", df_all_categories$Family, ignore.case = TRUE), "Mitochondria", "Bacteria"))

# Check if Treatment_agg exists, if not create it from S_Name
if(!"Treatment_agg" %in% colnames(df_all_categories)) {
  cat("Treatment_agg not found, creating from S_Name...\n")
  df_all_categories$Treatment_agg <- NA
  
  for(i in 1:nrow(df_all_categories)) {
    s_name <- as.character(df_all_categories$S_Name[i])
    
    if(grepl("optimized", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "Clamp_Optimized_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "Clamp_Optimized_Washed"
      }
    } else if(grepl("standard", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "Clamp_Standard_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "Clamp_Standard_Washed"
      }
    } else if(grepl("clamp", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "Clamp_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "Clamp_Washed"
      }
    } else if(grepl("no.?clamp", s_name, ignore.case = TRUE)) {
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "NoClamp_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        df_all_categories$Treatment_agg[i] <- "NoClamp_Washed"
      }
    } else {
      df_all_categories$Treatment_agg[i] <- "Other"
    }
  }
  
  # Set factor levels
  df_all_categories$Treatment_agg <- factor(df_all_categories$Treatment_agg, levels = final_order)
}

# Aggregate abundance values for categories
df_agg <- df_all_categories %>%
  group_by(Sample, S_Name, Treatment_agg, Plant, Category) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Calculate relative abundance
df_agg <- df_agg %>%
  group_by(Sample) %>%
  mutate(Abundance_relative = Abundance / sum(Abundance) * 100)

# Define custom colors for categories
category_palette <- c("Chloroplast" = '#009E73', "Mitochondria" = '#1F78B4', "Bacteria" = '#FF7F00')

# Ensure proper factor order
df_agg$Category <- factor(df_agg$Category, levels = c("Chloroplast", "Mitochondria", "Bacteria"))

# Add treatment type for faceting
df_agg$Treatment_Type <- sapply(df_agg$Treatment_agg, get_treatment_type)

# PLOT 2: Three-category composition (same style as family plot)
category_plot <- ggplot(df_agg, aes(x = Treatment_agg, y = Abundance_relative/100, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~ Treatment_Type + Plant, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Treatment Groups",
    y = "Relative Abundance",
    title = "Microbial Composition by Category",
    fill = "Category"
  ) +
  scale_fill_manual(values = category_palette) +
  scale_y_continuous(labels = scales::percent_format())

print(category_plot)

# Save category plot
ggsave("category_abundance_plot.png", category_plot, width = 16, height = 10, dpi = 700)

# PLOT 3: Three-category composition with percentage labels (alternative version)
# First, calculate the percentage for each category within each treatment
df_agg_percent <- df_agg %>%
  group_by(Treatment_agg, Plant, Category) %>%
  summarise(Mean_Abundance = mean(Abundance_relative), .groups = "drop")

# Add treatment type for faceting
df_agg_percent$Treatment_Type <- sapply(df_agg_percent$Treatment_agg, get_treatment_type)

category_plot_percent <- ggplot(df_agg_percent, aes(x = Treatment_agg, y = Mean_Abundance/100, fill = Category)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(
    aes(label = ifelse(Mean_Abundance > 5, sprintf("%.0f%%", Mean_Abundance), "")), 
    position = position_fill(vjust = 0.5), 
    size = 3, color = "black", fontface = "bold"
  ) +
  facet_grid(~ Treatment_Type + Plant, scales = "free_x", space = "free") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10,face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Treatment",
    y = "Relative Abundance",
    title = "",
    subtitle = "",
    fill = "Category"
  ) +
  scale_fill_manual(values = category_palette) +
  scale_y_continuous(labels = scales::percent_format())

print(category_plot_percent)

# Save percentage label version
ggsave("category_abundance_plot_with_labels.png", category_plot_percent, width = 16, height = 10, dpi = 700)

# Create a summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total samples:", nsamples(family_physeq_top20), "\n")
cat("Treatment order:", paste(final_order, collapse = " -> "), "\n")
cat("Plants analyzed:", paste(unique(df_all$Plant), collapse = ", "), "\n")
cat("Plots created:\n")
cat("- family_abundance_plot.png: Top 20 family-level composition\n")
cat("- category_abundance_plot.png: Three-category composition\n")
cat("- category_abundance_plot_with_labels.png: Three-category with percentage labels\n")