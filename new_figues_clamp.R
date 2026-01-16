# Load required libraries
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(vegan)

# Try to load ggtext for advanced text formatting, if not available use alternative method
if(!require(ggtext, quietly = TRUE)) {
  cat("ggtext package not available. Using alternative labeling method.\n")
  use_ggtext <- FALSE
} else {
  use_ggtext <- TRUE
  library(ggtext)
}

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

# METHOD 1: Using ggtext for colored labels (if available)
if(use_ggtext) {
  cat("Using ggtext for colored treatment labels...\n")
  
  create_treatment_labels <- function(treatment_levels) {
    labels <- c()
    for(treatment in treatment_levels) {
      treatment_type <- get_treatment_type(treatment)
      
      if(treatment == "Clamp_Optimized_Sterilized") {
        label_text <- "Clamp\nOptimized\n(Sterilized)"
      } else if(treatment == "Clamp_Standard_Sterilized") {
        label_text <- "Clamp\nStandard\n(Sterilized)"
      } else if(treatment == "Clamp_Sterilized") {
        label_text <- "Clamp\n(Sterilized)"
      } else if(treatment == "Clamp_Optimized_Washed") {
        label_text <- "Clamp\nOptimized\n(Washed)"
      } else if(treatment == "Clamp_Standard_Washed") {
        label_text <- "Clamp\nStandard\n(Washed)"
      } else if(treatment == "Clamp_Washed") {
        label_text <- "Clamp\n(Washed)"
      } else if(treatment == "NoClamp_Sterilized") {
        label_text <- "No Clamp\n(Sterilized)"
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
  
  main_plot <- ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~ Plant, scales = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      strip.text = element_text(face = "bold", size = 13),
      panel.spacing = unit(1.5, "lines")
    ) +
    labs(
      x = "Treatment Groups",
      y = "Relative Abundance",
      title = "Top 20 Family-level Microbial Composition by Plant Species",
      subtitle = "Sterilized treatments (light blue) vs Washed treatments (light orange)",
      fill = "Bacterial Family"
    ) +
    scale_fill_manual(values = family_colors) +
    scale_x_discrete(labels = create_treatment_labels(levels(df_all$Treatment_agg)))
  
} else {
  # METHOD 2: Alternative approach without ggtext - using faceting
  cat("Using alternative method with faceting...\n")
  
  # Add treatment type to the data
  df_all$Treatment_Type <- sapply(df_all$Treatment_agg, get_treatment_type)
  df_all$Treatment_Type <- factor(df_all$Treatment_Type, levels = c("Sterilized", "Washed", "Other"))
  
  # Create simple labels without colors
  create_simple_labels <- function(treatment_levels) {
    labels <- c()
    for(treatment in treatment_levels) {
      if(treatment == "Clamp_Optimized_Sterilized") {
        labels <- c(labels, "Clamp\nOptimized\n(Sterilized)")
      } else if(treatment == "Clamp_Standard_Sterilized") {
        labels <- c(labels, "Clamp\nStandard\n(Sterilized)")
      } else if(treatment == "Clamp_Sterilized") {
        labels <- c(labels, "Clamp\n(Sterilized)")
      } else if(treatment == "Clamp_Optimized_Washed") {
        labels <- c(labels, "Clamp\nOptimized\n(Washed)")
      } else if(treatment == "Clamp_Standard_Washed") {
        labels <- c(labels, "Clamp\nStandard\n(Washed)")
      } else if(treatment == "Clamp_Washed") {
        labels <- c(labels, "Clamp\n(Washed)")
      } else if(treatment == "NoClamp_Sterilized") {
        labels <- c(labels, "No Clamp\n(Sterilized)")
      } else if(treatment == "NoClamp_Washed") {
        labels <- c(labels, "No Clamp\n(Washed)")
      } else {
        labels <- c(labels, treatment)
      }
    }
    return(labels)
  }
  
  # Create plot with faceted treatment types
  main_plot <- ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_grid(Plant ~ Treatment_Type, scales = "free_x", space = "free") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 14),
      strip.text = element_text(face = "bold", size = 11),
      panel.spacing = unit(0.5, "lines"),
      strip.background = element_rect(fill = "lightgray")
    ) +
    labs(
      x = "Treatment Groups",
      y = "Relative Abundance",
      title = "Top 20 Family-level Microbial Composition by Plant Species",
      subtitle = "Grouped by Treatment Type (Sterilized vs Washed)",
      fill = "Bacterial Family"
    ) +
    scale_fill_manual(values = family_colors) +
    scale_x_discrete(labels = create_simple_labels(levels(df_all$Treatment_agg)))
}

# Print the plot
print(main_plot)

# Save the plot
if(use_ggtext) {
  ggsave("family_abundance_plot_colored_labels.png", main_plot, width = 18, height = 10, dpi = 300)
  cat("Plot saved as 'family_abundance_plot_colored_labels.png'\n")
} else {
  ggsave("family_abundance_plot_grouped_treatments.png", main_plot, width = 18, height = 10, dpi = 300)
  cat("Plot saved as 'family_abundance_plot_grouped_treatments.png'\n")
}

# Create a summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total samples:", nsamples(family_physeq_top20), "\n")
cat("Treatment order:", paste(final_order, collapse = " -> "), "\n")
cat("Plants analyzed:", paste(unique(df_all$Plant), collapse = ", "), "\n")
cat("Method used:", ifelse(use_ggtext, "ggtext with colored labels", "Faceting by treatment type"), "\n")