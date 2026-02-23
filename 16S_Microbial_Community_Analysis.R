##############################################################
# 16S rRNA Microbial Community Analysis – Visualization
# Author: Dr. Sobia Ajaz
# Description: Loads rarefied & subsetted phyloseq object from RDS,
#              then generates family-level and category-level plots.
##############################################################

# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(ggtext)

# Set working directory to current script location (RStudio only)
setwd(dirname(getActiveDocumentContext()$path))

# -------------------------------------------------------------
# Load the previously saved phyloseq object (rarefied and subsetted)
# -------------------------------------------------------------
plant_subset_physeq <- readRDS("plant_subset_16S.rds")

cat("Loaded phyloseq object:\n")
print(plant_subset_physeq)

# Verify that required sample variables exist
required_vars <- c("Plant", "Clamp", "S_Name")
missing_vars <- required_vars[!required_vars %in% sample_variables(plant_subset_physeq)]
if(length(missing_vars) > 0) {
  stop("Missing sample variables: ", paste(missing_vars, collapse = ", "))
} else {
  cat("All required sample variables found.\n")
}

# -------------------------------------------------------------
# Aggregate at Family level
# -------------------------------------------------------------
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

# -------------------------------------------------------------
# Create Treatment_agg variable based on S_Name and Clamp
# -------------------------------------------------------------
sample_data(family_physeq_top20)$Treatment_agg <- NA

for(i in 1:nsamples(family_physeq_top20)) {
  s_name <- as.character(sample_data(family_physeq_top20)$S_Name[i])
  clamp_status <- as.character(sample_data(family_physeq_top20)$Clamp[i])
  
  # Clamp samples
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
      # Generic clamp samples
      if(grepl("sterilized", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Sterilized"
      } else if(grepl("washed", s_name, ignore.case = TRUE)) {
        sample_data(family_physeq_top20)$Treatment_agg[i] <- "Clamp_Washed"
      }
    }
  } 
  # No Clamp samples
  else if(clamp_status == "no_clamp") {
    if(grepl("sterilized", s_name, ignore.case = TRUE)) {
      sample_data(family_physeq_top20)$Treatment_agg[i] <- "NoClamp_Sterilized"
    } else if(grepl("washed", s_name, ignore.case = TRUE)) {
      sample_data(family_physeq_top20)$Treatment_agg[i] <- "NoClamp_Washed"
    }
  }
  
  # Fallback
  if(is.na(sample_data(family_physeq_top20)$Treatment_agg[i])) {
    sample_data(family_physeq_top20)$Treatment_agg[i] <- "Other"
  }
}

# Define desired order: Sterilized first, then Washed
desired_order <- c(
  "Clamp_Optimized_Sterilized",
  "Clamp_Standard_Sterilized", 
  "Clamp_Sterilized",
  "NoClamp_Sterilized",
  "Clamp_Optimized_Washed",
  "Clamp_Standard_Washed",
  "Clamp_Washed",
  "NoClamp_Washed",
  "Other"
)

existing_treatments <- unique(sample_data(family_physeq_top20)$Treatment_agg)
final_order <- desired_order[desired_order %in% existing_treatments]

sample_data(family_physeq_top20)$Treatment_agg <- factor(
  sample_data(family_physeq_top20)$Treatment_agg,
  levels = final_order
)

# -------------------------------------------------------------
# Melt data for plotting (family-level)
# -------------------------------------------------------------
df_all <- psmelt(family_physeq_top20)

# Colorblind-friendly palette for families
colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#5D5D5D", "#66CCEE", 
  "#D55E00", "#999999", "#CC79A7", "#0072B2", "#882255",
  "#E69F99", "#F0E442", "#B15928", "#FF7F00", "#6A3D9A",
  "#FFFF99", "#1F78B4", "#33A02C", "#B2DF8A", "#FB9A99"
)

num_families <- length(unique(df_all$Family))
if(num_families <= 20) {
  family_colors <- colorblind_palette[1:num_families]
} else {
  family_colors <- colorRampPalette(colorblind_palette)(num_families)
}

# Helper to get treatment type (Sterilized / Washed)
get_treatment_type <- function(treatment) {
  if(grepl("Sterilized", treatment)) {
    return("Sterilized")
  } else if(grepl("Washed", treatment)) {
    return("Washed")
  } else {
    return("Other")
  }
}

df_all$Treatment_Type <- sapply(df_all$Treatment_agg, get_treatment_type)

# -------------------------------------------------------------
# PLOT 1: Family-level composition
# -------------------------------------------------------------
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
ggsave("family_abundance_plot.png", family_plot, width = 16, height = 10, dpi = 700)

# -------------------------------------------------------------
# Categorize families into Chloroplast, Mitochondria, Bacteria
# -------------------------------------------------------------
# Use full family object to overview
df_all_categories <- psmelt(family_physeq)

# --- ROBUST SAMPLE COLUMN HANDLING ---
if(!"Sample" %in% colnames(df_all_categories)) {
  possible_ids <- c("Sample", "sample", "SampleID", "sampleid", "id", "ID", "rowname")
  found_id <- possible_ids[possible_ids %in% colnames(df_all_categories)]
  if(length(found_id) == 0) {
    stop("No sample identifier column found. Columns: ", paste(colnames(df_all_categories), collapse=", "))
  }
  df_all_categories <- df_all_categories %>% rename(Sample = !!found_id[1])
  message("Renamed '", found_id[1], "' to 'Sample'.")
} else {
  message("Column 'Sample' already exists.")
}
# -------------------------------------

# Assign categories
df_all_categories$Category <- ifelse(
  grepl("Chloroplast", df_all_categories$Family, ignore.case = TRUE), "Chloroplast",
  ifelse(grepl("Mitochondria", df_all_categories$Family, ignore.case = TRUE), "Mitochondria", "Bacteria")
)

# Merge Treatment_agg and other metadata from family_physeq_top20
sample_df <- data.frame(sample_data(family_physeq_top20)) %>%
  select(S_Name, Treatment_agg, Plant, Clamp)

# Join by matching sample identifier (assuming Sample column contains S_Name values)
if(!"Treatment_agg" %in% colnames(df_all_categories)) {
  df_all_categories <- df_all_categories %>%
    left_join(sample_df, by = c("Sample" = "S_Name"))
}

# Ensure Plant exists
if(!"Plant" %in% colnames(df_all_categories) && "Plant" %in% colnames(sample_df)) {
  df_all_categories <- df_all_categories %>%
    left_join(sample_df %>% select(S_Name, Plant), by = c("Sample" = "S_Name"))
}

# Also ensure S_Name exists
if(!"S_Name" %in% colnames(df_all_categories)) {
  if(all(df_all_categories$Sample %in% sample_df$S_Name)) {
    df_all_categories$S_Name <- df_all_categories$Sample
    message("Created S_Name from Sample.")
  } else {
    stop("S_Name column not found and cannot be inferred.")
  }
}

# --- DIAGNOSTIC: Check columns before aggregation ---
cat("\nColumns in df_all_categories before aggregation:\n")
print(colnames(df_all_categories))
if(!"Sample" %in% colnames(df_all_categories)) {
  stop("ERROR: 'Sample' column is missing. Check output above.")
}
# ----------------------------------------------------

# Aggregate by sample and category using base R 
df_agg <- aggregate(Abundance ~ Sample + S_Name + Treatment_agg + Plant + Category,
                    data = df_all_categories, FUN = sum)

# Now compute relative abundance per sample using dplyr
df_agg <- df_agg %>%
  group_by(Sample) %>%
  mutate(Abundance_relative = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# Define category colors
category_palette <- c("Chloroplast" = '#009E73', "Mitochondria" = '#1F78B4', "Bacteria" = '#FF7F00')
df_agg$Category <- factor(df_agg$Category, levels = c("Chloroplast", "Mitochondria", "Bacteria"))
df_agg$Treatment_Type <- sapply(df_agg$Treatment_agg, get_treatment_type)

# -------------------------------------------------------------
# PLOT 2: Three-category composition
# -------------------------------------------------------------
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
ggsave("category_abundance_plot.png", category_plot, width = 16, height = 10, dpi = 700)

# -------------------------------------------------------------
# PLOT 3: Category plot with percentage labels 
# -------------------------------------------------------------
# Use aggregate to compute mean abundance per group (retains all grouping columns)
df_agg_percent <- aggregate(Abundance_relative ~ Treatment_agg + Treatment_Type + Plant + Category,
                            data = df_agg, FUN = mean)

# Rename the mean column for clarity
names(df_agg_percent)[names(df_agg_percent) == "Abundance_relative"] <- "Mean_Abundance"

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
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    panel.spacing = unit(0.5, "lines")
  ) +
  labs(
    x = "Treatment",
    y = "Relative Abundance",
    title = "",
    fill = "Category"
  ) +
  scale_fill_manual(values = category_palette) +
  scale_y_continuous(labels = scales::percent_format())

print(category_plot_percent)
ggsave("category_abundance_plot_with_labels.png", category_plot_percent, width = 16, height = 10, dpi = 700)

# -------------------------------------------------------------
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Total samples in top20 dataset:", nsamples(family_physeq_top20), "\n")
cat("Treatment order:", paste(final_order, collapse = " -> "), "\n")
cat("Plants analyzed:", paste(unique(df_all$Plant), collapse = ", "), "\n")
cat("Plots saved:\n")
cat("- family_abundance_plot.png\n")
cat("- category_abundance_plot.png\n")
cat("- category_abundance_plot_with_labels.png\n")
