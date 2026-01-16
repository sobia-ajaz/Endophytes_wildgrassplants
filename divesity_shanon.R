library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(ggpubr)

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

# Calculate Shannon diversity
shannon_data <- estimate_richness(clean_physeq, measures = "Shannon")
metadata_df <- as.data.frame(sample_data(clean_physeq))

# Combine with metadata
alpha_data <- cbind(metadata_df, shannon_data)

# Create a combined treatment variable focusing on clamp conditions
alpha_data <- alpha_data %>%
  mutate(
    Clamp_Condition = case_when(
      Clamp == "clamp" & sequencing == "Standard" ~ "Standard Clamp",
      Clamp == "clamp" & sequencing == "Optimized" ~ "Customized Clamp", 
      Clamp == "no_clamp" ~ "No Clamp",
      TRUE ~ "Other"
    )
  )

# Remove any "Other" category
alpha_data <- alpha_data %>% filter(Clamp_Condition != "Other")

# Set factor order
alpha_data$Clamp_Condition <- factor(alpha_data$Clamp_Condition,
                                     levels = c("Standard Clamp", "Customized Clamp", "No Clamp"))

# MANUAL APPROACH for summary statistics (no n() function)
plants <- unique(alpha_data$Plant)
clamp_conditions <- unique(alpha_data$Clamp_Condition)

summary_stats <- data.frame()

for (plant in plants) {
  for (condition in clamp_conditions) {
    subset_data <- alpha_data %>% 
      filter(Plant == plant, Clamp_Condition == condition)
    
    if (nrow(subset_data) > 0) {
      n_samples <- nrow(subset_data)
      mean_shannon <- round(mean(subset_data$Shannon, na.rm = TRUE), 3)
      sd_shannon <- round(sd(subset_data$Shannon, na.rm = TRUE), 3)
      se_shannon <- round(sd_shannon / sqrt(n_samples), 3)
      
      summary_stats <- rbind(summary_stats, data.frame(
        Plant = plant,
        Clamp_Condition = condition,
        n_samples = n_samples,
        mean_shannon = mean_shannon,
        sd_shannon = sd_shannon,
        se_shannon = se_shannon
      ))
    }
  }
}

print("Shannon Diversity Summary Statistics:")
print(summary_stats)

# Plot: Compare clamp conditions within each plant
p_shannon <- ggplot(alpha_data, aes(x = Clamp_Condition, y = Shannon, fill = Clamp_Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = Clamp_Condition)) +
  stat_compare_means(method = "kruskal.test", label = "p.format", 
                     label.y = max(alpha_data$Shannon) + 0.2) +
  facet_wrap(~ Plant, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Clamp Condition",
    y = "Shannon Diversity Index",
    title = "Shannon Diversity: Comparison of Clamp Conditions by Plant Species",
    subtitle = "Chloroplast and Mitochondria ASVs Removed"
  )

print(p_shannon)

# Add pairwise comparisons within each plant
p_shannon_pairwise <- ggplot(alpha_data, aes(x = Clamp_Condition, y = Shannon, fill = Clamp_Condition)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = Clamp_Condition)) +
  stat_compare_means(comparisons = list(c("Standard Clamp", "Customized Clamp"),
                                        c("Standard Clamp", "No Clamp"),
                                        c("Customized Clamp", "No Clamp")),
                     method = "wilcox.test", 
                     label = "p.signif",
                     step.increase = 0.1) +
  facet_wrap(~ Plant, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Clamp Condition",
    y = "Shannon Diversity Index",
    title = "Shannon Diversity with Pairwise Comparisons",
    subtitle = "Chloroplast and Mitochondria ASVs Removed"
  )

print(p_shannon_pairwise)

# Statistical tests for each plant
plants <- unique(alpha_data$Plant)
stat_results <- list()

for (plant in plants) {
  plant_data <- alpha_data %>% filter(Plant == plant)
  
  # Kruskal-Wallis test
  kruskal_test <- kruskal.test(Shannon ~ Clamp_Condition, data = plant_data)
  
  # Pairwise Wilcoxon tests
  pairwise_wilcox <- pairwise.wilcox.test(plant_data$Shannon, 
                                          plant_data$Clamp_Condition,
                                          p.adjust.method = "BH")
  
  stat_results[[plant]] <- list(
    kruskal_p = kruskal_test$p.value,
    pairwise = pairwise_wilcox$p.value
  )
  
  print(paste("Statistical results for", plant, ":"))
  print(paste("Kruskal-Wallis p-value:", round(kruskal_test$p.value, 4)))
  print("Pairwise Wilcoxon p-values:")
  print(pairwise_wilcox$p.value)
  print("---")
}

# Save the plot
ggsave("shannon_diversity_clamp_comparison.png", p_shannon_pairwise, 
       width = 12, height = 6, dpi = 300)

print("Analysis complete!")