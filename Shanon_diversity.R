library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(ggpubr)

# Set working directory
setwd(dirname(getActiveDocumentContext()$path))

# Load and prepare data
physeq <- readRDS("plant_subset_16S.rds")
physeq

# Remove chloroplast and mitochondria
clean_physeq <- subset_taxa(physeq, 
                            !grepl("chloroplast|Chloroplast", Order, ignore.case = TRUE) &
                              !grepl("chloroplast|Chloroplast", Family, ignore.case = TRUE) &
                              !grepl("mitochondria|Mitochondria", Family, ignore.case = TRUE))

print(paste("Final clean data:", ntaxa(clean_physeq), "ASVs after removing chloroplast/mitochondria"))

# Calculate Shannon diversity
shannon_data <- estimate_richness(clean_physeq, measures = "Shannon")
metadata_df <- as.data.frame(sample_data(clean_physeq))

# Combine with metadata
alpha_data <- cbind(metadata_df, shannon_data)

# Create treatment groups including no clamp - SEPARATE NO CLAMP BY TREATMENT
alpha_data <- alpha_data %>%
  mutate(
    TreatmentGroup = case_when(
      Clamp == "clamp" & Treatment == "Sterilized" & sequencing == "Standard" ~ "Sterilized Standard Clamp",
      Clamp == "clamp" & Treatment == "Sterilized" & sequencing == "Optimized" ~ "Sterilized Optimized Clamp", 
      Clamp == "clamp" & Treatment == "Washed" ~ "Washed Clamp",
      Clamp == "no_clamp" & Treatment == "Sterilized" ~ "Sterilized No Clamp",
      Clamp == "no_clamp" & Treatment == "Washed" ~ "Washed No Clamp",
      TRUE ~ "Other"
    )
  )

# Remove any "Other" category
alpha_data <- alpha_data %>% filter(TreatmentGroup != "Other")

# Set factor orders - STERILIZED FIRST, THEN WASHED
alpha_data$TreatmentGroup <- factor(alpha_data$TreatmentGroup,
                                    levels = c("Sterilized Standard Clamp", "Sterilized Optimized Clamp", 
                                               "Sterilized No Clamp", "Washed Clamp", "Washed No Clamp"))
alpha_data$Plant <- factor(alpha_data$Plant,
                           levels = c("Buttercup", "Holcus", "White clover"))

# SIMPLE APPROACH: Create summary statistics manually without dplyr grouping
treatment_groups <- unique(alpha_data$TreatmentGroup)
plants <- unique(alpha_data$Plant)

# Create plant ranking manually
plant_ranking <- data.frame()

for (treatment in treatment_groups) {
  for (plant in plants) {
    subset_data <- alpha_data[alpha_data$TreatmentGroup == treatment & alpha_data$Plant == plant, ]
    
    if (nrow(subset_data) > 0) {
      mean_shannon <- mean(subset_data$Shannon, na.rm = TRUE)
      n_samples <- nrow(subset_data)
      
      plant_ranking <- rbind(plant_ranking, data.frame(
        TreatmentGroup = treatment,
        Plant = plant,
        mean_shannon = round(mean_shannon, 3),
        n_samples = n_samples
      ))
    }
  }
}

# Calculate ranks manually
final_ranking <- data.frame()
for (treatment in treatment_groups) {
  treatment_data <- plant_ranking[plant_ranking$TreatmentGroup == treatment, ]
  treatment_data$rank <- rank(-treatment_data$mean_shannon)  # Rank from highest to lowest
  final_ranking <- rbind(final_ranking, treatment_data)
}

final_ranking <- final_ranking[order(final_ranking$TreatmentGroup, final_ranking$rank), ]

print("Plant Diversity Rankings Within Each Treatment:")
print(final_ranking)

# Summary statistics
summary_stats <- data.frame()

for (treatment in treatment_groups) {
  for (plant in plants) {
    subset_data <- alpha_data[alpha_data$TreatmentGroup == treatment & alpha_data$Plant == plant, ]
    
    if (nrow(subset_data) > 0) {
      n_samples <- nrow(subset_data)
      mean_shannon <- round(mean(subset_data$Shannon, na.rm = TRUE), 3)
      sd_shannon <- round(sd(subset_data$Shannon, na.rm = TRUE), 3)
      se_shannon <- round(sd_shannon / sqrt(n_samples), 3)
      
      summary_stats <- rbind(summary_stats, data.frame(
        TreatmentGroup = treatment,
        Plant = plant,
        n_samples = n_samples,
        mean_shannon = mean_shannon,
        sd_shannon = sd_shannon,
        se_shannon = se_shannon
      ))
    }
  }
}

print("Shannon Diversity Summary Statistics - Plants within Treatments:")
print(summary_stats)

# Plot: Compare plants within each treatment group - FIVE PANELS IN ONE ROW
p_plants_in_treatments <- ggplot(alpha_data, aes(x = Plant, y = Shannon, fill = Plant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = Plant)) +
  stat_compare_means(method = "kruskal.test", label = "p.format", 
                     label.y = max(alpha_data$Shannon) + 0.2, size = 4) +
  facet_wrap(~ TreatmentGroup, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(
    x = "Plant Species",
    y = "Shannon Diversity Index",
    title = "Shannon Diversity: Comparison of Plant Species Within Each Treatment",
    subtitle = "Chloroplast and Mitochondria ASVs Removed"
  )

print(p_plants_in_treatments)

# Add pairwise comparisons between plants within each treatment - FIVE PANELS IN ONE ROW
p_plants_pairwise <- ggplot(alpha_data, aes(x = Plant, y = Shannon, fill = Plant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2, aes(color = Plant)) +
  stat_compare_means(comparisons = list(c("Buttercup", "Holcus"),
                                        c("Buttercup", "White clover"),
                                        c("Holcus", "White clover")),
                     method = "wilcox.test", 
                     label = "p.signif",
                     step.increase = 0.1,
                     size = 4) +
  facet_wrap(~ TreatmentGroup, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(
    x = "Plant Species",
    y = "Shannon Diversity Index",
    title = "Shannon Diversity with Pairwise Comparisons Between Plant Species",
    subtitle = "Within Each Treatment Group - Chloroplast and Mitochondria Removed"
  )

print(p_plants_pairwise)

# Additional visualization: Bar plot showing mean diversity - FIVE PANELS IN ONE ROW
p_mean_diversity <- ggplot(summary_stats, aes(x = Plant, y = mean_shannon, fill = Plant)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), 
                width = 0.2) +
  geom_text(aes(label = paste0("n=", n_samples)), vjust = -0.5, size = 2.5) +
  facet_wrap(~ TreatmentGroup, nrow = 1) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(
    x = "Plant Species",
    y = "Mean Shannon Diversity Index (± SE)",
    title = "Mean Shannon Diversity of Plant Species Within Each Treatment",
    subtitle = "Error bars represent standard error"
  )

print(p_mean_diversity)

# Statistical tests for each treatment group comparing plants
treatment_groups <- unique(alpha_data$TreatmentGroup)
plant_results <- list()

for (treatment in treatment_groups) {
  treatment_data <- alpha_data[alpha_data$TreatmentGroup == treatment, ]
  
  # Only run statistical tests if we have at least 2 plants with data
  plants_with_data <- unique(treatment_data$Plant)
  if (length(plants_with_data) >= 2) {
    # Kruskal-Wallis test for plants within this treatment
    kruskal_test <- kruskal.test(Shannon ~ Plant, data = treatment_data)
    
    # Pairwise Wilcoxon tests between plants
    pairwise_wilcox <- pairwise.wilcox.test(treatment_data$Shannon, 
                                            treatment_data$Plant,
                                            p.adjust.method = "BH")
    
    plant_results[[treatment]] <- list(
      kruskal_p = kruskal_test$p.value,
      pairwise = pairwise_wilcox$p.value
    )
    
    print(paste("Plant comparison results for", treatment, ":"))
    print(paste("Kruskal-Wallis p-value:", round(kruskal_test$p.value, 4)))
    print("Pairwise Wilcoxon p-values between plants:")
    print(pairwise_wilcox$p.value)
    print("---")
  } else {
    print(paste("Insufficient data for statistical tests in", treatment))
    print(paste("Plants with data:", paste(plants_with_data, collapse = ", ")))
  }
}

# Focus on the significant findings in Sterilized Optimized Clamp
significant_data <- alpha_data[alpha_data$TreatmentGroup == "Sterilized Optimized Clamp", ]

p_significant <- ggplot(significant_data, aes(x = Plant, y = Shannon, fill = Plant)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3, aes(color = Plant)) +
  stat_compare_means(comparisons = list(c("Buttercup", "White clover"),
                                        c("Holcus", "White clover")),
                     method = "wilcox.test", 
                     label = "p.signif",
                     step.increase = 0.1,
                     size = 5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold")
  ) +
  labs(
    x = "Plant Species",
    y = "Shannon Diversity Index",
    title = "Significant Plant Differences in Sterilized Optimized Clamp",
    subtitle = "White clover shows significantly different microbial diversity"
  )

print(p_significant)

# Detailed analysis of the significant finding
optimized_data <- alpha_data[alpha_data$TreatmentGroup == "Sterilized Optimized Clamp", ]

print("Detailed analysis of Sterilized Optimized Clamp:")
print(paste("Buttercup mean Shannon:", round(mean(optimized_data$Shannon[optimized_data$Plant == "Buttercup"]), 3)))
print(paste("Holcus mean Shannon:", round(mean(optimized_data$Shannon[optimized_data$Plant == "Holcus"]), 3)))
print(paste("White clover mean Shannon:", round(mean(optimized_data$Shannon[optimized_data$Plant == "White clover"]), 3)))

# Determine if White clover has higher or lower diversity
white_clover_mean <- mean(optimized_data$Shannon[optimized_data$Plant == "White clover"])
others_mean <- mean(optimized_data$Shannon[optimized_data$Plant != "White clover"])

if (white_clover_mean > others_mean) {
  direction <- "HIGHER"
} else {
  direction <- "LOWER"
}

print(paste("White clover has", direction, "microbial diversity than other plants"))
print(paste("White clover mean:", round(white_clover_mean, 3)))
print(paste("Other plants mean:", round(others_mean, 3)))

# Create a summary table of key findings
key_findings <- data.frame(
  Treatment = c("Sterilized Standard Clamp", "Sterilized Optimized Clamp", "Sterilized No Clamp",
                "Washed Clamp", "Washed No Clamp"),
  Kruskal_P_value = c(0.1767, 0.0047, 0.1338, 0.7066, NA),  # Will update with actual values
  Significant = c("No", "Yes", "No", "No", "TBD"),
  Pattern = c("No clear pattern", 
              paste("White clover has", tolower(direction), "diversity"), 
              "No clear pattern", "No clear pattern", "TBD")
)

print("Key Findings Summary:")
print(key_findings)

# Save all plots with adjusted dimensions for five panels in one row
ggsave("shannon_diversity_plants_in_treatments.png", p_plants_pairwise, 
       width = 18, height = 6, dpi = 300)
ggsave("shannon_diversity_mean_comparison.png", p_mean_diversity, 
       width = 18, height = 6, dpi = 300)
ggsave("significant_plant_differences.png", p_significant, 
       width = 8, height = 6, dpi = 300)

print("Analysis complete! All FIVE treatment groups analyzed in order:")
print("1. Sterilized Standard Clamp")
print("2. Sterilized Optimized Clamp") 
print("3. Sterilized No Clamp")
print("4. Washed Clamp")
print("5. Washed No Clamp")
