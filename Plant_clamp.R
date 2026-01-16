library(qiime2R)
library(phyloseq)
library(readxl)
library(rstudioapi)
library(tidyr)
library(dplyr)
library(ggplot2)
library(maps)
library(sf)
library(nlme)
library(plyr)
library(microViz)
library(patchwork)

setwd(dirname(getActiveDocumentContext()$path))

ASVs<-read_qza("feature-table-16S-Plant-all.qza")
names(ASVs)
#To access the actual data stored within the object, access the data as below:
ASVs$data[1:5,1:2] #show first 2 samples and first 5 ASVs

metadata<-read_q2metadata("metadata_Plant.txt")
head(metadata) # show top lines of metadata

taxonomy<-read_qza("taxonomy_SILVA-16S-Plant-all.qza")
head(taxonomy$data)
taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)

physeq<-qza_to_phyloseq(
  features="feature-table-16S-Plant-all.qza",
  tree="rooted-tree-16S-Plant-all.qza",
  "taxonomy_SILVA-16S-Plant-all.qza",
  metadata = "metadata_Plant.txt"
)
physeq

###import env data
save.image(file = "Phyloseq-Plant.RData")



# setting the seed to one value in order to created reproducible results
set.seed(1)  

# scaling the data to the smallest samples. Note: rngseed is similar to set.seed
sample_sums(physeq)
sort(sample_sums(physeq))
S16_scaled <- rarefy_even_depth(physeq, sample.size=47289, replace=FALSE, rngseed = 1) 

Meta<-sample_data(S16_scaled)


Rich_16S<-unlist(estimate_richness(S16_scaled, measures="Observed"))
Shannon_16S<-unlist(estimate_richness(S16_scaled, measures="Shannon"))

df1<-data.frame(Meta, Rich_16S,Shannon_16S,SampleID=rownames(Meta))



# Check the taxonomic table structure to find the rank containing 'Chloroplast'
tax_table(physeq)
tax_table(S16_scaled)
# Filter chloroplast sequences
chloroplast_otus <- subset_taxa(physeq, Order == "Chloroplast" | Family == "Chloroplast")
chloroplast_otus2 <- subset_taxa(S16_scaled, Order == "Chloroplast" | Family == "Chloroplast")
# View the chloroplast OTU abundance
otu_table(chloroplast_otus) 
otu_table(chloroplast_otus2)
# View abundance table of chloroplast sequences
chloroplast_abundance <- otu_table(chloroplast_otus)
chloroplast_abundance2 <- otu_table(chloroplast_otus2)
# Calculate relative abundance (if required)
relative_abundance <- prop.table(as(otu_table(chloroplast_otus), "matrix"), margin = 2)
relative_abundance2 <- prop.table(as(otu_table(chloroplast_otus2), "matrix"), margin = 2)

# Add the calculated abundances back to the original phyloseq object
physeq <- transform_sample_counts(physeq, function(x) x / sum(x))
physeq_chloroplast <- subset_taxa(physeq, Order == "Chloroplast" | Family == "Chloroplast")
S16_scaled <- transform_sample_counts(S16_scaled, function(x) x / sum(x))
S16_scaled_chloroplast <- subset_taxa(S16_scaled, Order == "Chloroplast" | Family == "Chloroplast")
# Summarize results
chloroplast_summary <- sample_sums(physeq_chloroplast)
chloroplast_summary2 <- sample_sums(S16_scaled_chloroplast)

# Get sample metadata
sample_data_df <- as.data.frame(sample_data(physeq))
sample_data_df2 <- as.data.frame(sample_data(S16_scaled))
# Calculate total chloroplast abundance per sample
chloroplast_abundance <- sample_sums(chloroplast_otus)
chloroplast_abundance2 <- sample_sums(chloroplast_otus2)
# Add chloroplast abundance to the sample metadata
sample_data_df$Chloroplast_Abundance <- chloroplast_abundance
sample_data_df2$Chloroplast_Abundance <- chloroplast_abundance2
colnames(sample_data_df2)
relative_abundance_Chl <-(sample_data_df2$Chloroplast_Abundance)/47289
sample_data_df2$relchl <- relative_abundance_Chl

ggplot(sample_data_df2, aes(x = Clamp, y = relchl, fill = sequencing)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Chloroplast Abundance Across Metadata Categories",
    x = "Clamp",
    y = "Chloroplast Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(sample_data_df2, aes(x = Clamp, y = relchl, fill = sequencing)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(0.75)) +
  geom_jitter(aes(color = Plant, shape = Treatment), 
              position = position_dodge(0.75),  # Remove `width`
              alpha = 0.5, size = 1.5) +
  theme_minimal() +
 # scale_y_log10() +  # Optional: log scale
  labs(
    title = "Chloroplast Abundance Across Metadata Categories",
    x = "Clamp",
    y = "Chloroplast Abundance",
    color = "Plant",
    shape = "Treatment"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

sample_data_df2 <- as.data.frame(sample_data_df2)
subset_data <- sample_data_df2[sample_data_df2$Plant %in% c("Buttercup", "Holcus", "White clover"), ]


ggplot(subset_data , aes(x = Clamp, y = relchl, fill = sequencing)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Chloroplast Abundance Across Metadata Categories",
    x = "Clamp",
    y = "Chloroplast Abundance"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#########################################


ggplot(subset_data, aes(x = Plant, y = Chloroplast_Abundance, color = Treatment, shape = sequencing)) +
  # Individual data points with jitter and dodge
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  
  # Mean values (black diamond) for each Plant-Treatment group
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 4, 
    color = "black",
    position = position_dodge(0.5)
  ) +
  
  # Error bars for mean confidence interval
  stat_summary(
    fun.data = mean_cl_normal,  
    geom = "errorbar", 
    width = 0.2,
    position = position_dodge(0.5)
  ) +
  
  # Theme and labels
  theme_minimal() +
  scale_y_log10() +  # Optional: log scale
  labs(
    title = "Chloroplast Abundance Across Individual Plants",
    x = "Plant",
    y = "Chloroplast Abundance",
    color = "Treatment",
    shape = "sequencing"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Facet by Clamp for separate visualization
  facet_wrap(~ Clamp)


#########################################
#Mitochondria
#######################################
plant_subset_physeq <- subset_samples(S16_scaled, Plant %in% c("Buttercup", "Holcus", "White clover"))

# Check the number of samples after subsetting
print(sample_names(plant_subset_physeq))

# Function to find which taxonomic rank contains "mitochondria"
find_mito_rank <- function(plant_subset_physeq) {
  tax_table_df <- as.data.frame(tax_table(plant_subset_physeq))  # Convert tax_table to data.frame
  
  # Find which columns contain "mitochondria"
  mito_rank <- colnames(tax_table_df)[apply(tax_table_df, 2, function(x) any(grepl("mitochondria", x, ignore.case = TRUE)))]
  
  if (length(mito_rank) == 0) {
    stop("Mitochondria not found in any taxonomic rank.")
  } else {
    return(mito_rank)
  }
}

# Find the correct taxonomic rank for mitochondria
mito_ranks <- find_mito_rank(plant_subset_physeq)  
print(paste("Mitochondria found in:", paste(mito_ranks, collapse = ", ")))

# Filter phyloseq object for mitochondrial sequences (use all found ranks)
mito_physeq <- subset_taxa(plant_subset_physeq, apply(tax_table(plant_subset_physeq), 1, function(row) any(grepl("mitochondria", row, ignore.case = TRUE))))

# Extract OTU/ASV table from the filtered phyloseq object
mito_abundance <- otu_table(mito_physeq)

# Summarize abundance per sample
mito_abundance_per_sample <- colSums(mito_abundance)

# Convert to a data frame for easy visualization
mito_abundance_df <- data.frame(
  Sample = names(mito_abundance_per_sample),
  Mitochondria_Abundance = mito_abundance_per_sample
)

# Merge with metadata
metadata <- as.data.frame(sample_data(plant_subset_physeq))
mito_abundance_df <- merge(mito_abundance_df, metadata, by.x = "Sample", by.y = "row.names")

# Calculate relative abundance of mitochondria
total_abundance_per_sample <- sample_sums(plant_subset_physeq)
mito_abundance_df$Relative_Abundance <- mito_abundance_per_sample / total_abundance_per_sample

# View the data
head(mito_abundance_df)


# Plot mitochondrial abundance across samples


ggplot(mito_abundance_df, aes(x = Plant, y = Mitochondria_Abundance, color = Treatment, shape = sequencing)) +
  # Individual data points with jitter and dodge
  geom_point(alpha = 0.5, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  
  # Mean values (black diamond) for each Plant-Treatment group
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 18, 
    size = 4, 
    color = "black",
    position = position_dodge(0.5)
  ) +
  
  # Error bars for mean confidence interval
  stat_summary(
    fun.data = mean_cl_normal,  
    geom = "errorbar", 
    width = 0.2,
    position = position_dodge(0.5)
  ) +
  
  # Theme and labels
  theme_minimal() +
  scale_y_log10() +  # Optional: log scale
  labs(
    title = "Mitochondrial Abundance Across Individual Plants",
    x = "Plant",
    y = "Mitochondria_Abundance",
    color = "Treatment",
    shape = "sequencing"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  
  # Facet by Clamp for separate visualization
  facet_wrap(~ Clamp)

###########################################
#Bacterial families abundance
########################################
# Ensure your taxonomic table has a Family level
if (!"Family" %in% rank_names(plant_subset_physeq)) {
  stop("No 'Family' rank found in the taxonomic table.")
}

# Aggregate OTU counts at the Family level
family_physeq <- tax_glom(plant_subset_physeq, taxrank = "Family")
Phylum_physeq <- tax_glom(plant_subset_physeq, taxrank = "Phylum")


# Transform to relative abundance if needed
family_physeq_rel <- transform_sample_counts(family_physeq, function(x) x / sum(x))
top20_families <- names(sort(taxa_sums(family_physeq_rel), decreasing = TRUE)[1:20])
family_physeq_top20 <- prune_taxa(top20_families, family_physeq_rel)

plant1_physeq <- subset_samples(family_physeq_top20, Plant == "Holcus")
plant2_physeq <- subset_samples(family_physeq_top20, Plant == "Buttercup")
plant3_physeq <- subset_samples(family_physeq_top20, Plant == "White clover")


# Function to create a plot for a given plant species
create_plot <- function(physeq, plant_name) {
  # Convert to data frame
  df <- psmelt(physeq)
  
  # Plot
  ggplot(df, aes(x = interaction(Treatment, Clamp), y = Abundance, fill = Family)) +
    geom_bar(stat = "identity", position = "stack") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          plot.title = element_text(size = 12, face = "bold")) +
    labs(
      x = "Treatment x Clamp",
      y = "Relative Abundance",
      title = paste("Top 20 Family-level Microbial Abundance for", plant_name)
    ) +
    scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
    scale_color_viridis_d() # Use a colorblind-friendly palette
}

# Create plots for each plant species
plot_plant1 <- create_plot(plant1_physeq, "Holcus")
plot_plant2 <- create_plot(plant2_physeq, "Buttercup")
plot_plant3 <- create_plot(plant3_physeq, "White clover")

# Display plots
plot_plant1
plot_plant2
plot_plant3


# Combine all plant species into one data frame
df_all <- psmelt(family_physeq_top20)

# Combine all plant species into one data frame
df_all <- psmelt(family_physeq_top20)

# Plot with faceting
ggplot(df_all, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "Relative Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_colour_viridis_d() # Use a colorblind-friendly palette

# Plot with faceting
ggplot(df_all, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Sample names",
    y = "Relative Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_colour_brewer(palette = "Set1") # Use a colorblind-friendly palette







# Summarize the number of reads per family
family_reads <- taxa_sums(family_physeq)
family_data <- data.frame(Family = names(family_reads), Reads = family_reads)

# Order the families by number of reads
family_data <- family_data[order(family_data$Reads, decreasing = TRUE),]

# Select the top 10 families
top10_families <- head(family_data, 10)

# Create the pie chart
ggplot(top10_families, aes(x = "", y = Reads, fill = Family)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y") +
  labs(title = "Top 10 Families by Number of Reads", x = "", y = "") +
  theme_void() +  # Remove the background and axes
  scale_fill_brewer(palette = "Paired")  # Nice color palette





# Check the taxonomy table
tax_table(Phylum_physeq)

# Check sample data
sample_data(Phylum_physeq)

# Check OTU/ASV table
otu_table(Phylum_physeq)

colnames(sample_data(family_physeq))


# Fix taxonomic labels for unknown families
family_physeq <- tax_fix(family_physeq, unknowns = c("Unknown_Family", "uncultured"))


# Check available treatments
sample_data(family_physeq)$Treatment  # Replace "Treatment" with the actual metadata column name

# Subset phyloseq object for a specific treatment 
Sterlized_physeq <- subset_samples(family_physeq, Treatment == "Sterlized")
# Subset for another treatment
Washed_physeq <- subset_samples(family_physeq, Treatment == "Washed")

# Verify that subsetting worked
sample_data(Sterlized_physeq)
sample_data(Washed_physeq)

# Check available treatments
sample_data(Sterlized_physeq)$sequencing 
# Subset phyloseq object for a specific treatment 
Sterlized_before_optimization_physeq <- subset_samples(Sterlized_physeq, sequencing == "before_optimization")
Sterlized_before_optimization_clamp_physeq <- subset_samples(Sterlized_before_optimization_physeq, Clamp == "clamp")


sample_data(Sterlized_before_optimization_clamp_physeq)$Plant <- 
  factor(sample_data(Sterlized_before_optimization_clamp_physeq)$Plant, 
         levels = sort(unique(sample_data(Sterlized_before_optimization_clamp_physeq)$Plant)))



p1 <- comp_barplot(Sterlized_before_optimization_clamp_physeq, tax_level = "Family", n_taxa = 30, group_by = "Plant")
p1

# Wrap the plots and adjust layout
p1 <- p1 %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect")

# Flip coordinates
p1 <- p1 & coord_flip()

# Adjust axis labels
p1 <- p1 & labs(x = NULL, y = NULL) & theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Print final plot
p1



##########################
#after_optimization
##########################

# Check available treatments
sample_data(Sterlized_physeq)$sequencing 
# Subset phyloseq object for a specific treatment 
Sterlized_after_optimization_physeq <- subset_samples(Sterlized_physeq, sequencing == "after_optimization")

sample_data(Sterlized_after_optimization_physeq)$Clamp 
# Subset phyloseq object for a specific treatment 
Sterlized_after_optimization_clamp_physeq <- subset_samples(Sterlized_after_optimization_physeq, Clamp == "clamp")
Sterlized_after_optimization_no_clamp_physeq <- subset_samples(Sterlized_after_optimization_physeq, Clamp == "no_clamp")
p3 <- comp_barplot(Sterlized_after_optimization_clamp_physeq, tax_level = "Family", n_taxa = 20, group_by = "Plant")
p3

# Wrap the plots and adjust layout
p3 <- p3 %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect")

# Flip coordinates
p3 <- p3 & coord_flip()

# Adjust axis labels
p3 <- p3 & labs(x = NULL, y = NULL) & theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Print final plot
p3

p4 <- comp_barplot(Sterlized_after_optimization_no_clamp_physeq, tax_level = "Family", n_taxa = 20, group_by = "Plant")
p4

# Wrap the plots and adjust layout
p4 <- p4 %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect")

# Flip coordinates
p4 <- p4 & coord_flip()

# Adjust axis labels
p4 <- p4 & labs(x = NULL, y = NULL) & theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Print final plot
p4


p_combined <- (p1|p3) + 
  plot_annotation(title = "Comparison of Clamping before vs. After optimization") & 
  theme(legend.position = "bottom")  # Move legend to the bottom

print(p_combined)




library(ggplot2)
library(patchwork)
library(forcats)  # For reordering factors

# Ensure 'Plant' is in alphabetical order in both phyloseq objects
sample_data(Sterlized_after_optimization_clamp_physeq)$Plant <- 
  factor(sample_data(Sterlized_after_optimization_clamp_physeq)$Plant, 
         levels = sort(unique(sample_data(Sterlized_after_optimization_clamp_physeq)$Plant)))

sample_data(Sterlized_after_optimization_no_clamp_physeq)$Plant <- 
  factor(sample_data(Sterlized_after_optimization_no_clamp_physeq)$Plant, 
         levels = sort(unique(sample_data(Sterlized_after_optimization_no_clamp_physeq)$Plant)))

# Generate bar plots
p3 <- comp_barplot(Sterlized_after_optimization_clamp_physeq, tax_level = "Family", 
                   n_taxa = 20, group_by = "Plant")
p3
p3 <- p3 %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect")

# Flip coordinates
p3 <- p3 & coord_flip()

# Adjust axis labels
p3 <- p3 & labs(x = NULL, y = NULL) & theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Print final plot
p3
p4 <- comp_barplot(Sterlized_after_optimization_no_clamp_physeq, tax_level = "Family", 
                   n_taxa = 20, group_by = "Plant")
p4
p4 <- p4 %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect")

# Flip coordinates
p4 <- p4 & coord_flip()

# Adjust axis labels
p4 <- p4 & labs(x = NULL, y = NULL) & theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Print final plot
p4


p3 <-p3 + theme(legend.position = "none")  # Remove legend from p3
p4 <- p4 + theme(legend.position = "none")  # Remove legend from p4
# Combine p3 and p4 with a shared legend
p_combined <- (p1 | p3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Display the final plot
print(p_combined)



