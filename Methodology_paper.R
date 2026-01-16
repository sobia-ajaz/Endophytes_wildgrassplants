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
library(pairwiseAdonis)
library(vegan)
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
save.image(file = "Methodology-paper.RData")


# setting the seed to one value in order to created reproducible results
set.seed(1)  

S16_scaled <- rarefy_even_depth(physeq, sample.size=40000, replace=FALSE, rngseed = 1) 



Meta<-sample_data(S16_scaled)

Rich_16S<-unlist(estimate_richness(S16_scaled, measures="Observed"))
Shannon_16S<-unlist(estimate_richness(S16_scaled, measures="Shannon"))

df1<-data.frame(Meta, Rich_16S,Shannon_16S,SampleID=rownames(Meta))

# Check the taxonomic table structure to find the rank 
tax_table(S16_scaled)

#subset
plant_subset_physeq <- subset_samples(S16_scaled, Plant %in% c("Buttercup", "Holcus", "White clover"))

sample_sums(plant_subset_physeq)  # Lists total reads per sample
summary(sample_sums(plant_subset_physeq))

# Check the number of samples after subsetting
print(sample_names(plant_subset_physeq))
# Ensure your taxonomic table has a Family level
if (!"Family" %in% rank_names(plant_subset_physeq)) {
  stop("No 'Family' rank found in the taxonomic table.")
}
#set.seed(123) # Set a seed for reproducibility
#plant_subset_physeq_rar <- rarefy_even_depth(plant_subset_physeq, sample.size = min(sample_sums(plant_subset_physeq)))



##############################



# Aggregate OTU counts at the Family level
family_physeq <- tax_glom(plant_subset_physeq_rar, taxrank = "Family")

# Transform to relative abundance if needed

# Extract sample sums (total reads per sample)
Micro_abundance <- sample_sums(family_physeq)

# Add Micro_abundance as a new column in sample metadata
sample_data(family_physeq)$Micro_abundance <- Micro_abundance

# Compute relative abundance
relative_abundance_micro <- Micro_abundance / 40000

# Add relative abundance to sample metadata
sample_data(family_physeq)$relmicro <- relative_abundance_micro

# Optional: Convert phyloseq sample metadata to a data frame
sample_data_df2 <- as.data.frame(sample_data(family_physeq))



physeq_tss <- transform_sample_counts(family_physeq, function(x) x / sum(x) * 1e6)  # Converts to per million
family_physeq_rel <- transform_sample_counts(family_physeq, function(x) x / 40000)
top20_families <- names(sort(taxa_sums(family_physeq_rel), decreasing = TRUE)[1:20])
family_physeq_top20 <- prune_taxa(top20_families, family_physeq_rel)
sample_data(family_physeq_top20)$Clamp <- factor(sample_data(family_physeq_top20)$Clamp)
sample_data(family_physeq_top20)$Clamp <- factor(
  sample_data(family_physeq_top20)$Clamp,
  levels = c("clamp", "no_clamp")
)

# Combine all plant species into one data frame
df_all <- psmelt(family_physeq_top20)
# Ensure Clamp is ordered correctly
df_all$Clamp <- factor(df_all$Clamp, levels = c("clamp", "no_clamp"))

# colour
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73",'#5D5D5D', 
                        "#66CCEE", "#D55E00", "#999999","#CC79A7",  
                        "#0072B2", "#882255",'#E69F99', "#F0E442", 
                        '#B15928', '#FF7F00','#6A3D9A',
                        '#FFFF99','#1F78B4','#33A02C',
                        '#B2DF8A','#FB9A99')


# Plot with faceting
ggplot(df_all, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
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
  ) +  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels
# Use a colorblind-friendly palette

# Plot with faceting
ggplot(df_all, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
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
  scale_fill_manual(values = colorblind_palette) # Use a colorblind-friendly palette

################### comulative scale normalization ###############

library(metagenomeSeq)


# Function to perform CSS normalization on a phyloseq object
css_normalize_phyloseq <- function(physeq_obj) {
  # Convert to metagenomeSeq object
  physeq_css <- phyloseq_to_metagenomeSeq(physeq_obj)
  
  # Determine the cumulative normalization quantile
  p <- cumNormStat(physeq_css)
  
  # Apply cumulative sum scaling normalization
  physeq_css <- cumNorm(physeq_css, p = p)
  
  # Extract normalized counts
  norm_counts <- MRcounts(physeq_css, norm = TRUE)
  
  # Ensure taxa_are_rows matches the original format
  otu_table(physeq_obj) <- otu_table(norm_counts, taxa_are_rows = taxa_are_rows(physeq_obj))
  
  return(physeq_obj)
}

# Apply CSS normalization
physeq_normalized <- css_normalize_phyloseq(family_physeq)

# Extract top 20 most abundant families
top20_families_nor <- names(sort(taxa_sums(physeq_normalized), decreasing = TRUE)[1:20])

# Subset phyloseq object to keep only the top 20 families
family_physeq_top20_nor <- prune_taxa(top20_families_nor, physeq_normalized)

# Ensure 'Clamp' is a factor with correct levels
sample_data(family_physeq_top20_nor)$Clamp <- factor(sample_data(family_physeq_top20_nor)$Clamp, levels = c("clamp", "no_clamp"))

# Convert phyloseq object to data frame
df_all_nor <- psmelt(family_physeq_top20_nor)

# Ensure 'Clamp' is still correctly ordered
df_all_nor$Clamp <- factor(df_all_nor$Clamp, levels = c("clamp", "no_clamp"))

# Define colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#5D5D5D", 
                        "#66CCEE", "#D55E00", "#999999", "#CC79A7",  
                        "#0072B2", "#882255", "#E69F99", "#F0E442", 
                        "#B15928", "#FF7F00", "#6A3D9A", "#FFFF99", 
                        "#1F78B4", "#33A02C", "#B2DF8A", "#FB9A99")

# Create stacked bar plot with faceting by Plant species
ggplot(df_all_nor, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "CSS-Normalized Abundance",
    title = "Top 20 Family-Level Microbial Abundance by Plant Species"
  ) +  
  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels

# Plot with faceting
ggplot(df_all_nor, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Sample names",
    y = "CSS-Normalized Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_fill_manual(values = colorblind_palette) # Use a colorblind-friendly palette

####################Variance Stabilizing Transformation (VST) / Regularized Log (rlog)

library(DESeq2)


library(phyloseq)
library(DESeq2)
library(ggplot2)

# Step 1: Convert phyloseq object to DESeq2
dds <- phyloseq_to_deseq2(family_physeq, ~ 1)  # "~ 1" since no experimental design is specified

# Step 2: Apply Variance Stabilizing Transformation (VST)
dds_vst <- varianceStabilizingTransformation(dds)

# Step 3: Convert back to phyloseq object with transformed counts
physeq_vst <- phyloseq(otu_table(assay(dds_vst), taxa_are_rows = TRUE), 
                       tax_table(family_physeq), 
                       sample_data(family_physeq))

# Step 4: Extract top 20 most abundant families
top20_families_vst <- names(sort(taxa_sums(physeq_vst), decreasing = TRUE)[1:20])

# Step 5: Subset phyloseq object to keep only the top 20 families
physeq_top20_vst <- prune_taxa(top20_families_vst, physeq_vst)

# Step 6: Ensure 'Clamp' is a factor with correct levels
sample_data(physeq_top20_vst)$Clamp <- factor(sample_data(physeq_top20_vst)$Clamp, levels = c("clamp", "no_clamp"))

# Step 7: Melt phyloseq object into a data frame for ggplot
df_all_vst <- psmelt(physeq_top20_vst)

# Step 8: Ensure 'Clamp' is still correctly ordered
df_all_vst$Clamp <- factor(df_all_vst$Clamp, levels = c("clamp", "no_clamp"))

# Step 9: Define a colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#5D5D5D", 
                        "#66CCEE", "#D55E00", "#999999", "#CC79A7",  
                        "#0072B2", "#882255", "#E69F99", "#F0E442", 
                        "#B15928", "#FF7F00", "#6A3D9A", "#FFFF99", 
                        "#1F78B4", "#33A02C", "#B2DF8A", "#FB9A99")

# Step 10: Create stacked bar plot with faceting by Plant species
ggplot(df_all_vst, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "Variance Stabilized Abundance",
    title = "Top 20 Family-Level Microbial Abundance by Plant Species"
  ) +  
  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels


# Plot with faceting
ggplot(df_all_vst, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Sample names",
    y = "Variance Stabilized Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_fill_manual(values = colorblind_palette) # Use a colorblind-friendly palette
#########################total sum scaling#################

# Function to perform TSS normalization on a phyloseq object
tss_normalize_phyloseq <- function(physeq_obj) {
  transform_sample_counts(physeq_obj, function(x) x / sum(x))
}

# Apply TSS normalization
family_physeq_tss <- tss_normalize_phyloseq(family_physeq)

# Extract top 20 most abundant families
top20_families_tss <- names(sort(taxa_sums(family_physeq_tss), decreasing = TRUE)[1:20])

# Subset phyloseq object to keep only the top 20 families
family_physeq_top20_tss <- prune_taxa(top20_families_tss, family_physeq_tss)

# Ensure 'Clamp' is a factor with correct levels
sample_data(family_physeq_top20_tss)$Clamp <- factor(sample_data(family_physeq_top20_tss)$Clamp, levels = c("clamp", "no_clamp"))

# Convert phyloseq object to data frame
df_all_tss <- psmelt(family_physeq_top20_tss)

# Ensure 'Clamp' is still correctly ordered
df_all_tss$Clamp <- factor(df_all_tss$Clamp, levels = c("clamp", "no_clamp"))

# Define colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#5D5D5D", 
                        "#66CCEE", "#D55E00", "#999999", "#CC79A7",  
                        "#0072B2", "#882255", "#E69F99", "#F0E442", 
                        "#B15928", "#FF7F00", "#6A3D9A", "#FFFF99", 
                        "#1F78B4", "#33A02C", "#B2DF8A", "#FB9A99")

# Create stacked bar plot with faceting by Plant species
ggplot(df_all_tss, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme_minimal() +  # Clean theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "TSS-Normalized Abundance",
    title = "Top 20 Family-Level Microbial Abundance by Plant Species"
  ) +  
  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels


# Plot with faceting
ggplot(df_all_tss, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Sample names",
    y = "TSS-Normalized Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_fill_manual(values = colorblind_palette) # Use a colorblind-friendly palette



#################
#beta diversity
################
# Load required libraries
library(phyloseq)
library(vegan)          # For PERMANOVA
library(ggplot2)        # For visualization
library(pairwiseAdonis) # For pairwise comparisons
library(ape)            # For UniFrac

# Aggregate taxa at the Family level
family_physeq <- tax_glom(plant_subset_physeq_rar, taxrank = "Family")

# Compute Beta Diversity distances
bray_dist <- phyloseq::distance(family_physeq, method = "bray")
wunifrac_dist <- phyloseq::distance(family_physeq, method = "wunifrac")
unifrac_dist <- phyloseq::distance(family_physeq, method = "unifrac")

# Perform PCoA
pcoa_bray <- ordinate(family_physeq, method = "PCoA", distance = bray_dist)
pcoa_wunifrac <- ordinate(family_physeq, method = "PCoA", distance = wunifrac_dist)
pcoa_unifrac <- ordinate(family_physeq, method = "PCoA", distance = unifrac_dist)

# Extract metadata
metadata <- as.data.frame(sample_data(family_physeq))

# Function for PCoA plotting
plot_pcoa <- function(pcoa_result, title, color_var) {
  plot_ordination(family_physeq, pcoa_result, color = color_var) +
    geom_point(size = 4) +
    theme_minimal() +
    ggtitle(title)
}

# Generate PCoA plots
p1 <- plot_pcoa(pcoa_bray, "PCoA - Bray-Curtis (Family Level)", "Plant")
p2 <- plot_pcoa(pcoa_wunifrac, "PCoA - Weighted UniFrac (Family Level)", "Plant")
p3 <- plot_pcoa(pcoa_unifrac, "PCoA - Unweighted UniFrac (Family Level)", "Plant")

# Print plots
print(p1)
print(p2)
print(p3)

# Run PERMANOVA to test group differences
adonis_result <- adonis2(bray_dist ~ Plant, data = metadata, permutations = 999)
print(adonis_result)

rownames(metadata) <- sample_names(family_physeq)
# Convert sample_data to a standard data frame
metadata <- data.frame(as.matrix(sample_data(family_physeq)))

# Ensure 'Plant' is a factor
metadata$Plant <- as.factor(metadata$Plant)

# Set row names to match the distance matrix
rownames(metadata) <- sample_names(family_physeq)


# Pairwise PERMANOVA if 'Plant' has more than 2 groups
if (length(unique(metadata$Plant)) > 2) {
  library(pairwiseAdonis)  # Ensure pairwise.adonis2() is available
  pairwise_results <- pairwise.adonis2(bray_dist ~ Plant, data = metadata, permutations = 999)
  print(pairwise_results)
}

# Check for homogeneity of dispersion (betadisper)
if (length(unique(metadata$Plant)) > 1) {  # betadisper needs >1 group
  betadisper_test <- betadisper(bray_dist, metadata$Plant)
  anova_dispersion <- anova(betadisper_test)
  print(anova_dispersion)
  
  # If dispersion is significant, run Tukey HSD post-hoc test
  if (anova_dispersion$`Pr(>F)`[1] < 0.05) {
    print(TukeyHSD(betadisper_test))
  }
}



#####################FILTERING###############################
# Function to filter out unclassified taxa at multiple taxonomic levels
filter_unclassified_taxa <- function(plant_subset_physeq_rar, tax_levels = c("Phylum", "Class", "Order", "Family", "Genus")) {
  # Extract taxonomy table
  tax_table_df <- as.data.frame(tax_table(plant_subset_physeq_rar))
  
  # Define keywords that indicate unclassified taxa
  unclassified_patterns <- c("unclassified", "uncultured", "Incertae Sedis", "Ambiguous_taxa", "Unknown", "NA", "")
  
  # Identify rows where any specified taxonomic level contains unclassified patterns
  taxa_to_keep <- apply(tax_table_df[, tax_levels, drop = FALSE], 1, function(row) {
    !any(sapply(unclassified_patterns, function(pattern) any(grepl(pattern, row, ignore.case = TRUE))))
  })
  
  # Filter the taxonomy table
  filtered_taxa <- tax_table_df[taxa_to_keep, ]
  
  # Subset phyloseq object
  filtered_physeq <- prune_taxa(rownames(filtered_taxa), plant_subset_physeq_rar)
  
  return(filtered_physeq)
}

#  usage:
# Assume `plant_subset_physeq_rar` is  phyloseq object
physeq_filtered <- filter_unclassified_taxa(plant_subset_physeq_rar, tax_levels = c("Phylum", "Class", "Order", "Family", "Genus"))

# Check the number of taxa before and after filtering
cat("Number of taxa before filtering:", ntaxa(plant_subset_physeq_rar), "\n")
cat("Number of taxa after filtering:", ntaxa(physeq_filtered), "\n")