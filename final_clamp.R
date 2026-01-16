library(phyloseq)
library(qiime2R)
library(readxl)
library(rstudioapi)
library(tidyr)
library(ggplot2)
library(maps)
library(sf)
library(nlme)
library(plyr)
library(microViz)
library(vegan)  # For diversity calculations
library(ggplot2)  # For visualization
library(dplyr)  # For data manipulation
library(rstatix)  # For statistical tests
library(gridExtra)  # For arranging multiple plots
library(biomformat)
library(Biostrings)
library(pheatmap)
library(KEGGREST)
library(reshape2)
library(httr)
library(jsonlite)

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
save.image(file = "final-clamp.RData")


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
family_physeq <- tax_glom(plant_subset_physeq, taxrank = "Family")

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


family_physeq_rel <- transform_sample_counts(family_physeq, function(x) x / 40000)
top20_families <- names(sort(taxa_sums(family_physeq_rel), decreasing = TRUE)[1:20])
family_physeq_top20 <- prune_taxa(top20_families, family_physeq_rel)
sample_data(family_physeq_top20)$Clamp <- factor(sample_data(family_physeq_top20)$Clamp)
sample_data(family_physeq_top20)$Clamp <- factor(
  sample_data(family_physeq_top20)$Clamp,
  levels = c("clamp", "no_clamp")
)

# Combine all plant species into one data frame
df_all <- df_all %>%
  mutate(
    Treatment_agg = case_when(
      grepl("optimized", S_Name, ignore.case = TRUE) & grepl("sterilized", S_Name, ignore.case = TRUE) ~ "Clamp_Optimized_Sterilized",
      grepl("standard", S_Name, ignore.case = TRUE) & grepl("sterilized", S_Name, ignore.case = TRUE) ~ "Clamp_Standard_Sterilized",
      grepl("no.?clamp", S_Name, ignore.case = TRUE) & grepl("sterilized", S_Name, ignore.case = TRUE) ~ "NoClamp_Sterilized",
      grepl("clamp", S_Name, ignore.case = TRUE) & grepl("washed", S_Name, ignore.case = TRUE) ~ "Clamp_Washed",
      grepl("no.?clamp", S_Name, ignore.case = TRUE) & grepl("washed", S_Name, ignore.case = TRUE) ~ "NoClamp_Washed",
      TRUE ~ "Other"
    ),
    Treatment_agg = factor(Treatment_agg, levels = c(
      "Clamp_Optimized_Sterilized",
      "Clamp_Standard_Sterilized",
      "NoClamp_Sterilized",
      "Clamp_Washed",
      "NoClamp_Washed"
    ))
  )
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

# Create a consistent treatment column

# Plot with faceting
ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Treatments",
    y = "Relative Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +  scale_fill_manual(values = colorblind_palette) +
  #scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  +
  scale_x_discrete(labels = c(
    "Clamp_Optimized_Sterilized" = "Clamp Optimized (Sterilized)",
    "Clamp_Standard_Sterilized" = "Clamp Standard (Sterilized)",
    "NoClamp_Sterilized" = "No Clamp (Sterilized)",
    "Clamp_Washed" = "Clamp (Washed)",
    "NoClamp_Washed" = "No Clamp (Washed)"
  ))# Improve x-axis labels
# Use a colorblind-friendly palette

###########################################################################################
# --- Family aggregation and relative abundance ---
family_physeq <- tax_glom(plant_subset_physeq, taxrank = "Family")

family_physeq_rel <- transform_sample_counts(family_physeq, function(x) x / 40000)
top20_families <- names(sort(taxa_sums(family_physeq_rel), decreasing = TRUE)[1:20])
family_physeq_top20 <- prune_taxa(top20_families, family_physeq_rel)

# --- Combine all plant species into one data frame ---
df_all <- psmelt(family_physeq_top20)

# --- Create Treatment_agg column based on S_Name ---
df_all <- df_all %>%
  mutate(
    Treatment_agg = case_when(
      grepl("optimized", S_Name, ignore.case = TRUE) & grepl("sterilized", S_Name, ignore.case = TRUE) ~ "Clamp_Optimized_Sterilized",
      grepl("standard", S_Name, ignore.case = TRUE) & grepl("sterilized", S_Name, ignore.case = TRUE) ~ "Clamp_Standard_Sterilized",
      grepl("no.?clamp", S_Name, ignore.case = TRUE) & grepl("sterilized", S_Name, ignore.case = TRUE) ~ "NoClamp_Sterilized",
      grepl("clamp", S_Name, ignore.case = TRUE) & grepl("washed", S_Name, ignore.case = TRUE) ~ "Clamp_Washed",
      grepl("no.?clamp", S_Name, ignore.case = TRUE) & grepl("washed", S_Name, ignore.case = TRUE) ~ "NoClamp_Washed",
      TRUE ~ "Other"
    ),
    Treatment_agg = factor(Treatment_agg, levels = c(
      "Clamp_Optimized_Sterilized",
      "Clamp_Standard_Sterilized",
      "NoClamp_Sterilized",
      "Clamp_Washed",
      "NoClamp_Washed"
    ))
  )

# --- Color palette (20 colors) ---
colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#5D5D5D",
  "#66CCEE", "#D55E00", "#999999", "#CC79A7",
  "#0072B2", "#882255", "#E69F99", "#F0E442",
  "#B15928", "#FF7F00", "#6A3D9A",
  "#FFFF99", "#1F78B4", "#33A02C",
  "#B2DF8A", "#FB9A99"
)

# --- Plot ---
ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Treatments and Clamping",
    y = "Relative Abundance (%)\n16S rRNA Gene Sequences",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = c(
    "Clamp_Optimized_Sterilized" = "Clamp Optimized (Sterilized)",
    "Clamp_Standard_Sterilized" = "Clamp Standard (Sterilized)",
    "NoClamp_Sterilized" = "No Clamp (Sterilized)",
    "Clamp_Washed" = "Clamp (Washed)",
    "NoClamp_Washed" = "No Clamp (Washed)"
  ))

###########################################################################################

# Plot with faceting
ggplot(df_all, aes(x = S_Name , y = Abundance, fill = Family)) +
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
  #scale_fill_manual(values = colorblind_palette)+
  scale_x_discrete(labels = c(
    "Clamp_Optimized_Sterilized" = "Clamp Optimized (Sterilized)",
    "Clamp_Standard_Sterilized" = "Clamp Standard (Sterilized)",
    "NoClamp_Sterilized" = "No Clamp (Sterilized)",
    "Clamp_Washed" = "Clamp (Washed)",
    "NoClamp_Washed" = "No Clamp (Washed)"
  )) # Use a colorblind-friendly palette
########################Final Plot
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
##################
# Data preparation: calculate relative abundance per sample (S_Name)
df_agg2 <- df_agg %>%
  filter(Category %in% c("Chloroplast", "Mitochondria", "Bacteria")) %>%
  group_by(S_Name, Plant) %>%
  mutate(Abundance_relative = Abundance / sum(Abundance) * 100)  # Relative abundance within each sample (S_Name)

# Explicitly order the S_Name factor based on the numeric part, ensuring the order is consistent across all plants
df_agg2$S_Name <- factor(df_agg2$S_Name, 
                         levels = unique(df_agg2$S_Name[order(as.numeric(gsub("H", "", df_agg2$S_Name)))]))

# Plotting
ggplot(df_agg2, aes(x = S_Name, y = Abundance_relative, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +  # Use "stack" to show the categories stacked by their relative abundance
  geom_text(
    aes(label = sprintf("%.0f%%", Abundance_relative)),
    position = position_stack(vjust = 0.5),  # Position labels in the center of each stacked section
    size = 1.5, color = "black"
  ) +
  facet_wrap(~ Plant, scales = "free_x") +  # Faceting by plant
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
    x = "Samples, Clamping and Treatment",
    y = "Relative Abundance (%)\n16S rRNA Gene Sequences",
    fill = "Category"  # Updated to "Category" as it’s a grouping variable
  ) +
  scale_fill_manual(values = color_palette) +  # Use custom color palette
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Adjust x-axis labels if needed

##############










#########################relative abundance without chloroplast and mitochondria###########

# Transform to relative abundance if needed
family_physeq_without_chm <- subset_taxa(family_physeq, !grepl("Chloroplast", Family) & !grepl("Mitochondria", Family))
# Extract sample sums (total reads per sample)
Micro_abundance2 <- sample_sums(family_physeq_without_chm)

# Add Micro_abundance as a new column in sample metadata
sample_data(family_physeq_without_chm )$Micro_abundance2 <- Micro_abundance2

# Compute relative abundance
relative_abundance_micro2 <- Micro_abundance2 / 40000

# Add relative abundance to sample metadata
sample_data(family_physeq_without_chm)$relmicro2 <- relative_abundance_micro2

# Optional: Convert phyloseq sample metadata to a data frame
sample_data_df3 <- as.data.frame(sample_data(family_physeq_without_chm))


family_physeq_rel2 <- transform_sample_counts(family_physeq_without_chm , function(x) x / 40000)
top20_families2 <- names(sort(taxa_sums(family_physeq_rel2), decreasing = TRUE)[1:20])
family_physeq_top20_2 <- prune_taxa(top20_families2, family_physeq_rel2)
sample_data(family_physeq_top20_2)$Clamp <- factor(sample_data(family_physeq_top20_2)$Clamp)
sample_data(family_physeq_top20_2)$Clamp <- factor(
  sample_data(family_physeq_top20_2)$Clamp,
  levels = c("clamp", "no_clamp")
)

# Combine all plant species into one data frame
df_all2 <- psmelt(family_physeq_top20_2)
# Ensure Clamp is ordered correctly
df_all2$Clamp <- factor(df_all2$Clamp, levels = c("clamp", "no_clamp"))

# colour
colorblind_palette2 <- c("#E69F00", "#56B4E9", "#009E73", '#CAB2D6',
                        "#66CCEE", "#D55E00", "#999999","#CC79A7",  
                        '#8DD3C7', "#882255",'#E69F99', "#F0E442", 
                        '#B15928', '#FF7F00','#6A3D9A',
                        '#FFFF99','#1F78B4','#33A02C',
                        '#B2DF8A','#FB9A99' )




# Plot with faceting
ggplot(df_all2, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
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
  ) 
 +
  labs(
    x = "Treatment and Clamping",
    y = "Relative Abundance",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +  scale_fill_manual(values = colorblind_palette2) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels
# Use a colorblind-friendly palette

# Plot with faceting

ggplot(df_all2, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
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
    x = "Treatment and Clamping",
    y = "Relative Abundance",
    title = ""
  ) +
  scale_fill_manual(values = colorblind_palette2) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))

################### total read counts ###############

# Identify the top 20 most abundant Families
top_families <- names(sort(taxa_sums(family_physeq), decreasing = TRUE)[1:20])

# Prune the phyloseq object to keep only the top 20 Families
family_top_physeq <- prune_taxa(top_families, family_physeq)

#Convert to a melted dataframe for ggplot2
family_melt <- psmelt(family_top_physeq)

 #Plot absolute read counts for top 20 families

# Plot with faceting
ggplot(family_melt, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "Treatment and Clamping",
    y = "Read Count",
    title = " Family-level Microbial Abundance by Plant Species"
  ) +  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels
# Use a colorblind-friendly palette

# Plot with faceting
ggplot(family_melt, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "Sample names",
    y = "Read count",
    title = "Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_fill_manual(values = colorblind_palette) # Use a colorblind-friendly palette

##########################Read counts without chloroplast and mitochondria#################



# 2. Remove chloroplast and mitochondrial sequences
family_physeq_without_chm <- subset_taxa(family_physeq, !grepl("Chloroplast", Family) & !grepl("Mitochondria", Family))

# 3. Identify the top 20 most abundant remaining Families
top_families2 <- names(sort(taxa_sums(family_physeq_without_chm), decreasing = TRUE)[1:20])

# 4. Prune the phyloseq object to keep only these top 20 Families
family_physeq_without_chm<- prune_taxa(top_families2, family_physeq_without_chm)

# 5. Convert to a melted dataframe for ggplot2
family_melt2 <- psmelt(family_physeq_without_chm)

#Plot absolute read counts for top 20 families

# Plot with faceting
ggplot(family_melt2, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "Read Count",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +  scale_fill_manual(values = colorblind_palette2) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x))  # Improve x-axis labels
# Use a colorblind-friendly palette

# Plot with faceting
ggplot(family_melt2, aes(x = Sample , y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Sample names",
    y = "Read count",
    title = "Top 20 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\.", " x ", x)) + # Improve x-axis labels
  scale_fill_manual(values = colorblind_palette2) # Use a colorblind-friendly palette

##############shanon diversity############


# Compute Shannon diversity index
shannon_div <- estimate_richness(family_physeq, measures = "Shannon")
meta_data <- as.data.frame(sample_data(family_physeq))

# Combine diversity index with metadata
data <- cbind(meta_data, shannon_div)

# Perform Wilcoxon rank sum test
wilcox_test_result <- wilcox_test(Shannon ~ Plant, data = data)

# Determine significance level
sig_label <- case_when(
  wilcox_test_result$p < 0.001 ~ "***",
  wilcox_test_result$p < 0.01 ~ "**",
  wilcox_test_result$p < 0.05 ~ "*",
  TRUE ~ "ns"
)


# Boxplot with significance annotation
# Subset data for each plot condition
plot_data1 <- filter(data, Treatment == "Sterilized", sequencing == "Optimized") 
plot_data2 <- filter(data, Treatment == "Sterilized", sequencing == "Standard")
plot_data3 <- filter(data, Treatment == "Washed", sequencing == "Optimized")

# Function to create plots
plot_shannon <- function(df, title) {
  ggplot(df, aes(x = Plant, y = Shannon, fill = Clamp)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    theme_minimal() +
    labs(title = title, y = "Shannon Index", x = "Plant")
}

# Generate plots
plot1 <- plot_shannon(plot_data1, "Sterilized-Optimized") 
plot2 <- plot_shannon(plot_data2, "Sterilized-Standard")
plot3 <- plot_shannon(plot_data3, "Washed-Optimized")

# Display plots

grid.arrange(plot1, plot2, plot3, ncol = 3)
#######new shahnon index
# Load required libraries
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(patchwork)
library(cowplot)  # For get_legend and plot_grid

# -----------------------------
# 1. Calculate Shannon diversity
# -----------------------------
shannon_div <- estimate_richness(family_physeq, measures = "Shannon")
meta_data <- as.data.frame(sample_data(family_physeq))
data <- cbind(meta_data, shannon_div)

# Ensure Shannon is numeric
data$Shannon <- as.numeric(as.character(data$Shannon))

# Ensure grouping variables are factors
data$Plant <- as.factor(data$Plant)
data$Clamp <- as.factor(data$Clamp)
data$Treatment <- as.factor(data$Treatment)
data$sequencing <- as.factor(data$sequencing)

# -----------------------------
# 2. Subset function
# -----------------------------
subset_and_fix <- function(df, treatment, seq_method) {
  d <- filter(df, Treatment == treatment, sequencing == seq_method)
  d$Shannon <- as.numeric(as.character(d$Shannon))
  return(d)
}

plot_data1 <- subset_and_fix(data, "Sterilized", "Optimized")
plot_data2 <- subset_and_fix(data, "Sterilized", "Standard")
plot_data3 <- subset_and_fix(data, "Washed", "Optimized")

# -----------------------------
# 3. Generate pairwise comparisons
# -----------------------------
get_comparisons <- function(df) {
  groups <- levels(df$Plant)
  combn(groups, 2, simplify = FALSE)
}

# -----------------------------
# 4. Plotting function
# -----------------------------
create_shannon_plot <- function(df, title) {
  comps <- get_comparisons(df)
  
  y_max <- max(df$Shannon, na.rm = TRUE)
  y_min <- min(df$Shannon, na.rm = TRUE)
  step <- (y_max - y_min) * 0.15
  y_positions <- seq(from = y_max + step, by = step, length.out = length(comps))
  
  ggplot(df, aes(x = Plant, y = Shannon, fill = Clamp)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    stat_compare_means(
      comparisons = comps,
      method = "wilcox.test",
      label = "p.signif",
      hide.ns = FALSE,
      size = 5,
      tip.length = 0.01,
      bracket.size = 0.7,
      label.y = y_positions
    ) +
    labs(title = title, x = "Plant", y = "Shannon Index") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text = element_text(face = "bold", size = 12)
    )
}

# -----------------------------
# 5. Generate plots
# -----------------------------
p1 <- create_shannon_plot(plot_data1, "Sterilized - Optimized")
p2 <- create_shannon_plot(plot_data2, "Sterilized - Standard")
p3 <- create_shannon_plot(plot_data3, "Washed - Optimized")

# -----------------------------
# 6. Extract bold legend
# -----------------------------
legend <- get_legend(
  ggplot(data, aes(x = Plant, y = Shannon, fill = Clamp)) +
    geom_boxplot() +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(face = "bold", size = 12)
    )
)

# -----------------------------
# 7. Combine plots with shared legend
# -----------------------------
final_plot <- plot_grid(
  plot_grid(p1, p2, p3, ncol = 3, align = "v"),
  legend,
  ncol = 1,
  rel_heights = c(4, 0.5)
)

# -----------------------------
# 8. Display the final plot
# -----------------------------
print(final_plot)




  ##########################################
library(phyloseq)
library(ggplot2)
library(dplyr)
library(rstatix)
library(ggpubr)
library(patchwork)

# 1. Calculate Shannon diversity and combine with metadata
shannon_div <- estimate_richness(family_physeq, measures = "Shannon")
meta_data <- as.data.frame(sample_data(family_physeq))
data <- cbind(meta_data, shannon_div)

# Make sure Shannon is numeric
data$Shannon <- as.numeric(as.character(data$Shannon))

# Ensure factors
data$Plant <- as.factor(data$Plant)
data$Clamp <- as.factor(data$Clamp)
data$Treatment <- as.factor(data$Treatment)
data$sequencing <- as.factor(data$sequencing)

# 2. Subset function
subset_and_fix <- function(df, treatment, seq_method) {
  d <- filter(df, Treatment == treatment, sequencing == seq_method)
  d$Shannon <- as.numeric(as.character(d$Shannon))
  return(d)
}

plot_data1 <- subset_and_fix(data, "Sterilized", "Optimized")
plot_data2 <- subset_and_fix(data, "Sterilized", "Standard")
plot_data3 <- subset_and_fix(data, "Washed", "Optimized")

# 3. Function to create boxplots per Plant comparing Clamp groups
plot_plant_clamp <- function(df, treatment_seq_label) {
  
  plants <- levels(df$Plant)
  
  plot_list <- list()
  
  for (pl in plants) {
    df_sub <- df %>% filter(Plant == pl)
    
    # Wilcoxon test comparing Clamp groups for this Plant subset
    stat_test <- wilcox_test(df_sub, Shannon ~ Clamp) %>% add_significance()
    
    y_max <- max(df_sub$Shannon, na.rm = TRUE)
    
    p <- ggplot(df_sub, aes(x = Clamp, y = Shannon, fill = Clamp)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
      stat_pvalue_manual(stat_test, label = "p.signif", tip.length = 0.01, size = 5,
                         y.position = y_max + 0.2) +
      labs(title = paste(pl, treatment_seq_label, sep = " - "),
           x = "Clamp",
           y = "Shannon Index") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    
    plot_list[[pl]] <- p
  }
  
  combined <- patchwork::wrap_plots(plot_list, ncol = length(plot_list))
  return(combined)
}

# 4. Generate plots for each subset
p1 <- plot_plant_clamp(plot_data1, "Sterilized - Optimized")
p2 <- plot_plant_clamp(plot_data2, "Sterilized - Standard")
p3 <- plot_plant_clamp(plot_data3, "Washed - Optimized")

# 5. Create a plot to extract the legend
legend_plot <- ggplot(data, aes(x = Clamp, y = Shannon, fill = Clamp)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "bottom")

legend <- get_legend(legend_plot)

# 6. Combine the plots vertically and add shared legend at bottom
final_plot <- p1 / p2 / p3 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# 7. Print the final plot
print(final_plot)


#######################################


################Beta diversity#################
# Beta diversity calculation

# Compute Bray-Curtis Beta Diversity
beta_dist <- phyloseq::distance(family_physeq, method = "bray")

# Perform PCoA Ordination
ordination <- ordinate(family_physeq, method = "PCoA", distance = beta_dist)



# Convert sample_data to a dataframe explicitly
sample_metadata <- data.frame(sample_data(family_physeq)) %>%
  mutate(
    Treatment = as.character(Treatment),
    sequencing = as.character(sequencing),
    Plant = as.factor(Plant),  # Ensure Plant remains a factor
    Treatment_Sequencing = as.character(paste0(Treatment, "_", sequencing))
  )

# Check for missing values
print(table(sample_metadata$Treatment_Sequencing, useNA = "ifany"))
print(table(sample_metadata$Plant, useNA = "ifany"))

# Convert sample names into a column for grouping
sample_metadata$SampleID <- rownames(sample_metadata)

# Define line types
line_types <- c(
  "Sterilized_Optimized" = "solid",
  "Sterilized_Standard" = "dotted",
  "Washed_Optimized" = "dotdash"
)

# Ensure Treatment_Sequencing levels match line_types
sample_metadata$Treatment_Sequencing <- factor(sample_metadata$Treatment_Sequencing, levels = names(line_types))

# Reassign updated metadata back to phyloseq object
sample_data(family_physeq) <- sample_data(sample_metadata)

# Plot Beta diversity with convex hulls for Treatment & Sequencing within Plant
p_beta <- plot_ordination(family_physeq, ordination, color = "Plant", shape = "Clamp") +
  geom_point(aes(color = Plant, shape = Clamp), size = 3, alpha = 0.7) +
  stat_ellipse(aes(linetype = Treatment_Sequencing, color = Treatment_Sequencing, 
                   group = interaction(Plant, Treatment_Sequencing)), 
               type = "t", size = 1) +
  scale_linetype_manual(values = line_types) +
  theme_minimal() +
  labs(title = "Beta Diversity (PCoA - Bray-Curtis)", color = "Plant", linetype = "Treatment & Sequencing")

print(p_beta)

p_beta <- plot_ordination(family_physeq, ordination, color = "Plant", shape = "Clamp") +
  geom_point(aes(color = Plant, shape = Clamp), size = 3, alpha = 0.7) +
  stat_ellipse(aes(linetype = Treatment_Sequencing, group = interaction(Plant, Treatment_Sequencing), color = Plant), 
               type = "t", size = 1) +
  scale_linetype_manual(values = line_types) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    strip.text = element_text(size = 118, face = "bold")  # facet labels
  ) +
  labs(
    title = "",
    color = "Plant",
    linetype = "Treatment & Sequencing",
    shape = "Clamp"
  )
print(p_beta)
# Extract ordination scores correctly
ordination_scores <- as.data.frame(ordination$vectors[, 1:2])  # Extract first two axes
colnames(ordination_scores) <- c("PCoA1", "PCoA2")  # Rename for clarity
ordination_scores$SampleID <- rownames(ordination_scores)  # Ensure SampleID for merging

# Merge ordination scores with metadata
plot_data <- sample_metadata %>%
  left_join(ordination_scores, by = "SampleID")  # Ensure Plant is included

# Check if Plant exists after merging
print(colnames(plot_data))

plot_data <- plot_data %>%
  mutate(PCoA1 = as.numeric(PCoA1),
         PCoA2 = as.numeric(PCoA2))

table(plot_data$Plant, plot_data$Treatment_Sequencing)
plot_data <- plot_data %>%
  drop_na(PCoA1, PCoA2)

# Create convex hulls
 # Get convex hull
hull_data <- plot_data %>%
  group_by(Plant, Treatment_Sequencing) %>%
  filter(n() >= 3) %>%
  group_modify(~ .x[chull(.x$PCoA1, .x$PCoA2), ])

# Plot Beta diversity with convex hulls for Treatment & Sequencing within Plant
p_beta2 <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2)) +
  geom_point(aes(color = Plant, shape = Clamp), size = 3, alpha = 0.7) +
  
  # Draw convex hulls
  geom_polygon(data = hull_data, aes(x = PCoA1, y = PCoA2, 
                                     group = interaction(Plant, Treatment_Sequencing), 
                                     fill = Treatment_Sequencing), 
               alpha = 0.2, color = "grey") +
  
  scale_linetype_manual(values = line_types) +
  theme_minimal() +
  labs(title = "Beta Diversity (PCoA - Bray-Curtis)", 
       color = "Plant", fill = "Treatment & Sequencing")

print(p_beta2)

############################
# Install PICRUSt2 if not installed
# Convert phyloseq object to BIOM format

# Load your phyloseq object
ps <- plant_subset_physeq  

# Extract OTU table
otu_table <- as(otu_table(ps), "matrix")

# Convert to a BIOM format
biom_file <- "asv_table.biom"
write_biom(make_biom(data = otu_table), biom_file)

# Load sequences from exported FASTA file
asv_seqs <- readDNAStringSet("exported_rep_seqs/dna-sequences.fasta")

# Assign to the phyloseq object
ps <- merge_phyloseq(ps, asv_seqs)

# Verify sequences
refseq(ps)


# Convert to FASTA
writeXStringSet(refseq(ps), "rep_seqs.fasta")
##########################################
# Run PICRUSt2 externally, then re-import the results
###################pathway prediction#################
#Import MetaCyc Pathway Prediction
pathway_data <- read.table("Path_abun_unstrat.tsv", 
                           header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(pathway_data)


# Create a phyloseq object for pathways
pathway_physeq <- phyloseq(otu_table(as.matrix(pathway_data), taxa_are_rows = TRUE))
#####pathway visualization
pheatmap(pathway_data)

# Assume pathway_data is already loaded

# Ensure rownames are stored as a column before merging
if (!"Pathway_ID" %in% colnames(pathway_data)) {
  pathway_data$Pathway_ID <- rownames(pathway_data)
}

# Create a mapping table for MetaCyc pathway names
metacyc_map <- data.frame(
  Pathway_ID = c("PWYO-1298", "PWYO-1586", "PWY66-389", "PWYO-862", "PWY-5973", "PWY6282", "PWY-1042",  
                 "PHOASLIPSYN-PWY", "PWY4FS-7", "PWY4FS-8", "SER-GLYSYN-PWY", "FASYN-INITAL-PWY",
                 "PWY-5690", "PWY-5667", "PWY01319", "NONOXIPENT-PWY", "PWY-8178",  
                 "VALSYN-PWY", "ILEUSYN-PWY", "BRANCHED-CHAIN-AA-SYN-PWY",   
                 "PWY0-1298", "PWY-6185", "PWY-7527", "PCPDEG-PWY", "FAO-PWY", "PWY-922"),
  
  Function = c("Pyrimidine Deoxyribonucleotides Salvage","Unknown", "Phytol Degradation", "Unknown", "Adenosylcobalamin Salvage",  
               "L-Histidine Biosynthesis", "Pyrimidine Ribonucleosides Salvage",  
               "Phosphatidylserine & Phosphatidylethanolamine Biosynthesis",
               "Fatty Acid Biosynthesis (Initiation)", "Fatty Acid Biosynthesis (Elongation)",
               "L-Serine and L-Glycine Biosynthesis", "Fatty Acid Biosynthesis (Initial Reactions)",
               "Glutathione Biosynthesis", "UDP-N-acetylmuramoyl-pentapeptide Biosynthesis",
               "Lysine Biosynthesis via DAP Pathway", "Non-Oxidative Phase of Pentose Phosphate Pathway",
               "Superpathway of Sulfur Amino Acid Biosynthesis", "L-Valine Biosynthesis",
               "L-Isoleucine Biosynthesis", "Branched-Chain Amino Acid Biosynthesis",
               "Pyrimidine Deoxyribonucleotides Salvage", "Unknown", "Nucleotide Salvage Pathway",
               "Pentachlorophenol Degradation", "Fatty Acid Beta-Oxidation",
               "Superpathway of Glutamate Biosynthesis")
)

# Convert Pathway_ID columns to character type for consistent merging
pathway_data$Pathway_ID <- as.character(pathway_data$Pathway_ID)
metacyc_map$Pathway_ID <- as.character(metacyc_map$Pathway_ID)

# Merge pathway mapping with pathway abundance data
pathway_data <- merge(pathway_data, metacyc_map, by = "Pathway_ID", all.x = TRUE)

# Replace missing function names with the original pathway IDs
pathway_data$Function[is.na(pathway_data$Function)] <- pathway_data$Pathway_ID

# Set function names as rownames for visualization
rownames(pathway_data) <- pathway_data$Function
pathway_data$Pathway_ID <- NULL  # Remove the Pathway_ID column

# Keep only numeric columns (abundance data)
numeric_pathway_data <- pathway_data[, sapply(pathway_data, is.numeric)]

# Select the top 20 most abundant pathways
top_n <- min(20, nrow(numeric_pathway_data))  # Ensures no error if <20 pathways exist
top_pathways <- numeric_pathway_data[order(rowSums(numeric_pathway_data), decreasing = TRUE), ][1:top_n, ]

# Plot heatmap with pathway names
heatmap_plot <- pheatmap(top_pathways, fontsize_row = 8, fontsize_col = 8, main = "Top 20 Functional Pathways")

####################################################


# Import KEGG Ortholog (KO) Predictions
ko_data <- read.table("pred_metagenome_unstrat_KO.tsv", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Check structure
head(ko_data)

# Convert to Phyloseq object
ko_physeq <- phyloseq(otu_table(as.matrix(ko_data), taxa_are_rows = TRUE))

# Select Top 20 KO Functions
ko_top <- ko_data[order(rowSums(ko_data), decreasing = TRUE), ][1:20, ]
ko_top$KO <- rownames(ko_top)  # Ensure KO IDs are stored

# Reshape for ggplot
ko_melt <- melt(ko_top, id.vars = "KO", variable.name = "Sample", value.name = "Abundance")

# Plot KEGG Ortholog Abundance
ggplot(ko_melt, aes(x = Sample, y = Abundance, fill = KO)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "Sample", y = "Abundance", fill = "KEGG Orthologs", title = "Top 20 Predicted Functions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



############################################

# Import KEGG Ortholog (KO) Predictions
ko_data <- read.table("pred_metagenome_unstrat_KO.tsv", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Convert to Phyloseq object
ko_physeq <- phyloseq(otu_table(as.matrix(ko_data), taxa_are_rows = TRUE))

# Select Top 10 KO Functions
ko_top <- ko_data[order(rowSums(ko_data), decreasing = TRUE), ][1:20, ]
ko_top$KO <- rownames(ko_top)  # Ensure KO IDs are stored

# Reshape data for ggplot
ko_melt <- melt(ko_top, id.vars = "KO", variable.name = "Sample", value.name = "Abundance")

# Extract Sample Metadata from Phyloseq Object
sample_metadata <- data.frame(sample_data(family_physeq))

# Ensure row names match sample IDs
sample_metadata$Sample <- rownames(sample_metadata)

# Merge KO data with sample metadata
ko_melt <- merge(ko_melt, sample_metadata, by = "Sample", all.x = TRUE)

# Convert metadata columns to factors
ko_melt$Treatment <- as.factor(ko_melt$Treatment)
ko_melt$Clamp <- as.factor(ko_melt$Clamp)
ko_melt$sequencing <- as.factor(ko_melt$sequencing)
ko_melt$Plant <- as.factor(ko_melt$Plant)


# Generate Faceted KO Abundance Plot
ggplot(ko_melt, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = KO)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "Abundance",
    title = "Top 20 Predicted KO Functions"
  ) +
  scale_fill_manual(values = colorblind_palette2) +
  scale_x_discrete(labels = function(x) gsub("_", " x ", x))  # Improve x-axis labels



# Import KEGG Ortholog (KO) Predictions
ko_data <- read.table("pred_metagenome_unstrat_KO.tsv", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Convert to Phyloseq object
ko_physeq <- phyloseq(otu_table(as.matrix(ko_data), taxa_are_rows = TRUE))

# Select Top 10 KO Functions
ko_top <- ko_data[order(rowSums(ko_data), decreasing = TRUE), ][1:20, ]
ko_top$KO <- rownames(ko_top)  # Ensure KO IDs are stored

# Reshape data for ggplot
ko_melt <- melt(ko_top, id.vars = "KO", variable.name = "Sample", value.name = "Abundance")

# Extract Sample Metadata from Phyloseq Object
sample_metadata <- data.frame(sample_data(family_physeq))

# Ensure row names match sample IDs
sample_metadata$Sample <- rownames(sample_metadata)

# Merge KO data with sample metadata
ko_melt <- merge(ko_melt, sample_metadata, by = "Sample", all.x = TRUE)

# Convert metadata columns to factors
ko_melt$Treatment <- as.factor(ko_melt$Treatment)
ko_melt$Clamp <- as.factor(ko_melt$Clamp)
ko_melt$sequencing <- as.factor(ko_melt$sequencing)
ko_melt$Plant <- as.factor(ko_melt$Plant)

# Extract unique KO identifiers
ko_ids <- unique(ko_melt$KO)

# Function to Query KEGG with Error Handling
fetch_kegg_name <- function(ko_id) {
  result <- tryCatch({
    kegg_info <- keggGet(ko_id)
    if (!is.null(kegg_info) && !is.na(kegg_info[[1]]$NAME)) {
      return(kegg_info[[1]]$NAME[1])  # Use the first name
    } else {
      return(NA)
    }
  }, error = function(e) NA)
  return(result)
}

# Query KEGG for KO Function Names (Batch Query)
ko_names <- sapply(ko_ids, fetch_kegg_name)

# Create KO Mapping Table
ko_mapping <- data.frame(KO = ko_ids, Function = ko_names)

# Merge KO names into the melted dataset
ko_melt1 <- merge(ko_melt, ko_mapping, by = "KO", all.x = TRUE)


# Generate Faceted KO Abundance Plot with Function Names
ggplot(ko_melt1, aes(x = interaction(Treatment, Clamp, sequencing), y = Abundance, fill = Function)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 12, face = "bold")) +
  labs(
    x = "Treatment x Clamp x Sequencing",
    y = "Abundance",
    title = "Top Predicted KEGG Functions"
  ) +
  scale_fill_manual(values = colorblind_palette2) +
  scale_x_discrete(labels = function(x) gsub("_", " x ", x))  # Improve x-axis labels



##############
#Import Enzyme Commission (EC) Predictions
ec_data <- read.table("pred_metagenome_unstrat.tsv", 
                      header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
head(ec_data)


#########################
#  FAPROTAXR

# Extract taxonomy table
taxonomy_table <- as.data.frame(tax_table(plant_subset_physeq))

# Save to CSV
write.csv(taxonomy_table, "taxonomy_table_methodpaper.csv", row.names = TRUE)
write.table(taxonomy_table, "taxonomy_table_methodpaper.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

########python3 collapse_table.py -i taxonomy_table_methodpaper.tsv -g FAPROTAX.txt -o faprotax_output_methodpaper.csv


faprotax_results <- read.csv("faprotax_output_methodpaper.csv", row.names = 1)
head(faprotax_results)



