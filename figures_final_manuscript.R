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

library(phyloseq)
library(ggplot2)
library(dplyr)

# Transform counts to relative abundance
family_physeq_rel <- transform_sample_counts(family_physeq, function(x) x / 40000)

# Select top 40 families
top40_families <- names(sort(taxa_sums(family_physeq_rel), decreasing = TRUE)[1:40])
family_physeq_top40 <- prune_taxa(top40_families, family_physeq_rel)

# Order 'Clamp' factor
sample_data(family_physeq_top40)$Clamp <- factor(
  sample_data(family_physeq_top40)$Clamp,
  levels = c("clamp", "no_clamp")
)

# Convert to data frame
df_all <- psmelt(family_physeq_top40)

# Define the treatment aggregation categories
df_all <- df_all %>%
  mutate(Treatment_agg = case_when(
    grepl("clamp_optimized_sterlized", S_Name, ignore.case = TRUE) ~ "Clamp_Optimized_Sterilized",
    grepl("clamp_standard_sterlized", S_Name, ignore.case = TRUE) ~ "Clamp_Standard_Sterilized",
    grepl("no_clamp_sterlized", S_Name, ignore.case = TRUE) ~ "NoClamp_Sterilized",
    grepl("clamp_washed", S_Name, ignore.case = TRUE) ~ "Clamp_Washed",
    grepl("no_clamp_washed", S_Name, ignore.case = TRUE) ~ "NoClamp_Washed",
    TRUE ~ "Other"
  ))

# ✅ Reorder x-axis so all Sterilized come first, then Washed
df_all$Treatment_agg <- factor(
  df_all$Treatment_agg,
  levels = c(
    "Clamp_Optimized_Sterilized",
    "Clamp_Standard_Sterilized",
    "NoClamp_Sterilized",
    "Clamp_Washed",
    "NoClamp_Washed"
  )
)

colorblind_palette <- c(
  # Original 20
  "#E69F00", "#56B4E9", "#009E73", "#5D5D5D",
  "#66CCEE", "#D55E00", "#999999", "#CC79A7",
  "#0072B2", "#882255", "#E69F99", "#F0E442",
  "#B15928", "#FF7F00", "#6A3D9A",
  "#FFFF99", "#1F78B4", "#33A02C",
  "#B2DF8A", "#FB9A99",
  
  # Additional 20 new distinct colors
  "#A6761D", "#F781BF", "#BC80BD", "#CAB2D6",
  "#FDB462", "#80B1D3", "#B3DE69", "#FB8072",
  "#8DD3C7", "#CCEBC5", "#FFD92F", "#E41A1C",
  "#984EA3", "#4DAF4A", "#FF33CC", "#A6CEE3",
  "#1B9E77", "#D95F02", "#7570B3", "#E7298A"
)
# Plot
ggplot(df_all, aes(x = Treatment_agg, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Plant, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = "Treatments",
    y = "Relative Abundance",
    title = "Top 40 Family-level Microbial Abundance by Plant Species"
  ) +
  scale_fill_manual(values = colorblind_palette) +
  scale_x_discrete(labels = c(
    "Clamp_Optimized_Sterilized" = "Clamp Optimized (Sterilized)",
    "Clamp_Standard_Sterilized" = "Clamp Standard (Sterilized)",
    "NoClamp_Sterilized" = "No Clamp (Sterilized)",
    "Clamp_Washed" = "Clamp (Washed)",
    "NoClamp_Washed" = "No Clamp (Washed)"
  ))

