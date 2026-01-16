library(phyloseq)
library(microbiomeMarker)
library(ggplot2)

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



########################. LEfSe Multi-Level Species Analysis
lefse_results <- run_lefse(family_physeq, group = "Plant")


