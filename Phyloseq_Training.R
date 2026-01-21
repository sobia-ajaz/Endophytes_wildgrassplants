library(phyloseq)
library(qiime2R)
library(dplyr)
library(ggplot2)
library(data.table)
library(tibble)
library(readxl)
library(rstudioapi)
library(tidyr)
library(vegan)  # For diversity calculations


setwd(dirname(getActiveDocumentContext()$path))

ASVs<-read_qza("id-filtered-table.qza")
names(ASVs)
#To access the actual data stored within the object, access the data as below:
ASVs$data[1:5,1:2] #show first 2 samples and first 5 ASVs

metadata<-read_q2metadata("metadata-16S-Training.txt")
head(metadata) # show top lines of metadata

taxonomy<-read_qza("taxonomy_SILVA.qza")
head(taxonomy$data)
taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)

physeq<-qza_to_phyloseq(
  features="id-filtered-table.qza",
  tree="rooted-tree.qza",
  "taxonomy_SILVA.qza",
  metadata = "metadata-16S-Training.txt"
)
physeq

###import env data
save.image(file = "Training.RData")


# setting the seed to one value in order to created reproducible results
set.seed(1)  

S16_scaled <- rarefy_even_depth(physeq, sample.size=50000, replace=FALSE, rngseed = 1) 




