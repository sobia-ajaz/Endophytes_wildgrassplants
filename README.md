# Endophytes_wildgrassplants
# Optimised Workflow for Detecting Root Endophytic Bacterial Communities
This repository contains scripts used in the study:
# “Unveiling Hidden Endophytes by Optimising Identification of Endophytic Bacterial Communities from Wild Grassland Plant Roots.”

This project uses a reproducible bioinformatics pipeline to process 16S rRNA amplicon sequencing data for identifying bacterial endophytic communities from plant roots. The workflow integrates QIIME2 processing, downstream analysis in R using phyloseq, and comparative analyses with and without host-derived sequences. <br>



Step 1 — QIIME2 Processing <br>
All primary sequence processing was conducted in QIIME2.
Key Steps:
Import paired-end FASTQ files
Quality filtering and trimming
Denoising using DADA2
ASV inference
Chimera removal
Generation of feature table
Taxonomic assignment using reference database
Outputs:
Feature table (ASV abundance)
Representative sequences
Taxonomy table
QIIME2 artefacts (.qza files)
Phylogenetic Tree


