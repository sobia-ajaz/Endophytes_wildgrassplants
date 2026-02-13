# Endophytes_wildgrassplants
# Optimised Workflow for Detecting Root Endophytic Bacterial Communities
This repository contains scripts used in the study:
# “Unveiling Hidden Endophytes by Optimising Identification of Endophytic Bacterial Communities from Wild Grassland Plant Roots.”

This project uses a reproducible bioinformatics pipeline to process 16S rRNA amplicon sequencing data for identifying bacterial endophytic communities from plant roots. The workflow integrates QIIME2 processing, downstream analysis in R using phyloseq, and comparative analyses with and without host-derived sequences. <br>



**Step 1 — QIIME2 Processing** <br>
All primary sequence processing was conducted in QIIME2.<br>
**Key Steps:** <br>
Import paired-end FASTQ files<br>
Quality filtering and trimming<br>
Denoising using DADA2<br>
ASV inference<br>
Chimera removal<br>
Generation of feature table<br>
Taxonomic assignment using reference database<br>
**Outputs:**<br>
Feature table (ASV abundance)<br>
Representative sequences<br>
Taxonomy table<br>
QIIME2 artefacts (.qza files)<br>
Phylogenetic Tree <br>

**Step 2 — Export and Import into R ** <br>
QIIME2 outputs were exported and imported into R <br>
**Actions:** <br>
Convert feature table to biom format <br>
Import taxonomy and metadata <br>
Create phyloseq object <br>
**Output:** <br>
Integrated phyloseq dataset <br>


