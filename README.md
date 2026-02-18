# Endophytes_wildgrassplants
# Optimised Workflow for Detecting Root Endophytic Bacterial Communities
This repository contains scripts used in the study:
# “Unveiling Hidden Endophytes by Optimising Identification of Endophytic Bacterial Communities from Wild Grassland Plant Roots.”

This project uses a reproducible bioinformatics pipeline to process 16S rRNA amplicon sequencing data for identifying bacterial endophytic communities from plant roots. The workflow integrates QIIME2 processing, downstream analysis in R using phyloseq, and comparative analyses with and without host-derived sequences. <br>



**Step 1 — QIIME2 Processing** <br>
All primary sequence processing was conducted in QIIME2.<br>
This Bash script implements a Qiime2 pipeline for processing paired-end 16S rRNA amplicon data (V4 region). It begins by importing raw sequencing files, trimming forward and reverse primers with Cutadapt, and denoising the reads using DADA2 to generate a feature table and representative sequences. Taxonomic assignment is performed with a pre-trained SILVA classifier, and a phylogenetic tree is constructed via MAFFT alignment and FastTree. The pipeline concludes with alpha rarefaction and core diversity metrics, all consistently referencing the same metadata file (metadata-16S.txt) throughout the workflow.
**16s_Plant_Endophytes.sh** <br> script used for this step <br>
**Outputs:** <br>
Feature table (ASV abundance)<br>
Representative sequences<br>
Taxonomy table<br>
QIIME2 artefacts (.qza files)<br>
Phylogenetic Tree <br>

**Step 2 — Export and Import into R** <br>
QIIME2 outputs were exported and imported into R <br>
**Actions:** <br>
Convert feature table to biom format <br>
Import taxonomy and metadata <br>
Create phyloseq object <br>
**Output:** <br>
Integrated phyloseq dataset <br>


