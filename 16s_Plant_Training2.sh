#!/bin/bash 


#activate qiime2
conda activate qiime2-amplicon-2024.2


#Filter feature table to remove feature IDs belongs to singletons (prepare metadata file that contain feature IDS to retain)
qiime feature-table filter-features --i-table feature-table.qza --m-metadata-file filtered_feature-table.txt --o-filtered-table id-filtered-table.qza

#visualisation

qiime feature-table summarize --i-table id-filtered-table.qza --o-visualization id-filtered-table.qzv --m-sample-metadata-file metadata-16S-Training.txt


#Removing sequences from repseqs that belong to singletons
 
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table id-filtered-table.qza --o-filtered-data id-filtered-rep-seqs.qza

#visualisation

qiime feature-table tabulate-seqs --i-data id-filtered-rep-seqs.qza --o-visualization id-filtered-rep-seqs.qzv



#Taxonomic assignment
#Download the pre-trained reference taxonomy classifier from https://docs.qiime2.org/2024.2/data-resources/
#or you can train your own classifier
 
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-515-806-nb-classifier.qza --i-reads id-filtered-rep-seqs.qza --o-classification taxonomy_SILVA.qza



#For visualization

qiime metadata tabulate --m-input-file taxonomy_SILVA.qza --o-visualization taxonomy_SILVA.qzv


#Taxabarplots

qiime taxa barplot --i-table id-filtered-table.qza  --i-taxonomy taxonomy_SILVA.qza --m-metadata-file metadata-16S-Training.txt --o-visualization taxa-bar-plots.qzv




#Creating a tree to do phylogenetic diversity analyzes. For more info follow the link
https://docs.qiime2.org/2020.2/plugins/available/phylogeny/


qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#Alpha and beta diversity analysis
qiime diversity alpha-rarefaction --i-table feature-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 95230 --m-metadata-file metadata1.txt --o-visualization alpha-rarefaction.qzv

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table feature-table.qza  --p-sampling-depth 95230 --m-metadata-file metadata1.txt --output-dir core-metrics-results
