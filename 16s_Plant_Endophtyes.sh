#!/bin/bash

#activate qiime2
conda activate qiime2-amplicon-2024.2

#Import the files
qiime tools import --input-path Training_Dataset --type 'SampleData[PairedEndSequencesWithQuality]' --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end_16S-Plant-Training.qza

#To see the imported data
qiime demux summarize --i-data demux-paired-end_16S-Plant-Training.qza --o-visualization demux-paired-end_16S-Plant-Training.qzv

#RRIMERS (16S) V4
#Amplicon size: ~250bp
#Primer forward (515F): (19 bp)GTGCCAGCMGCCGCGGTAA
#Primer reverse (806R): (20 bp)GGACTACNVGGGTHTCTAAT

#In the following command, --p-front-f is the forward primer (5’-3’),
#--p-front-r the reverse primer (5’-3’)
#You can also included --p-adapter-f, which is the reverse complement of the reverse primer, and --p-adapter-r, which is the reverse complement of the forward primer
#These adapters would be present in the reads if the amplicon region is shorter than the read length - in this case both primers are typically present in each forward and reverse read
# and you should treat them like “linked adapters”, where the forward primer is linked with the reverse complement of the reverse primer, and vice versa.
#--p-discard-untrimmed discards any untrimmed reads. In this command, reads will only be discarded if they don’t have either the normal primer or the linked adapter (just having one is ok)
#--verbose gives you a detailed output of how many reads and base pairs were trimmed off in each sample (by default if this was run as a script on the ccv it would be stored in the slurm file)

#Cutadapt: 
#https://docs.qiime2.org/2024.2/plugins/available/cutadapt/trim-paired/
#https://cutadapt.readthedocs.io/en/stable/guide.html#five-prime-adapters

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end_16S-Plant-Training.qza --p-front-f GTGCCAGCMGCCGCGGTAA --p-front-r GGACTACNVGGGTHTCTAAT --p-discard-untrimmed --o-trimmed-sequences demux-cutadapt_16S-Plant-Training.qza --verbose

#To see the trimmed data
qiime demux summarize --i-data demux-cutadapt_16S-Plant-Training.qza --o-visualization demux-cutadapt_16S-Plant-Training.qzv

#DADA2

qiime dada2 denoise-paired --i-demultiplexed-seqs demux-cutadapt_16S-Plant-Training.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 220 --p-trunc-len-r 180 --o-table feature-table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoizing_stats.qza

#Visualization

qiime metadata tabulate --m-input-file denoizing_stats.qza --o-visualization denoizing_stats.qzv

qiime feature-table summarize --i-table feature-table.qza --o-visualization feature-table.qzv --m-sample-metadata-file metadata-16S.txt 

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

qiime metadata tabulate --m-input-file denoizing_stats.qza --o-visualization denoizing_stats.qzv

#Export to tsv for Visualization 

qiime tools export --input-path feature-table.qza --output-path exported-table

cd exported-table

biom convert -i feature-table.biom -o feature-table.tsv --to-tsv



#Taxonomic assignment
#Download the pre-trained reference taxonomy classifier from https://docs.qiime2.org/2024.2/data-resources/
#or you can train your own classifier
 
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy_SILVA.qza



#For visualization

qiime metadata tabulate --m-input-file taxonomy_SILVA.qza --o-visualization taxonomy_SILVA.qzv


#Taxabarplots

qiime taxa barplot --i-table feature-table.qza  --i-taxonomy taxonomy_SILVA.qza --m-metadata-file metadata-16S.txt --o-visualization taxa-bar-plots.qzv




#Creating a tree to do phylogenetic diversity analyzes. For more info follow the link
#https://docs.qiime2.org/2020.2/plugins/available/phylogeny/


qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

#Alpha and beta diversity analysis
qiime diversity alpha-rarefaction --i-table feature-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 95230 --m-metadata-file metadata-16S.txt --o-visualization alpha-rarefaction.qzv

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table feature-table.qza  --p-sampling-depth 95230 --m-metadata-file metadata-16S.txt --output-dir core-metrics-results


