#creating phyloseq object using qiime2 output files directly

#you need 4 files:

#1) ASV table (feature-table.qza)
#2) Metadata (metadata1.txt)
#3) Taxonomic composition (taxonomy_unite_ITS.qza)
#4) Phylogenetic tree (rooted-tree.qza)

#You can import these files into R directly. 

#Note: be aware that you need to install qiime2R using following link:
  

#if
#(!requireNamespace("devtools",
#                  quietly =
#                 TRUE)){install.packages("devtools")}
#devtools::install_github("jbisanz/qiime2R")

#Once, it’s installed run following commands in R to create phyloseq object
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
######################
#install.packages(
#  "microViz",
#  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
#)


######
####



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

setwd(dirname(getActiveDocumentContext()$path))

ASVs<-read_qza("feature-table-16S-Plant-all.qza")
names(ASVs)
#To access the actual data stored within the object, access thenano  data as below:
ASVs$data[1:5,1:2] #show first 2 samples and first 5 ASVs

metadata<-read_q2metadata("metadata-16S-Plant.txt")
head(metadata) # show top lines of metadata

taxonomy<-read_qza("taxonomy_SILVA-16S-Plant-all.qza")
head(taxonomy$data)
taxonomy<-parse_taxonomy(taxonomy$data)
head(taxonomy)

physeq<-qza_to_phyloseq(
  features="feature-table-16S-Plant-all.qza",
  tree="rooted-tree-16S-Plant-all.qza",
  "taxonomy_SILVA-16S-Plant-all.qza",
  metadata = "metadata-16S-Plant.txt"
)
physeq
####
###
library(readxl)
Env_data <- read_excel("Metadata.xlsx", sheet = "Sheet1")
                                
###import env data
save.image(file = "Phyloseq-Plant.RData")



# setting the seed to one value in order to created reproducible results
set.seed(1)  

# scaling the data to the smallest samples. Note: rngseed is similar to set.seed
sample_sums(physeq)
sort(sample_sums(physeq))
S16_scaled <- rarefy_even_depth(physeq, sample.size=47289, replace=FALSE, rngseed = 1) 
# rarefying reads to even depth (Note that this is not recommended by the phyloseq authors)


Meta<-sample_data(S16_scaled)


Rich_16S<-unlist(estimate_richness(S16_scaled, measures="Observed"))
Shannon_16S<-unlist(estimate_richness(S16_scaled, measures="Shannon"))

df1<-data.frame(Meta, Rich_16S,Shannon_16S,SampleID=rownames(Meta))

df1_joint<- merge(Env_data, df1, by='SampleID',all.x=TRUE)
names(df1_joint)

Data_set<-data.frame(
  ID=df1_joint$SampleID,
  Lat=df1_joint$Lat,
  Long=df1_joint$Long,
  Rich=df1_joint$Rich_16S,
  Shannon=df1_joint$Shannon_16S,
  Treat=as.factor(df1_joint$Treatment.x),
  Site=as.factor(df1_joint$Site.x),
  Moist=df1_joint$`Mean Soil Moisture (Gravimetric %)`,
  Yield=df1_joint$`Grass Yield`,
  Prec=df1_joint$`Mean 10yr Precip (mm)`,
  BD=df1_joint$`Mean Bulk Density`,
  DNA=df1_joint$`Soil DNA TOT (ng/µL in 100µL elution volume)`,
  PrecRange=df1_joint$`Range 10yr Precip`,
  FloodRisk=df1_joint$`Flood Risk`
  )



ggplot(Data_set, aes(x=Treat, y=Rich)) +
  geom_boxplot(aes(color=Treat),outlier.shape = NA, size=0.5, show.legend = FALSE) +
  geom_point(aes(color=Treat, shape = Site),
             position=position_jitterdodge(dodge.width=0.4), size=3,show.legend = FALSE) +
  scale_color_manual(values=c("gold3","lightskyblue3"), name="") +
  scale_shape_manual(name="",values=seq(1:20))+
  scale_x_discrete(labels = c("Control","Drought"))+
  labs(y= "Richness (16S)")+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


################################
##########plot with Risk of Flood vs shelter
######################
#####create combined factor

Shelter<-as.factor(Data_set$Treat)
levels(Shelter) <- list(no="Control", yes = "Drought")

RiskF<-as.factor(Data_set$FloodRisk)

Type<-interaction(Shelter,RiskF)


df<-data.frame(cbind(Type,Data_set))


Data_set<-cbind(Data_set,Type)
rownames(Data_set)<-Data_set$ID
sample_data(S16_scaled) <- sample_data(Data_set)

######
##

##################
########
ggplot(df, aes(x=Type, y=Moist)) +
  geom_boxplot(aes(color=Type),outlier.shape = NA, size=0.5, show.legend = FALSE) +
  geom_point(aes(color=Type, shape = Type),
             position=position_jitterdodge(dodge.width=0.9), size=3,show.legend = FALSE) +
  scale_color_manual(values=c("gold1", "gold4","lightskyblue1","lightskyblue3"), name="") +
  scale_shape_manual(name="",values=c(19, 18,19,18))+
  scale_x_discrete(labels = c("No Shelter 
  High Flood Risk","Shelter 
  High Flood Risk","No Shelter 
  Low Flood Risk", "Shelter 
  Low Flood Risk" ))+
  labs(y= "Soil Moisture (%)")+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


###########################
######
ggplot(df, aes(x=Type, y=Rich)) +
  geom_boxplot(aes(color=Type),outlier.shape = NA, size=0.5, show.legend = FALSE) +
  geom_point(aes(color=Type, shape = Type),
             position=position_jitterdodge(dodge.width=0.9), size=3,show.legend = FALSE) +
  scale_color_manual(values=c("gold1", "gold4","lightskyblue1","lightskyblue3"), name="") +
  scale_shape_manual(name="",values=c(19, 18,19,18))+
  scale_x_discrete(labels = c("No Shelter 
  High Flood Risk","Shelter 
  High Flood Risk","No Shelter 
  Low Flood Risk", "Shelter 
  Low Flood Risk" ))+
  labs(y= "16S ASV Richness")+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


##########################
#######shannon
ggplot(df, aes(x=Type, y=Shannon)) +
  geom_boxplot(aes(color=Type),outlier.shape = NA, size=0.5, show.legend = FALSE) +
  geom_point(aes(color=Type, shape = Type),
             position=position_jitterdodge(dodge.width=0.9), size=3,show.legend = FALSE) +
  scale_color_manual(values=c("gold1", "gold4","lightskyblue1","lightskyblue3"), name="") +
  scale_shape_manual(name="",values=c(19, 18,19,18))+
  scale_x_discrete(labels = c("No Shelter 
  High Flood Risk","Shelter 
  High Flood Risk","No Shelter 
  Low Flood Risk", "Shelter 
  Low Flood Risk" ))+
  labs(y= "16S ASV Shannon index")+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

###########################
##########

plot(df$Shannon,df$Moist)
plot(df$Moist,log(df$Rich),ylab="Log of 16S Richness",xlab="Moisture",main="r = 0.48, p-value < 0.01")
cor.test(df$Moist,log(df$Rich))
cor.test(log(df$Moist),log(df$Rich))


plot(df$Moist,log(df$Rich))
#cor.test(df$Moist,df$Shannon)


No_Shelter_HF_Risk<-subset(df,Type=="no.High")
points(No_Shelter_HF_Risk$Moist,log(No_Shelter_HF_Risk$Rich),cex=1.5,pch=19,col="lightskyblue3")

No_Shelter_LF_Risk<-subset(df,Type=="no.Low")
points(No_Shelter_LF_Risk$Moist,log(No_Shelter_LF_Risk$Rich),cex=1.5,pch=19,col="lightskyblue1")


Yes_Shelter_HF_Risk<-subset(df,Type=="yes.High")
points(Yes_Shelter_HF_Risk$Moist,log(Yes_Shelter_HF_Risk$Rich),cex=1.5,pch=19,col="gold1")

Yes_Shelter_LF_Risk<-subset(df,Type=="yes.Low")
points(Yes_Shelter_LF_Risk$Moist,log(Yes_Shelter_LF_Risk$Rich),cex=1.5,pch=19,col="gold4")


legend("topleft", legend=c("Y_Shelter_LRF","Y_Shelter_HRF","No_Shelter_LRF","No_Shelter_HRF"),
       col=c("gold4","gold1","lightskyblue1","lightskyblue3"), pch=c(19,19,19,19), cex=0.8)

#################
##########
plot(df$Moist,df$Shannon,ylab="16S Shannon index",xlab="Moisture",main="r = 0.57, p-value < 0.0001")
cor.test(df$Moist,df$Shannon)


No_Shelter_HF_Risk<-subset(df,Type=="no.High")
points(No_Shelter_HF_Risk$Moist,No_Shelter_HF_Risk$Shannon,cex=1.5,pch=19,col="lightskyblue3")

No_Shelter_LF_Risk<-subset(df,Type=="no.Low")
points(No_Shelter_LF_Risk$Moist,No_Shelter_LF_Risk$Shannon,cex=1.5,pch=19,col="lightskyblue1")


Yes_Shelter_HF_Risk<-subset(df,Type=="yes.High")
points(Yes_Shelter_HF_Risk$Moist,Yes_Shelter_HF_Risk$Shannon,cex=1.5,pch=19,col="gold1")

Yes_Shelter_LF_Risk<-subset(df,Type=="yes.Low")
points(Yes_Shelter_LF_Risk$Moist,Yes_Shelter_LF_Risk$Shannon,cex=1.5,pch=19,col="gold4")



legend("topleft", legend=c("Y_Shelter_LRF","Y_Shelter_HRF","No_Shelter_LRF","No_Shelter_HRF"),
       col=c("gold4","gold1","lightskyblue1","lightskyblue3"), pch=c(19,19,19,19), cex=0.8)

###############
#####

##############
#####
Data_set<-cbind(Data_set,Type)
rownames(Data_set)<-Data_set$ID
sample_data(S16_scaled) <- sample_data(Data_set)

########
##
Phylum1 <- tax_glom(S16_scaled, taxrank = "Phylum")


#####################
#########Plot Phyla by
#sample_names(S16_scaled) %>% head()
#taxa_names(S16_scaled) %>% head()
#S16_scaled %>%tax_fix(unknowns = NULL) %>%
 # tax_table() %>%
  #head(50)

#S16_scaled%>% tax_fix(unknowns = c("uncultured"))
#S16_scaled %>%
 # tax_fix() %>%
  #tax_agg(rank = "Genus")



p2 <- comp_barplot(Phylum1,
                   tax_level = "Phylum", n_taxa = 20, group_by = "Type",
                   sample_order = "euclidean", bar_outline_colour = NA
) %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect") &
  coord_flip() & labs(x = NULL, y = NULL) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2




p2 <- comp_barplot(Phylum1,
                   tax_level = "Phylum", n_taxa = 20, group_by = "Treat",
                   sample_order = "euclidean", bar_outline_colour = NA
) %>%
  patchwork::wrap_plots(nrow = 3, guides = "collect") &
  coord_flip() & labs(x = NULL, y = NULL) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2


otu_table(Phylum1)
taxa<-as.data.frame(tax_table(Phylum1))
PhyNames<-taxa$Phylum
otu<-as.data.frame(otu_table(Phylum1))

B_table<-as.matrix(t(otu_table(Phylum1)))@ .Data #transform in order to work in vegan (extracts tables). on scaled data
colnames(B_table)<-PhyNames
Env<-as.data.frame(sample_data(Phylum1))


#############################
##############Hellinger transformation to use linear models
##############which is good for variance partitioning
##########
library(vegan)
B_table_hel<-decostand(B_table, method="hellinger") 
####Legendre, P. & Gallagher, E.D. (2001) Ecologically meaningful transformations for ordination of species data. Oecologia 129, 271–280.

dim(B_table_hel)


#Fungi_full_model<-rda(Fungi_table_hel~YLD+Prot+AMF+Condition(Y+PCNM1+PCNM5+PCNM8),data=e_data) #rda response is AMf table , constraining variables are p and pH conditional on spatial model meaning that it first plots the spatial varienility then takes out the residuals and plots them ahgainst the p and pH (ie taknginto account all that space already explains)
#Fungi_full_model<-rda(Fungi_table_hel~YLD+Prot+AMF,data=e_data) #rda response is AMf table , constraining variables are p and pH conditional on spatial model meaning that it first plots the spatial varienility then takes out the residuals and plots them ahgainst the p and pH (ie taknginto account all that space already explains)
B_PCoA<-rda(B_table_hel~Moist+Yield+Prec+BD+DNA,data=Data_set) 
summary(B_PCoA)

####Plot constrained model for pure effect of abiotic gradients_Supplementary Materials
plot(B_PCoA,type="n",cex=0.001,xlab="RDA 1 (30%)",ylab="RDA2 (1%)")
text(B_PCoA, dis="cn", scaling="sites")
text(B_PCoA, dis="sp", cex=0.4,scaling="sites")

PCoA_ordSc<-scores(B_PCoA,choices=c(1,2))
PCoA_ordSc$sites
PCoA_ordScAxes<-as.data.frame(PCoA_ordSc$sites)



no_Low<-subset(PCoA_ordScAxes,Type=="no.Low")
points(no_Low$RDA1,no_Low$RDA2,cex=1.8,pch=18,col="gold1")


yes_Low<-subset(PCoA_ordScAxes,Type=="yes.Low")
points(yes_Low$RDA1,yes_Low$RDA2,cex=1.8,pch=19,col="gold4")


no_High<-subset(PCoA_ordScAxes,Type=="no.High")
points(no_High$RDA1,no_High$RDA2,cex=1.8,pch=18,col="lightblue1")


yes_High<-subset(PCoA_ordScAxes,Type=="yes.High")
points(yes_High$RDA1,yes_High$RDA2,cex=1.8,pch=19,col="lightblue4")


legend("topleft", legend=c("Y_Shelter_LRF","Y_Shelter_HRF","No_Shelter_LRF","No_Shelter_HRF"),
       col=c("gold4","lightskyblue4","gold1","lightskyblue1"), pch=c(19,19,18,18), cex=0.8)

