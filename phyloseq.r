# Loading required packages
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(microbiome)

# Setting the working directory
setwd("D:/Tharuka Research/Methodology")

# Loading the dataset
# metadata table
metadata<-read_csv("datafiles/metadata_table.txt")
metadata <- data.frame(metadata)
rownames(metadata) <- metadata[,1]
metadata = metadata[,-1]
samples <- sample_data(metadata)

# OTU table
SVs<-read_qza("datafiles/table.qza")
otutable <- as.matrix(SVs$data)
OTU <- otu_table(otutable, taxa_are_rows = TRUE)

# taxonomy table 
taxonomy <-read_qza("datafiles/taxonomy.qza")
taxtable<-taxonomy$data[,1:2] 
# convert the table into a tabular split version
taxtable<- taxtable %>% separate(Taxon, sep="; ", c("Kingdom","Phylum","Class","Order","Family","Genus","Species")) 
# remove the first column
rownames(taxtable) <- taxtable[,1]
taxtable = taxtable[,-1]
TAX = tax_table(as.matrix(taxtable))

# Creating the phyloseq object of complete data
complete_physeq <- phyloseq(OTU,TAX,samples)

#aggregating taxa in the genus level
complete_agg <- tax_glom(complete_physeq, 'Genus')
complete_agg <- aggregate_taxa(complete_agg, 'Genus')

# keep only taxa that were observed at least twice
minTotRelAbun <- 5e-5
sums <- taxa_sums(complete_agg)
keepTaxa <- taxa_names(complete_agg)[which((sums / sum(sums)) > minTotRelAbun)]
filt_complete_agg <- prune_taxa(keepTaxa, complete_agg)

# select data of Different Agricultural zones.
Dry_physeq <- subset_samples(filt_complete_agg, AgriculturalZone =="Dry")
Int_physeq <- subset_samples(filt_complete_agg, AgriculturalZone =="Intermediate")
Wet_physeq <- subset_samples(filt_complete_agg, AgriculturalZone =="Wet")

# Saving the phyloseq objects
saveRDS(complete_physeq, "Complete data.rds")
saveRDS(complete_agg, "Aggregated data.rds")
saveRDS(filt_complete_agg, "Filtered Aggregated data.rds")
saveRDS(Dry_physeq, "Dry Zone.rds")
saveRDS(Int_physeq, "Intermediate Zone.rds")
saveRDS(Wet_physeq, "Wet Zone.rds")
