
#######################################################
 ## DATA LOAD + FORMATTING 
#######################################################

setwd("/AWI_MPI/FRAM/CTD/Rstats")
setwd("/Users/mwietz/ownCloud - mwietz@owncloud.mpi-bremen.de/AWI_MPI/FRAM/IceExp/Rstats")

library(gtools)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggridges)
library(vegan)
library(ampvis2)
library(FEAST)
library(UpSetR)
library(ANCOMBC)
library(phyloseq)
library(iNEXT)
library(olsrr)
library(coin)
#load("IceExp.Rdata")


#######################################################
  ### LOAD 16S ASVs + TAX ###
#######################################################

# Load ASV table
# Filter 
ASV.bac <- read.table(
  "../dada/IceExp_16S_asv.txt",
  h = T, sep = "\t",
  check.names = F) %>% 
  filter(rowSums(.>= 3) >= 2)

# Load taxonomy table
TAX.bac <- read.table(
  "../dada/IceExp_16S_tax.txt",
  h = T, sep = "\t")

# Match TAX after abundance filter
TAX.bac <- TAX.bac[row.names(ASV.bac),]

# Subtract negative counts from each column
ASV.bac = ASV.bac - ASV.bac$file032_16S_clip

# Set neg. values to zero
ASV.bac[ASV.bac < 0] <- 0

# Remove Mitochondria and Chloroplasts
# Remove more potential contaminants
TAX.bac <- TAX.bac %>% 
  filter(
  !grepl('Chloroplast|Propionibacteriales|Lactobacillales', Order) &
  !grepl('Corynebacteriaceae|Bacillaceae|Streptococcaceae|Weeksellaceae|Enterococcaceae|Staphylococcaceae|Erysipelotrichaceae|Carnobacteriaceae|Mitochondria', Family))

# Match TAX after contaminant removal
ASV.bac <- round(ASV.bac[row.names(TAX.bac), ], 0)

#################################################

# Rename NAs with last known taxrank + "uc"
k <- ncol(TAX.bac)-1
for (i in 2:k) {
if (sum(is.na(TAX.bac[, i])) >1) {
  temp <- TAX.bac[is.na(TAX.bac[, i]), ]
  for (j in 1:nrow(temp)) {
    if (sum(is.na(
      temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
      temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
      temp[j, (i+1):(k+1)] <- temp[j, i]
  }
}
  TAX.bac[is.na(TAX.bac[, i]), ] <- temp}
  if (sum(is.na(TAX.bac[, i]))==1) {
    temp <- TAX.bac[is.na(TAX.bac[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
    temp[i] <- paste(temp[(i-1)], " uc", sep="")
    temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.bac[is.na(TAX.bac[, i]),] <- temp
  }
}
TAX.bac[is.na(TAX.bac[, (k+1)]), (k+1)] <- paste(
  TAX.bac[is.na(TAX.bac[, (k+1)]), k], " uc", sep="")

# shorten/modify names
TAX.bac <- TAX.bac %>%
  mutate(across(everything(),~gsub("Clade ","SAR11 Clade ", .))) %>%
  mutate(across(everything(),~gsub("_clade","", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_NAC11-7_lineage","NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("Marinimicrobia (SAR406 clade)","SAR406", ., fixed=T)))


#######################################################
### LOAD 18S ASVs + TAX ###
#######################################################

# Load ASV table
# Filter 
ASV.euk <- read.table(
  "../dada/IceExp_18S_asv.txt",
  h = T, sep = "\t",
  check.names = F) %>% 
  filter(rowSums(.>= 3) >= 2)

# Load taxonomy table
# Rename PR2 taxranks 
# Remove Supergroup to align with 16S taxonomy
TAX.euk <- read.table(
  "../dada/IceExp_18S_tax.txt",
  h = T, sep = "\t") %>% 
 setNames(c(
  "Kingdom","Supergroup","Phylum","Class",
  "Order","Family","Genus","Species")) %>% 
  dplyr::select(-Supergroup)

# Match TAX after abundance filter
TAX.euk <- TAX.euk[row.names(ASV.euk),]

# Subtract negative counts from each column
ASV.euk = ASV.euk - ASV.euk$file032_18S_clip

# Set neg. values to zero
ASV.euk[ASV.euk < 0] <- 0

# Remove Metazoa and potential contaminants
TAX.euk <- TAX.euk %>% 
  filter(!grepl('Basidiomycota|Embryophyceae', Class) & Phylum != "Metazoa")

# Match TAX after contaminant removal
ASV.euk <- round(ASV.euk[row.names(TAX.euk), ], 0)

#################################################

# Rename NAs with last known taxrank + "uc"
k <- ncol(TAX.euk)-1
for (i in 2:k) {
  if (sum(is.na(TAX.euk[, i])) >1) {
    temp <- TAX.euk[is.na(TAX.euk[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], " uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX.euk[is.na(TAX.euk[, i]), ] <- temp}
  if (sum(is.na(TAX.euk[, i]))==1) {
    temp <- TAX.euk[is.na(TAX.euk[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], " uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX.euk[is.na(TAX.euk[, i]),] <- temp
  }
}
TAX.euk[is.na(TAX.euk[, (k+1)]), (k+1)] <- paste(
  TAX.euk[is.na(TAX.euk[, (k+1)]), k], " uc", sep="")

## shorten/modify names 
TAX.euk <- TAX.euk %>%
  mutate(across(everything(),~gsub("Dino-Group-I-Clade-","Dino-I-", .))) %>%
  mutate(across(everything(),~gsub("Dino-Group-II-Clade-","Dino-II-", .))) %>%
  mutate(across(everything(),~gsub("Polar-centric-","", .))) %>%
  mutate(across(everything(),~gsub("Radial-centric-basal-","", .))) %>%
  mutate(across(everything(),~gsub("Chrysophyceae_Clade-","Chrysophyceae ", .))) %>%
  mutate(across(everything(),~gsub("Stephanoecidae_Group_","Stephanoecidae ", .))) %>%
  mutate(across(everything(),~gsub("Pirsonia_Clade_","Pirsonia", .))) %>%
  mutate(across(everything(),~gsub("_X|_XX|_XXX|_XXXX"," uc", .))) 


#######################################################
   ###  SAMPLE INFO / METADATA
#######################################################

ENV <- read.table(
  "metadata.txt", h=T, sep="\t", check.names=F,
  stringsAsFactors=F, skipNul=T) 

# Factorize parameters of interest
ENV$condition <- factor(ENV$condition, levels=c(
  "OrgSW","OrgIce","IceMelt","Ctr"))

# Separate BAC/EUK 
ENV.bac <- filter(ENV, locus_tag=="16S")
ENV.euk <- filter(ENV, locus_tag=="18S")


#######################################################
  ###  PROCESSING / FORMATTING
#######################################################

# Calculate relative abundances
ASV.bac.rel = as.data.frame(
  apply(ASV.bac, 2, function(x) x / sum(x) * 100)) 
ASV.euk.rel = as.data.frame(
  apply(ASV.euk, 2, function(x) x / sum(x) * 100)) 

# Prepare Ampvis objects
amp.bac <- data.frame(
  OTU = rownames(ASV.bac),
  ASV.bac, TAX.bac, check.names=F)
amp.bac <- amp_load(amp.bac, ENV.bac)

amp.euk <- data.frame(
  OTU = rownames(ASV.euk),
  ASV.euk, TAX.euk, check.names=F)
amp.euk <- amp_load(amp.euk, ENV.euk)

# Find most abundant ASV
topBac <- ASV.bac.rel %>%
  rownames_to_column("asv") %>%
  right_join(TAX.bac %>% rownames_to_column("asv"), by="asv") %>%
  dplyr::select(asv, Genus, where(is.numeric)) %>%
  reshape2::melt() %>%
  drop_na() %>%
  group_by(Genus, asv) %>%
  summarise(mean=mean(value)) %>%
  slice_max(mean, n = 1) %>%
  ungroup()

topEuk <- ASV.euk.rel %>%
  rownames_to_column("asv") %>%
  right_join(TAX.euk %>% rownames_to_column("asv"), by="asv") %>%
  dplyr::select(asv, Genus, where(is.numeric)) %>%
  reshape2::melt() %>%
  drop_na() %>%
  group_by(Genus, asv) %>%
  summarise(mean=mean(value)) %>%
  slice_max(mean, n = 1) %>%
  ungroup()

# Define function for presence-absence analysis
# Function for presence-absence matrix
presAbs <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)}


###############################

# Remove temporary data
rm(temp,i,j,k)

# Save data
save.image("IceExp.Rdata")
