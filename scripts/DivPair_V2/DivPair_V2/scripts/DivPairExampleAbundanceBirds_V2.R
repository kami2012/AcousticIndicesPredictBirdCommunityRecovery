# DivPair Example, version 2
# 2024-11-04

# DivPair Demo
# by Joerg and Anne
# updated by Oliver, 15.04.2024
# update, Oliver, 15.9.2024
rm(list=ls())

# Install this version of tidytree to get rid of lots of phylogenetic tree warnings
# library(devtools)
# install_version("tidytree", version = "0.4.2", repos = "http://cran.us.r-project.org")
library(tidytree)

library(iNEXT.beta3D)
library(snowfall)
library(dplyr)

# My working directory
#setwd("C:/Users/mitesser/mxinfected/Wald/DivPair/05_ToDO-PDFD")
#setwd("C:/Users/mitesser/mxinfected/Wald/DivPair/06_DivPairMareike/iNEXTbeta3D_incidence/iNEXTbeta3D_incidence")

# Load script by Anne Chao (updated by Oliver)
source("scripts/DivPairCoverage_V2.R")

## ================================================ Example birds ================================================ ##
## Use abundance data "Birds" that Joreg sent as an example.
## There are 177 species (rows) and 5 assemblages (columns) in "Birds" data.
## The code below is to compute observed and coverage-based similarity between any two assemblages.

# Load "Birds" data
Birds <- read.csv2("data/Birds.csv", header = T, row.names = 1)
# Flip bird data
tab<-t(Birds)

# Load tree
Birds_tree <- ape::read.tree("data/Birds_phylo.tre")

# Load traits and rearrange data
Birds_traits <- read.csv2("data/Birds_traits.csv")
rownames(Birds_traits) <- Birds_traits[,1]
Birds_traits$Habitat<-as.factor(Birds_traits$Habitat)
Birds_traits$TrophicNiche<-as.factor(Birds_traits$TrophicNiche)
Birds_traits<-Birds_traits[,-1]
# Calcualte trait distance matrix
library(cluster)
Birds_distances <- as.matrix(daisy(Birds_traits, metric = "gower"))

# Compute Cmax in the joint assemblages for alpha reference samples by the function "DataInfobeta3D" in the package "iNEXT.beta3D"
# For Birds data, Cmax_joint = 0.9103252
Cmax_joint = sapply(1:ncol(Birds), function(i) 
  sapply(1:ncol(Birds), function(j) DataInfobeta3D( Birds[,c(i,j)] )$`SC(2n)`[4])) %>% min  

## Use function iNEXTbeta3D_pair to compute similarity index.
tt_TD  = iNEXTbeta3D_pair3D(tab, div0="TD", SC = Cmax_joint, parallel=T)
tt_PD  = iNEXTbeta3D_pair3D(tab, div0="PD", SC = Cmax_joint, PDTree0 = Birds_tree, parallel=T)
tt_FD  = iNEXTbeta3D_pair3D(tab, div0="FD", SC = Cmax_joint, FDdistM0 = Birds_distances, parallel=T)

# Just to check - without parallelisation
tt_TD2  = iNEXTbeta3D_pair3D(tab, div0="TD", SC = Cmax_joint)
tt_PD2  = iNEXTbeta3D_pair3D(tab, div0="PD", SC = Cmax_joint, PDTree0 = Birds_tree)
tt_FD2  = iNEXTbeta3D_pair3D(tab, div0="FD", SC = Cmax_joint, FDdistM0 = Birds_distances)

# Save data
save(tt_TD, tt_PD, tt_FD , file= "results/tt_org_3.rda")
load("results/tt_org_3.rda")

# at the following step select which species distance goes into the analysis
# TD as an example
colnames(tt_TD) # jac, sor - > q = 0, hor -> q = 1, mor_hor -> q = 2

# Results
dis_birds_q0obs <- matrix(tt_TD$ResultTable[,"sor_obs"], ncol = nrow(tab), byrow = T)
# Convert to dissimilarity convertieren
dis_birds_q0obs <- (1 - dis_birds_q0obs)
row.names(dis_birds_q0obs)<-row.names(tab)

dis_birds_q0obs
