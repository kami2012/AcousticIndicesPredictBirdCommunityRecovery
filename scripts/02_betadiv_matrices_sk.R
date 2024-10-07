rm(list=ls())

# Install this version of tidytree to get rid of lots of phylogenetic tree warnings
# library(devtools)
# install_version("tidytree", version = "0.4.2", repos = "http://cran.us.r-project.org")
library(tidytree)
library(iNEXT.beta3D)
library(snowfall)
library(dplyr)

# DivPair Demo by Joerg and Anne
# updated by Oliver, 15.04.2024
source("scripts/div_pair_coverage_OM.R")

## ================================================================================================ ##
## Use abundance data / community composition "BirdCommunityNamePhylogeny_Indices.csv" from Joerg.
## There are 85 assemblages (rows) and 334 birds (columns), identified by experts.
## The code below is to compute observed and coverage-based similarity and dissimilarity between any two assemblages.



# Load taxonomic diversity data
com_td <- read.csv("taxonomic_diversity/com_td.csv", row.names = 1)

# Load phylogenetic diversity data
com_pd <- read.csv("phylogenetic_diversity/com_pd.csv", row.names = 1)
load(file = "phylogenetic_diversity/bird_tree.rda")

# Load functional diversity data
com_fd <- read.csv("functional_diversity/com_fd.csv", row.names = 1)
load(file = "functional_diversity/trait_matrix.rda")

########################## Observed (with default SC) ################################

## Use function iNEXTbeta3D_pair to compute similarity index (with default (min) sample coverage (SC))
tt_TD  = iNEXTbeta3D_pair2(com_td, div0="TD", SC = NULL, parallel = T, cpus = 4)
save(tt_TD, file= "taxonomic_diversity/tt_tax.rda")

tt_PD  = iNEXTbeta3D_pair2(com_pd, div0="PD", SC = NULL, PDTree0 = tr)
save(tt_PD, file= "phylogenetic_diversity/tt_phy.rda")

tt_FD  = iNEXTbeta3D_pair2(com_fd, div0="FD", SC = NULL, FDdistM0 = distM)
save(tt_FD, file= "functional_diversity/tt_func.rda")


save(tt_TD, tt_PD, tt_FD , file= "data/tt_org.rda")
#load("data/tt_org.rda")

########################## Coverage-based (with user-specific SC) ################################
# !!!!!!!!!!!!!! NOT USED !!!!!!!!!!!!!!!!!


# Use function "DataInfobeta3D" in the package "iNEXT.beta3D"

# Convert columns to rows (birds) and rows to columns (plots)
com_td_transposed <- t(com_td)

# Compute sample size, observed species richness and sample coverage estimates (SC(n), SC(2n))
inextinfo <- DataInfobeta3D(com_td_transposed, diversity = "TD", datatype = "abundance") 

hist(inextinfo[1:86,"SC(n)"])
abline(v=0.75, col="red")

hist(inextinfo[1:86,"SC(2n)"])
abline(v=0.75, col="red")

## Use function iNEXTbeta3D_pair AND sample coverage to compute similarity index.

tt_TD93  = iNEXTbeta3D_pair2(birds, div0="TD", SC = 0.93)
save(tt_TD93, file= "data/tt_tax93.rda")

tt_PD93  = iNEXTbeta3D_pair2(birds, div0="PD", SC = 0.93, PDTree0 = tr)
save(tt_PD93, file= "data/tt_phy93.rda")

tt_FD93  = iNEXTbeta3D_pair2(birds_avo, div0="FD", SC = 0.93, FDdistM0 = distM)
save(tt_FD93, file= "data/tt_func93.rda")


save(tt_TD93, tt_PD93, tt_FD93 , file= "data/tt_org93.rda")
#load("data/tt_org93.rda")


colnames(tt_TD) # jac, sor - > q = 0, hor -> q = 1, mor_hor -> q = 2

# !!!!!!!!!!!!!! NOT USED !!!!!!!!!!!!!!!!!



## Compute dissimilarity matrices for the different diversities (taxonomic, pyhlogenetic and functional) and Hill numbers (q0, q1, q2)

#  Dissimilarity matrices TD 
# q0
sim_com_exp_q0est <- matrix(tt_TD[,"sor_est"], ncol = nrow(com_td), byrow = T)
dis_com_exp_q0est <- (1 - sim_com_exp_q0est) # convert to dissimilarity ('0' = completely similar, '1' = completely dissimilar)
row.names(dis_com_exp_q0est)<-row.names(com_td)

# q1
sim_com_exp_q1est <- matrix(tt_TD[,"hor_est"], ncol = nrow(com_td), byrow = T)
dis_com_exp_q1est <- (1 - sim_com_exp_q1est)
row.names(dis_com_exp_q1est)<-row.names(com_td)

# q2
sim_com_exp_q2est <- matrix(tt_TD[,"mor_hor_est"], ncol = nrow(com_td), byrow = T)
dis_com_exp_q2est <- (1 - sim_com_exp_q2est)
row.names(dis_com_exp_q2est)<-row.names(com_td)



#  Dissimilarity matrices PD
# q0
sim_com_exp_q0est_pd <- matrix(tt_PD[,"sor_est"], ncol = nrow(com_pd), byrow = T)
dis_com_exp_q0est_pd <- (1 - sim_com_exp_q0est_pd)
row.names(dis_com_exp_q0est_pd)<-row.names(com_pd)

# q1
sim_com_exp_q1est_pd <- matrix(tt_PD[,"hor_est"], ncol = nrow(com_pd), byrow = T)
dis_com_exp_q1est_pd <- (1 - sim_com_exp_q1est_pd)
row.names(dis_com_exp_q1est_pd)<-row.names(com_pd)

# q2
sim_com_exp_q2est_pd <- matrix(tt_PD[,"mor_hor_est"], ncol = nrow(com_pd), byrow = T)
dis_com_exp_q2est_pd <- (1 - sim_com_exp_q2est_pd)
row.names(dis_com_exp_q2est_pd)<-row.names(com_pd)



#  Dissimilarity matrices FD
# q0
sim_com_exp_q0est_fd <- matrix(tt_FD[,"sor_est"], ncol = nrow(com_fd), byrow = T)
dis_com_exp_q0est_fd <- (1 - sim_com_exp_q0est_fd)
row.names(dis_com_exp_q0est_fd)<-row.names(com_fd)

# q1
sim_com_exp_q1est_fd <- matrix(tt_FD[,"hor_est"], ncol = nrow(com_fd), byrow = T)
dis_com_exp_q1est_fd <- (1 - sim_com_exp_q1est_fd)
row.names(dis_com_exp_q1est_fd)<-row.names(com_fd)

# q2
sim_com_exp_q2est_fd <- matrix(tt_FD[,"mor_hor_est"], ncol = nrow(com_fd), byrow = T)
dis_com_exp_q2est_fd <- (1 - sim_com_exp_q2est_fd)
row.names(dis_com_exp_q2est_fd)<-row.names(com_fd)



save(dis_com_exp_q0est, dis_com_exp_q1est, dis_com_exp_q2est , file= "taxonomic_diversity/distances_com_exp_tax.rda")
save(dis_com_exp_q0est_pd, dis_com_exp_q1est_pd, dis_com_exp_q2est_pd , file= "phylogenetic_diversity/distances_com_exp_phy.rda")
save(dis_com_exp_q0est_fd, dis_com_exp_q1est_fd, dis_com_exp_q2est_fd , file= "functional_diversity/distances_com_exp_func.rda")





