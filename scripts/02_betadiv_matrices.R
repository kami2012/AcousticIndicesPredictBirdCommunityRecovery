# DivPair Demo
# by Joerg and Anne
# updated by Oliver, 15.04.2024
rm(list=ls())

# Install this version of tidytree to get rid of lots of phylogenetic tree warnings
# library(devtools)
# install_version("tidytree", version = "0.4.2", repos = "http://cran.us.r-project.org")
library(tidytree)

library(iNEXT.beta3D)
library(snowfall)
library(dplyr)

# Load script by Anne Chao (updated by Oliver)
source("scripts/div_pair_coverage_OM.R")

## ================================================ Example birds ================================================ ##
## Use abundance data "Birds" that Joreg sent as an example.
## There are 177 species (rows) and 5 assemblages (columns) in "Birds" data.
## The code below is to compute observed and coverage-based similarity between any two assemblages.

# Load "Birds" data
load(file = "data/bird_matrix_taxonomic.rda")
birds1 <- t(birds)

# Load tree
load(file = "data/bird_tree.rda")

# Load traits and rearrange data
load(file = "data/bird_matrix_functional.rda")
load(file = "data/trait_matrix.rda")

birds_avo1 <- t(birds_avo)

# Compute Cmax in the joint assemblages for alpha reference samples by the function "DataInfobeta3D" in the package "iNEXT.beta3D"
# For Birds data, Cmax_joint = 0.9103252

inextinfo <- DataInfobeta3D(birds1, diversity = "TD", datatype = "abundance") 
hist(inextinfo[1:86,"SC(n)"])
abline(v=0.75, col="red")

hist(inextinfo[1:86,"SC(2n)"])
abline(v=0.75, col="red")


Cmax_joint = sapply(1:ncol(birds1), function(i) 
  sapply(1:ncol(birds1), function(j) DataInfobeta3D( birds1[,c(i,j)] )$`SC(2n)`[4])) %>% min  

## Use function iNEXTbeta3D_pair to compute similarity index.
tt_TD  = iNEXTbeta3D_pair2(birds, div0="TD", SC = 0.83)
save(tt_TD, file= "data/tt_tax.rda")

tt_PD  = iNEXTbeta3D_pair2(birds, div0="PD", SC = 0.83, PDTree0 = tr)
save(tt_PD, file= "data/tt_phy.rda")

tt_FD  = iNEXTbeta3D_pair2(birds_avo, div0="FD", SC = 0.83, FDdistM0 = distM)
save(tt_FD, file= "data/tt_func.rda")


save(tt_TD, tt_PD, tt_FD , file= "data/tt_org.rda")
#load("data/tt_org.rda")


################################################################################
## Use function iNEXTbeta3D_pair to compute similarity index.
tt_TD93  = iNEXTbeta3D_pair2(birds, div0="TD", SC = 0.93)
save(tt_TD93, file= "data/tt_tax93.rda")

tt_PD93  = iNEXTbeta3D_pair2(birds, div0="PD", SC = 0.93, PDTree0 = tr)
save(tt_PD93, file= "data/tt_phy93.rda")

tt_FD93  = iNEXTbeta3D_pair2(birds_avo, div0="FD", SC = 0.93, FDdistM0 = distM)
save(tt_FD93, file= "data/tt_func93.rda")


save(tt_TD93, tt_PD93, tt_FD93 , file= "data/tt_org93.rda")
#load("data/tt_org93.rda")


colnames(tt_TD) # jac, sor - > q = 0, hor -> q = 1, mor_hor -> q = 2

# Results TD
dis_birds_q0est <- matrix(tt_TD[,"sor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q0est <- (1 - dis_birds_q0est)
row.names(dis_birds_q0est)<-row.names(birds)

dis_birds_q1est <- matrix(tt_TD[,"hor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q1est <- (1 - dis_birds_q1est)
row.names(dis_birds_q1est)<-row.names(birds)

dis_birds_q2est <- matrix(tt_TD[,"mor_hor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q2est <- (1 - dis_birds_q2est)
row.names(dis_birds_q2est)<-row.names(birds)



# Results PD
dis_birds_q0est_pd <- matrix(tt_PD[,"sor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q0est_pd <- (1 - dis_birds_q0est_pd)
row.names(dis_birds_q0est_pd)<-row.names(birds)

dis_birds_q1est_pd <- matrix(tt_PD[,"hor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q1est_pd <- (1 - dis_birds_q1est_pd)
row.names(dis_birds_q1est_pd)<-row.names(birds)

dis_birds_q2est_pd <- matrix(tt_PD[,"mor_hor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q2est_pd <- (1 - dis_birds_q2est_pd)
row.names(dis_birds_q2est_pd)<-row.names(birds)



# Results FD
dis_birds_q0est_fd <- matrix(tt_FD[,"sor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q0est_fd <- (1 - dis_birds_q0est_fd)
row.names(dis_birds_q0est_fd)<-row.names(birds)

dis_birds_q1est_fd <- matrix(tt_FD[,"hor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q1est_fd <- (1 - dis_birds_q1est_fd)
row.names(dis_birds_q1est_fd)<-row.names(birds)

dis_birds_q2est_fd <- matrix(tt_FD[,"mor_hor_est"], ncol = nrow(birds), byrow = T)
dis_birds_q2est_fd <- (1 - dis_birds_q2est_fd)
row.names(dis_birds_q2est_fd)<-row.names(birds)



save(dis_birds_q0est, dis_birds_q1est, dis_birds_q2est , file= "data/distances_tax83.rda")
save(dis_birds_q0est_pd, dis_birds_q1est_pd, dis_birds_q2est_pd , file= "data/distances_phy83.rda")
save(dis_birds_q0est_fd, dis_birds_q1est_fd, dis_birds_q2est_fd , file= "data/distances_func83.rda")





