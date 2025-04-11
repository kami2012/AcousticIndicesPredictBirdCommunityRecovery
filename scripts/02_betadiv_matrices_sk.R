rm(list=ls())

# Install this version of tidytree to get rid of lots of phylogenetic tree warnings
# library(devtools)
# install_version("tidytree", version = "0.4.2", repos = "http://cran.us.r-project.org")
library(tidytree)
library(iNEXT.beta3D)
library(snowfall)
library(dplyr)

# Working directory
setwd("E:/Manuscripts/1_3D_ecuador")

# DivPair Demo by Joerg and Anne
# updated by Oliver, 15.04.2024
source("scripts/DivPairCoverage_V2.R")

## ================================================================================================ ##
## Use incidence data / community composition
## The code below is to compute observed and coverage-based/estimated dissimilarity between any two assemblages.


# Load taxonomic diversity data
com_td <- read.csv("taxonomic_diversity/com_td.csv", row.names = NULL)

# Load phylogenetic diversity data
com_pd <- read.csv("phylogenetic_diversity/com_pd.csv", row.names = NULL)
load(file = "phylogenetic_diversity/bird_tree.rda")

# Load functional diversity data
com_fd <- read.csv("functional_diversity/com_fd.csv", row.names = NULL)
load(file = "functional_diversity/trait_matrix.rda")


########################## Compute Sampling Coverage (SC) ############################

### Use function "DataInfobeta3D" in the package "iNEXT.beta3D"
## Change the data format of com_td to fit to iNEXT.3D 
# Transform into list with 1 data frame (rows = Files, columns = Birds) per Plot 
birds_inci = lapply(unique(com_td$plot), function(i) com_td %>% filter(plot == i) %>% .[,-(1:2)] %>% t)
names(birds_inci) = unique(com_td$plot)   ## incidence raw data for iNEXT.3D format

#info <- DataInfo3D(birds_inci, diversity = 'TD', datatype = "incidence_raw")

info <- DataInfobeta3D(birds_inci, diversity = "TD", datatype = "incidence_raw",
               PDtree = NULL, PDreftime = NULL, FDdistM = NULL, FDtype = "AUC", FDtau = NULL) 

## Take the median of the extrapolated Sampling Coverages (2T) as Sampling Coverage (SC)
SC <- median(info$`SC(2T)`) # 0.913082021...
#SC <- 0.913082021

# Histogram of the extrapolated Sampling Coverages (2T)
# Open a PNG graphics device

# Supplementary Figure S1 

png("plots/histogram_sc.png")

hist(info[1:85, "SC(2T)"], 
     main = "Extrapolated Sampling Coverages (SC(2T)) of the 85 plots", 
     xlab = "SC(2T)", 
     col = "lightblue", 
     border = "black")

# Line for the chosen Sampling Coverage (SC)
abline(v = SC, col = "blue", lwd = 2, lty = 2)

# Label the SC line
text(x = SC, y = par("usr")[4] * 0.8, labels = "median = SC", col = "blue", pos = 4)

dev.off()


########################## Observed (with default SC) ################################
# !!!!!!!!!!!!!! NOT USED !!!!!!!!!!!!!!!!!

## Use function iNEXTbeta3D_pair to compute similarity index (with default (min) sample coverage (SC))
pairwise_TD  = iNEXTbeta3D_pair3D(com_td, div0="TD", SC = NULL, datatype0="incidence_raw", parallel = T, cpus = 4)
save(pairwise_TD, file= "taxonomic_diversity/tt_tax.rda")

pairwise_PD  = iNEXTbeta3D_pair3D(com_pd, div0="PD", SC = NULL, datatype="incidence_raw", PDTree0 = tr, 
                            parallel = T, cpus = 4)
save(pairwise_PD, file= "phylogenetic_diversity/tt_phy.rda")

pairwise_FD  = iNEXTbeta3D_pair3D(com_fd, div0="FD", SC = NULL, datatype="incidence_raw", FDdistM0 = distM)
save(pairwise_FD, file= "functional_diversity/tt_func.rda")


save(pairwise_TD, pairwise_PD, pairwise_FD , file= "data/tt_org.rda")
#load("data/tt_org.rda")

# !!!!!!!!!!!!!! NOT USED !!!!!!!!!!!!!!!!!



########################## Coverage-based (with computed SC) ################################

## Use function iNEXTbeta3D_pair3D AND sample coverage to compute similarity index.
 
pairwise_TD91 = iNEXTbeta3D_pair3D(com_td, div0="TD", SC = SC, datatype0 = "incidence_raw", parallel = T, cpus = 8)
save(pairwise_TD91, file= "taxonomic_diversity/pairwise_TD91.rda")


pairwise_PD91 = iNEXTbeta3D_pair3D(com_pd, div0="PD", SC = SC, datatype0 = "incidence_raw", PDTree0 = tr, 
                                   parallel = T, cpus = 8)
save(pairwise_PD91, file= "phylogenetic_diversity/pairwise_PD91.rda")


pairwise_FD91 = iNEXTbeta3D_pair3D(com_fd, div0="FD", SC = SC, datatype0 = "incidence_raw", FDdistM0 = distM, 
                                  parallel = T, cpus = 8)
save(pairwise_FD91, file= "functional_diversity/pairwise_FD91.rda")


save(pairwise_TD91, pairwise_PD91, pairwise_FD91 , file= "data/pairwise_TD_PD_FD.rda")
#load("data/pairwise_TD_PD_FD.rda")




## Compute dissimilarity matrices for the different diversities (taxonomic, pyhlogenetic and functional) and Hill numbers (q0, q1, q2)

#  Dissimilarity matrices TD 
## observed
dis_com_q0obs_td <- 1- pairwise_TD91[["Matrices"]][["jac_obs"]]
dis_com_q1obs_td <- 1- pairwise_TD91[["Matrices"]][["hor_obs"]]
dis_com_q2obs_td <- 1- pairwise_TD91[["Matrices"]][["mor_hor_obs"]]

## estimated
dis_com_q0est <- 1- pairwise_TD91[["Matrices"]][["jac_est"]]
dis_com_q1est <- 1- pairwise_TD91[["Matrices"]][["hor_est"]]
dis_com_q2est <- 1- pairwise_TD91[["Matrices"]][["mor_hor_est"]]



#  Dissimilarity matrices PD
## observed
dis_com_q0obs_pd <- 1- pairwise_PD91[["Matrices"]][["jac_obs"]]
dis_com_q1obs_pd <- 1- pairwise_PD91[["Matrices"]][["hor_obs"]]
dis_com_q2obs_pd <- 1- pairwise_PD91[["Matrices"]][["mor_hor_obs"]]

## estimated
dis_com_q0est_pd <- 1- pairwise_PD91[["Matrices"]][["jac_est"]]
dis_com_q1est_pd <- 1- pairwise_PD91[["Matrices"]][["hor_est"]]
dis_com_q2est_pd <- 1- pairwise_PD91[["Matrices"]][["mor_hor_est"]]



#  Dissimilarity matrices FD
##observed
dis_com_q0obs_fd <- 1- pairwise_FD91[["Matrices"]][["jac_obs"]]
dis_com_q1obs_fd <- 1- pairwise_FD91[["Matrices"]][["hor_obs"]]
dis_com_q2obs_fd <- 1- pairwise_FD91[["Matrices"]][["mor_hor_obs"]]

## estimated
dis_com_q0est_fd <- 1- pairwise_FD91[["Matrices"]][["jac_est"]]
dis_com_q1est_fd <- 1- pairwise_FD91[["Matrices"]][["hor_est"]]
dis_com_q2est_fd <- 1- pairwise_FD91[["Matrices"]][["mor_hor_est"]]




save(dis_com_q0est, dis_com_q1est, dis_com_q2est , file= "taxonomic_diversity/distances_com_exp_tax.rda")
save(dis_com_q0est_pd, dis_com_q1est_pd, dis_com_q2est_pd , file= "phylogenetic_diversity/distances_com_exp_phy.rda")
save(dis_com_q0est_fd, dis_com_q1est_fd, dis_com_q2est_fd , file= "functional_diversity/distances_com_exp_func.rda")





