# DivPair Example, version 2
# 2024-11-04

rm(list=ls(all=TRUE))

library(iNEXT.beta3D)
library(snowfall)
library(dplyr)
library(vegan)

source(file = "scripts/DivPairCoverage_V2.R")

load(file = "data/bird_data_incidence_phylo.rda")

#### TD with parallel cpus
pairwise_TD <- iNEXTbeta3D_pair3D(data, SC = 0.9, div0="TD", datatype0 = "incidence_raw", parallel = T, cpus = 4)

range(pairwise_TD$Matrices$jac_obs)
range(pairwise_TD$Matrices$sor_obs)

### now we have negative values, but they should be between 0 and 1 shouldnt they? 


#### PD with parallel cpus
pairwise_PD <- iNEXTbeta3D_pair3D(data, SC = 0.9, div0="PD",PDTree0= tr, datatype0 = "incidence_raw", parallel = T, cpus = 6)








