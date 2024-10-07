
rm(list=ls(all=TRUE))
Sys.setenv(LANG = "en")

library(ape)
library(picante)
library(vegan)
library(dplyr)


#### loading species data and preparing with correct names for phylogenetic analysis 
all_G <- read.table("data/CommunityAllplots334BirdSpecies.csv", header=T, check.names = F,sep=";",encoding="latin1")

names(all_G)
birds <- all_G
row.names(birds) <- birds$`Name Phylogeny Latin`

birds$Avonet_Latin <- NULL
birds$Latin_correct <- NULL
birds$`Name Phylogeny Latin` <- NULL
birds<-as.data.frame(t(birds))



#reading environmental data for 66 plots
environ<-read.csv2("data/EnvironAll86plots.csv", header=T)
#environ$Category7_ord<-ordered(environ$Category7, c("Pasture","Cacao","PRegI","CRegI","PRegII","CRegII","Oldgr"))


#matching community and environmental data
birds<-birds[as.vector(environ$Plot),]

#adjusting format of species names
#names(birds) <- gsub(" ","_",names(birds), fixed = TRUE)




##### load phylo tree and compare with bird data 

tr <- read.tree(file = "data/TreeBirds334Canande.tre")

names(birds)[!names(birds)%in%tr$tip.label]
tr <- drop.tip(tr, setdiff(tr$tip.label, names(birds)))

birds <- birds[,match(tr$tip.label, names(birds))]
all(names(birds) == tr$tip.label)


save(birds, file = "data/bird_matrix_taxonomic.rda")
save(tr, file = "data/bird_tree.rda")
#################################################################################
#################################################################################

###### species matrix with avonet names for functional analysis

birds_avo <- all_G
row.names(birds_avo) <- birds_avo$Avonet_Latin

birds_avo$Avonet_Latin <- NULL
birds_avo$Latin_correct <- NULL
birds_avo$`Name Phylogeny Latin` <- NULL
birds_avo<-as.data.frame(t(birds_avo))

#matching community and environmental data
birds_avo<-birds_avo[as.vector(environ$Plot),]



#### load traits and compare species 
traits <- read.table("data/CommunityAllplots334BirdSpecies_TraitsAvonet.csv", header=T, check.names = F,sep=";",dec = ",",encoding="latin1")
names(birds_avo)[names(birds_avo)%in%traits$Avonet_Latin]


##### clean up traits 
row.names(traits) <- traits$Avonet_Latin

####### model body size dependent traits and take residuals 
####### body mass in log(bodymass) 

traits$log_mass <- log(traits$Mass)

mod.beak.l <- lm(log(BeakLength_Culmen) ~ log(Mass), data = traits)
traits$beak_length_res <- mod.beak.l$residuals

mod.beak.w <- lm(log(BeakWidth) ~ log(Mass), data = traits)
traits$beak_width_res <- mod.beak.w$residuals

mod.tarsus <- lm(log(TarsusLength) ~ log(Mass), data = traits)
traits$tarsus_res <- mod.tarsus$residuals

mod.wing <- lm(log(WingLength) ~ log(Mass), data = traits)
traits$wing_res <- mod.wing$residuals

mod.tail <- lm(log(TailLength) ~ log(Mass), data = traits)
traits$tail_res <- mod.tail$residuals


keep <- c("Hand-WingIndex","HabitatDensity","Migration","TrophicNiche","PrimaryLifestyle",
          "beak_length_res","beak_width_res","tarsus_res","wing_res","tail_res","log_mass")


traits <- traits[,names(traits) %in% keep]



for (i in 1:ncol(traits)) {
  if (class(traits[,i]) == "character") traits[,i] <- factor(traits[,i], levels = unique(traits[,i]))
}

distM <- cluster::daisy(x = traits, metric = "gower") %>% as.matrix()




save(birds_avo, file = "data/bird_matrix_functional.rda")
save(distM, file = "data/trait_matrix.rda")
