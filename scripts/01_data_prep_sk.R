
rm(list=ls(all=TRUE))
Sys.setenv(LANG = "en")

library(ape)
library(picante)
library(vegan)
library(dplyr)


#### loading expert community composition data; take plots as rownames 
data_com <- read.csv2("data/BirdCommunityNamePhylogeny_Indices.csv", row.names=1,header=T)

rownames(data_com)
names(data_com)

####################################################################
#####################Taxonomic Diversity############################
####################################################################

com_td<-data_com[,14:347]
# 85 334
dim(com_td)

write.csv(com_td, file = "taxonomic_diversity/com_td.csv")

####################################################################
#####################Phylogenetic Diversity#########################
####################################################################

com_pd <- com_td

##### load phylo tree and compare with data_com 
tr <- read.tree(file = "data/TreeBirds334Canande.tre")

# check which column names in com_pd are not present in the tip labels of the phylogenetic tree tr.
names(com_pd)[!names(com_pd)%in%tr$tip.label]
# find and drop tips (species) in the tree that are not in the com_pd 
tr <- drop.tip(tr, setdiff(tr$tip.label, names(com_pd)))
# reorder the columns of com_pd to match the order of the tree tip labels
com_pd <- com_pd[,match(tr$tip.label, names(com_pd))]
# check if all column names of com_pd now exactly match the tip labels of the phylogenetic tree tr
all(names(com_pd) == tr$tip.label)


write.csv(com_pd, file = "phylogenetic_diversity/com_pd.csv")
save(tr, file = "phylogenetic_diversity/bird_tree.rda")

####################################################################
#####################Functional Diversity###########################
####################################################################

com_fd <- com_td

##### load bird traits dataframe

data_traits <- read.csv("data/CommunityAllplots334BirdSpecies_TraitsAvonet.csv", sep=";")

# check which column names in com_fd are not present in the Name.Phylogeny.Latin column of data_traits
names(com_fd)[!names(com_fd) %in% data_traits$Name.Phylogeny.Latin]
# reorder the columns of com_fd to match the order of data_traits
com_fd <- com_fd[,match(data_traits$Name.Phylogeny.Latin, names(com_fd))]

##### clean up traits 
row.names(data_traits) <- data_traits$Name.Phylogeny.Latin

# Create a copy of the dataframe to work on
traits_cleaned <- data_traits

# List of columns that should be numeric
numeric_cols <- c("BeakLength_Culmen", "BeakLength_Nares", "BeakWidth", "BeakDepth", 
                  "TarsusLength", "WingLength", "KippsDistance", "Secondary1", 
                  "Hand.WingIndex", "TailLength", "Mass", "MinLatitude", 
                  "MaxLatitude", "CentroidLatitude", "CentroidLongitude", "RangeSize")

# Function to convert columns to numeric
convert_to_numeric <- function(column) {
  # Replace commas with dots
  column <- gsub(",", ".", column)
  # Convert to numeric
  as.numeric(column)
}

# Apply conversion to specified columns
traits_cleaned[numeric_cols] <- lapply(traits_cleaned[numeric_cols], convert_to_numeric)

####### model body size dependent traits and take residuals 
####### body mass in log(bodymass) 

traits <- traits_cleaned

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



write.csv(com_fd, file = "functional_diversity/com_fd.csv")
save(distM, file = "functional_diversity/trait_matrix.rda")
