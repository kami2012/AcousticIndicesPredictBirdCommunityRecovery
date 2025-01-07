rm(list=ls())

library(ape)
library(picante)
library(vegan)
library(dplyr)


## ================================================================================================ ##
## Prepare bird incidence data 

### loading bird incidence data with the following format: 
 
##  filename | plot | bird1 | bird2 |...
##                  |  0    |  1    |... 

data_complete <- read.table("data/detections_freile_gelis_2021_2022_birds_dummy_pivot_reordered.csv", header=T, check.names = F,sep=";")

# Remove plots with all zero values (ACN1) and plots without context information (CR16)
data_filtered <- data_complete[!data_complete$plot %in% c("ACN1", "CR16"), ] # => 85 plots



####################################################################
#####################Taxonomic Diversity############################
####################################################################

com_td <- data_filtered

write.csv(com_td, file = "taxonomic_diversity/com_td.csv", row.names = FALSE, fileEncoding = "UTF-8")


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
# reorder the columns of com_pd to match the order of the tree tip labels, but keep the first two 
com_pd <- com_pd[, c(1:2, match(tr$tip.label, names(com_pd)[-c(1:2)]) + 2)]
# check if all column names of com_pd (except plot and filename) now exactly match the tip labels of the phylogenetic tree tr
all(names(com_pd)[-c(1, 2)] == tr$tip.label)


write.csv(com_pd, file = "phylogenetic_diversity/com_pd.csv", row.names = FALSE, fileEncoding = "UTF-8")
save(tr, file = "phylogenetic_diversity/bird_tree.rda")

####################################################################
#####################Functional Diversity###########################
####################################################################

com_fd <- com_td

##### load bird traits dataframe

traits <- read.csv("data/CommunityAllplots334BirdSpecies_TraitsAvonet.csv", sep=";")

# check which column names in com_fd are not present in the Name.Phylogeny.Latin column of traits
names(com_fd)[!names(com_fd) %in% traits$Name.Phylogeny.Latin]
# reorder the columns of com_fd to match the order of traits
com_fd <- com_fd[, c(1:2, match(traits$Name.Phylogeny.Latin, names(com_fd)[-c(1:2)]) + 2)]

##### clean up traits 
row.names(traits) <- traits$Name.Phylogeny.Latin

# Create a copy of traits to work on
traits_cleaned <- traits

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



write.csv(com_fd, file = "functional_diversity/com_fd.csv", row.names = FALSE, fileEncoding = "UTF-8")
save(distM, file = "functional_diversity/trait_matrix.rda")
