rm(list=ls())

library(vegan)

## ================================================================================================ ##
## Use abundance data / community composition "BirdCommunityNamePhylogeny_Indices.csv" from Joerg.
## There are 85 assemblages (rows) and 334 birds (columns), identified by experts.
## The table contains also 8 different sound indices and plot categories (Pasture/Cacao -> Old Growth)
## The table contains also the first and second nmds Axis for td, pd, and fd (q0, q1, q2)

## Creating the distance by sound

# load the dataframe with all information (plots, sound indices, bird abundance data, nmds axis 1 and 2)
data_com <- read.csv("data/BirdCommunityNamePhylogeny_Indices_nmds.csv", row.names = 1) 

# load distance_com_exp data for consistency check
load("taxonomic_diversity/distances_com_exp_tax.rda")

# Store the sound indices columns 
# (AcousticDiversity, AcousticEveness, BioAcoustic, SoundscapeSaturation,
# EntropyOfVarianceSpectrum, AcousticComplexity, TemporalEntropy, EventsPerSecond)
sound_indices_columns <- names(data_com)[6:13]

# Create a list to store the sound index distance matrices
distances_soundIndices <- list()

# Loop through each sound index
for (sound_index in sound_indices_columns) {
  
  # Convert column to dataframe
  SoundIndex <- as.data.frame(data_com[[sound_index]])
  
  # Set Plots as row names
  rownames(SoundIndex) <- rownames(data_com)
  
  # Standardize (subtracting the mean and dividing by the standard deviation)
  SoundIndex_s <- decostand(SoundIndex, "standardize")
  
  # Compute Euclidean Distance Matrix (Euclidean distance measures the straight-line distance between points in multidimensional space.)
  SoundIndex_s_dist <- dist(SoundIndex_s, method = "euclidean")
  
  # Convert the distance matrix to a data frame ...
  tt <- data.frame(as.matrix(SoundIndex_s_dist))
  
  # ... reorder it to ensure consistency ...
  tt <- tt[(row.names(tt)),]
  tt <- tt[,(colnames(tt))] 
  
  # ... and then convert it back to a distance object
  dist_name <- paste0("dis_", sound_index, "_s")
  distances_soundIndices[[dist_name]] <- as.dist(tt, upper = TRUE, diag = TRUE)
  
  # convert dis_com_exp to distance matrix
  dis_com_exp_q0est <- as.dist(dis_com_exp_q0est, upper = TRUE, diag = TRUE)
  
  # check consistency
  print(all(rownames(as.matrix(distances_soundIndices[[dist_name]])) == rownames(as.matrix(dis_com_exp_q0est))))
}
  
save(distances_soundIndices, file = "data/distances_soundIndices.RData")





