
# the script changes two types of incidence data
# a data frame with species and several samples/files per plot
# a data frame with already aggregated incidence values 
# plots can have different numbers of analysed samples/files

# example data contains bird species and audio files 


rm(list=ls(all=TRUE))

library(dplyr)
library(iNEXT.3D)

## raw data with several samples/files per plot 
data <- read.table("data/ExampleBirdsEcuadorIncidences.csv", header=T, check.names = F,sep=";",encoding="latin1")

## transform into list with 1 data frame per plot 
inci.raw = lapply(unique(data$Plot), function(i) data %>% filter(Plot == i) %>% .[,-(1:2)] %>% t)
names(inci.raw) = unique(data$Plot)   ## incidence raw data for iNEXT.3D format


DataInfo3D(inci.raw, datatype = "incidence_raw")  ## Data information

inci.output = iNEXT3D(inci.raw, datatype = "incidence_raw")
ggiNEXT3D(inci.output, facet.var = "Order.q", type = 1)  ## Size-based iNEXT
ggiNEXT3D(inci.output, facet.var = "Order.q", type = 3)  ## Coverage-based iNEXT


###############################################################################

#### transforms already aggregated matrix

### column with number of samples per plot needs to be named "Samples"
### column with plot names needs to be named "Plot"


data2 <- read.table("data/ExampleBirdsEcuadorIncidencesAggregated.csv", 
                    header=T, check.names = F,sep=";",encoding="latin1")


### function to transform your data 
transform_species_matrix <- function(species_matrix) {
  # Initialize an empty list to hold the data frames for each plot
  plot_dataframes <- list()
  
  # Identify the column names for "Plot" and "Samples"
  plot_col <- which(names(species_matrix) == "Plot")
  samples_col <- which(names(species_matrix) == "Samples")
  
  # Check if the required columns exist
  if (length(plot_col) == 0 || length(samples_col) == 0) {
    stop("The species matrix must contain 'Plot' and 'Samples' columns.")
  }
  
  # Iterate over each plot (row) in the species matrix
  for (i in 1:nrow(species_matrix)) {
    # Get the number of samples for the current plot
    num_samples <- species_matrix[i, samples_col]
    
    # Create an empty data frame for the current plot
    plot_df <- data.frame(matrix(nrow = ncol(species_matrix) - 2, ncol = num_samples, 
                                 dimnames = list(names(species_matrix)[-c(plot_col, samples_col)], NULL)))
    
    # Fill the data frame based on frequency values
    species_cols <- setdiff(1:ncol(species_matrix), c(plot_col, samples_col))
    for (j in species_cols) {
      frequency <- species_matrix[i, j]
      # Create a sample vector with 1s for the frequency and 0s for the rest
      if (frequency > num_samples) {
        stop("Frequency cannot exceed the number of samples.")
      }
      samples <- c(rep(1, frequency), rep(0, num_samples - frequency))
      samples <- sample(samples)  # Shuffle to randomize the placement of 1s and 0s
      
      # Assign the samples to the data frame
      plot_df[rownames(plot_df) == names(species_matrix)[j], ] <- samples
    }
    
    # Append the data frame for the current plot to the list
    plot_dataframes[[i]] <- plot_df
    names(plot_dataframes)[i] <- species_matrix[i, plot_col]  # Name the list element
  }
  
  return(plot_dataframes)
}


### transform aggregated data into list 
inci.raw2 <- transform_species_matrix(data2)

## test if it worked properly
all(data2[3,3:length(names(data2))] == rowSums(inci.raw2[["APN2"]]))
all(data2[1,3:length(names(data2))] == rowSums(inci.raw2[["ACN1"]]))


inci.output2 = iNEXT3D(inci.raw2, datatype = "incidence_raw")

ggiNEXT3D(inci.output2, facet.var = "Order.q", type = 1)  ## Size-based iNEXT
ggiNEXT3D(inci.output2, facet.var = "Order.q", type = 3)  ## Coverage-based iNEXT




