
# the script converts pivot incidence data to "incidence_raw" (= a list of matrices/data.frames, with each matrix
# representing a species-by-sampling-unit incidence matrix for one of the assemblages)


rm(list=ls(all=TRUE))

library(dplyr)
library(iNEXT.3D)

## raw data 
## plot | filename | bird1 | bird2 |...
##                 |  0    |  1    |... 

data <- read.table("data/detections_freile_gelis_2021_2022_birds_dummy_pivot.csv", header=T, check.names = F,sep=";",encoding="latin1")

## transform into list with 1 data frame (rows = Files, columns = Species) per Plot 
inci.raw = lapply(unique(data$plot), function(i) data %>% filter(plot == i) %>% .[,-(1:2)] %>% t)
names(inci.raw) = unique(data$plot)   ## incidence raw data for iNEXT.3D format

# Check and remove plots with all zero values (=> ACN1)
inci.raw <- inci.raw[sapply(inci.raw, function(x) any(x != 0))]

# Save inci.raw as a RDS file
saveRDS(inci.raw, file = "data/inci_raw.rds")


