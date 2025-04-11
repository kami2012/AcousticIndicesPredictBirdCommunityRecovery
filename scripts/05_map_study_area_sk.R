rm(list = ls())

#install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")

library(scales)
library(maptools)
library(raster)
library(sf)
library(dplyr)
library(rworldmap)
library(mapview)
library(maptiles)
library(prettymapr)

# Working directory
setwd("E:/Manuscripts/1_3D_ecuador")


# Read the plot information
plc <- read.table("data/BirdCommunityNamePhylogeny_Indices.csv", header=T, sep = ";", dec = ",")

# Keep only the plot information
plc <- plc[, 1:6]

# Creating a new data frame with color and point shapes for every plot category
clr_data <- data.frame(
  clr = c("orange", "orange", "sienna", "#ffd500", "#ffd500", "#92cb11", "#92cb11", "#359619", "#359619", "darkgreen"),
  Category10 = c("A_Caca", "A_Past", "A_Old", "F_Caca", "F_Past", "F_CReg1", "F_PReg1", "F_CReg2", "F_PReg2", "F_Old"),
  shp = c(24, 22, 21, 24, 22, 24, 22, 24, 22, 21), # pch values for symbols with borders (24=triangle, 22=square, 21=circle)
  border = "white"  # White border for all shapes
)


# Add the color and point shape to all plots
plc$clrs <- clr_data[match(plc$Category10, clr_data$Category10), "clr"]
plc$shp <- clr_data[match(plc$Category10, clr_data$Category10), "shp"]
plc$border <- clr_data[match(plc$Category10, clr_data$Category10), "border"]


# Convert plc into a spatial object (sf); Assigns the Coordinate Reference System (CRS) as EPSG:4326 (WGS 84)
plc <- st_as_sf(plc, coords = c("Longitude","Latitude"), crs=4326)

# Preview the data
mapview(plc)


# Define the area of interest correctly
bboxsq <- st_bbox(st_buffer(st_centroid(st_combine(plc)), 12000), crs=4326)

# Retrieve the background map
tiles <- get_tiles(x = bboxsq, provider = "Esri.WorldImagery", crop = TRUE, 
                   cachedir = tempdir(), verbose = TRUE, zoom = 13, project = TRUE)  


# Get a map from Ecuador to locate the study area
df <- getMap(resolution = "low")
newmap <- df


# Save figure
tiff("plots/studyarea_ecuador_satellite.tiff", width = 220, height = 220, units = "mm", res = 500, compression = "lzw")

par(mar = c(0.1, 0.1, 0.1, 0.1))

# Create the base plot first
plot_tiles(tiles, adjust = FALSE, smooth = TRUE)

# Plot points with colored fill and white border
plot(st_geometry(plc), 
     pch = plc$shp,        # Shape from your data (21 for circle, 22 for square, 24 for triangle)
     col = plc$border,     # Border color (white)
     bg = plc$clrs,        # Fill color
     cex = 1.5,            # Point size
     add = TRUE)           # Add to existing plot


box()



# Add the scale bar
barpar <- scalebarparams()
barpart <- (barpar$widthplotunit / barpar$widthu) / (par("usr")[2] - par("usr")[1])
barbar <- 1 * (barpart * 10)

par(new = T)
plot(1:10, xlab = "", ylab = "", type = "n", axes = F)
segments(1.25, 1.25, 1.25 + barbar, 1.25, col = "white")
text(1.25 + (barbar / 2), 1.65, "1 km", col = "white")

# Add a legend
legend("bottomright",
       pch = c(NA, 19, 19, 19, 
               NA, 19, 19, 19, 19, 19, 
               2, 0, 1),
       col = c("black", "orange", "sienna", "white", 
               "black", "#ffd500", "#92cb11", "#359619", "darkgreen", "white",
               "black", "black", "black"),
       legend = c(expression(bold("Agricultural Matrix")), "Agriculture", "Old-growth", "",
                  expression(bold("Forest Matrix")), "Agriculture" , "Regeneration I", "Regeneration II", "Old-growth", "",
                  "(former) Cacao", "(former) Pasture", "Old-growth"),
       cex = 1,2,  # Text size
       bty = "o",  # Box type: "o" for box, "n" for none
       bg = "white")  # Background color

# Add location map
par(new = TRUE, fig = c(0, 0.35, 0.65, 1), mar = c(0, 0, 0, 0))
plot(newmap, 
     xlim = c(-100, -30), 
     ylim = c(-30, 23), 
     col = "#ededed", 
     border = "white")

# Add a circle at the exact position of the study area
plot(st_geometry(st_centroid(st_combine(plc))), 
     col = "darkgreen", 
     cex = 1.6, 
     add = TRUE)

box()

dev.off()




























