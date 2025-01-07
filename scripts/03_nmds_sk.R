rm(list=ls())

library(vegan)



## ================================================================================================ ##
## Compute and plot non-metric multidimensional scaling (nmds)

## Load the dissimilarity matrices
# Load taxonomic dissimilarity matrices
load(file = "taxonomic_diversity/distances_com_exp_tax.rda")

# Load phylogenetic dissimilarity matrices
load(file = "phylogenetic_diversity/distances_com_exp_phy.rda")

# Load functional dissimilarity matrices
load(file = "functional_diversity/distances_com_exp_func.rda")


## Load the plot and category information
# Load community composition data (old!!!) containing plot, categories and sound indices
data_com_complete <- read.csv2("data/BirdCommunityNamePhylogeny_Indices.csv", 
                               header = TRUE)
# Use only plot, categories and sound indices
data_com <- data_com_complete[, 1:14]

# Reorder data_com to match the rownames in dis_com_q0est
data_com <- data_com[match(rownames(dis_com_q0est), data_com$Plot_ID), ]


######################################################################################################
######################## Perform Nonmetric Multidimensional Scaling (NMDS) ###########################

# List of dissimilarity matrices
dissimilarity_matrices <- list(
  TD_q0 = dis_com_q0est,
  TD_q1 = dis_com_q1est,
  TD_q2 = dis_com_q2est,
  FD_q0 = dis_com_q0est_fd,
  FD_q1 = dis_com_q1est_fd,
  FD_q2 = dis_com_q2est_fd,
  PD_q0 = dis_com_q0est_pd,
  PD_q1 = dis_com_q1est_pd,
  PD_q2 = dis_com_q2est_pd
)

# Initialize a list to store NMDS results
nmds_results <- list()

# Initialize a dataframe to store stress levels
stress_levels <- data.frame(
  Matrix = character(),
  Stress = numeric(),
  stringsAsFactors = FALSE
)

# Compute NMDS for each matrix
for (name in names(dissimilarity_matrices)) {
  
  # Adjust the value of k based on the matrix name
  #k <- if (name %in% c("TD_q0", "PD_q0")) 3 else 2
  
  set.seed(1)
  nmds <- metaMDS(dissimilarity_matrices[[name]], k = 2, maxit = 999, trymax = 100)
  nmds_results[[name]] <- nmds
  
  # Save axis scores to the data_com dataframe
  axis1 <- paste0(name, "_Axis1")
  data_com[[axis1]] <- nmds$points[, 1]
  
  # Save the stress level
  stress_levels <- rbind(stress_levels, data.frame(Matrix = name, Stress = nmds$stress))
}


# Save the NMDS results to a file
saveRDS(nmds_results, "data/nmds_results.rds")

# Save the updated data_com dataframe with axis scores
write.csv(data_com, file = "data/plots_categories_indices_nmds.csv", row.names = FALSE)

# Save the stress levels to a CSV file
write.csv(stress_levels, file = "data/nmds_stress_levels.csv", row.names = FALSE)


##########################################################################################
################################ Plot nmds ###############################################

################################ Plot all in one #########################################

# Load NMDS results
nmds_results <- readRDS("data/nmds_results.rds")

# Manual color mapping
color_mapping <- c(
  "A_Caca" = "orange",
  "A_Past" = "orange",
  "A_Old" = "sienna",
  "F_Caca" = "yellow",
  "F_Past" = "yellow",
  "F_CReg1" = "greenyellow",
  "F_PReg1" = "greenyellow",
  "F_CReg2" = "chartreuse3",
  "F_PReg2" = "chartreuse3",
  "F_Old" = "darkgreen"
)

data_com$Color <- color_mapping[data_com$Category10]

# Helper function to plot NMDS, add hulls and spiders
plot_nmds <- function(nmds, name) {
  ordiplot(nmds, display = "sites", type = "n")
  
  # Plot points with their respective colors
  points(nmds, display = "sites", col = data_com$Color, pch = 19, cex = 0.8)
  
  # Add hulls for each group with correct color mapping
  unique_groups <- unique(data_com$Category10)
  for (group in unique_groups) {
    ordihull(nmds, groups = data_com$Category10, display = "sites", 
             show.groups = group, draw = "polygon", col = color_mapping[group], alpha = 0.5)
  }
  
  # Add spiders for each group with correct color mapping
  for (group in unique_groups) {
    ordispider(nmds, groups = data_com$Category10, display = "sites", 
               show.groups = group, col = color_mapping[group], label = TRUE)
  }
}

# Create the plot layout
png(filename = "plots/nmds_new.png", width = 170, height = 200, units = "mm", res = 1000)
layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol = 3, byrow = TRUE), 
       widths = c(1, 1, 1), heights = c(1, 1, 1))
par(mar = c(2, 1.5, 0.5, 1.5), oma = c(1.5, 4, 3, 4))

# Plot each NMDS result
plot_order <- c("TD_q0", "TD_q1", "TD_q2", "FD_q0", "FD_q1", "FD_q2", "PD_q0", "PD_q1", "PD_q2")
for (name in plot_order) {
  plot_nmds(nmds_results[[name]], name)
}

# Add labels for rows and columns
mtext("q = 0", side = 3, line = 0.5, outer = TRUE, at = 0.167, cex = 1.2)
mtext("q = 1", side = 3, line = 0.5, outer = TRUE, at = 0.5, cex = 1.2)
mtext("q = 2", side = 3, line = 0.5, outer = TRUE, at = 0.833, cex = 1.2)

mtext("Taxonomic Diversity", side = 2, line = 2, outer = TRUE, at = 0.83, cex = 1.2)
mtext("Functional Diversity", side = 2, line = 2, outer = TRUE, at = 0.5, cex = 1.2)
mtext("Phylogenetic Diversity", side = 2, line = 2, outer = TRUE, at = 0.17, cex = 1.2)

dev.off()

################################ Plot only  PD q = 0 #########################################

# Select the right nmds object
nmds <- nmds_results[["PD_q0"]]

# Create the plot layout
png(filename = "plots/nmds_PD_q0.png", width = 180, height = 200, units = "mm", res = 1000)

par(mar = c(5, 5, 4, 4)) # bottom, left, top, right

ordiplot(nmds, display = "sites", type = "n")

# Add a title
#title(main = "NMDS of Phylogenetic Diversity, q = 0", cex.main = 1.5)

# Plot points with their respective colors
points(nmds, display = "sites", col = data_com$Color, pch = 19, cex = 0.8)

# Add hulls for each group with correct color mapping
unique_groups <- unique(data_com$Category10)
for (group in unique_groups) {
  ordihull(nmds, groups = data_com$Category10, display = "sites", 
           show.groups = group, draw = "polygon", col = color_mapping[group], alpha = 0.5)
}

# Add spiders for each group with correct color mapping
for (group in unique_groups) {
  ordispider(nmds, groups = data_com$Category10, display = "sites", 
             show.groups = group, col = color_mapping[group], label = TRUE)
}

# Add a legend with adjusted inner margins
legend("topleft", 
       legend = c("A = \"Agricultural Matrix\"", "F = \"Forest Matrix\""), 
       bty = "o",          # Box around the legend
       cex = 0.9,          # Adjust text size
       text.col = "black", 
       inset = 0.02)       # Move the legend slightly inside
       

dev.off()
