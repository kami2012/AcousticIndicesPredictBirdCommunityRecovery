rm(list=ls())

library(gt)
library(webshot2)

## ================================================================================================ ##
## Test sound indices as predictors for the first axis of ordination of td, pd and fd  (q = 0, q = 1, q = 2)



######################################################################################################
######################## Define the data, response variable, and predictors ##########################


# Load data
data_com <- read.csv("data/plots_categories_indices_nmds.csv")
# Sort by nmds axis 1
data_com <- data_com[order(data_com$TD_q0_Axis1), ]

# Define the response variables
nmds_axis1 <- c("TD_q0_Axis1", "TD_q1_Axis1", "TD_q2_Axis1",
                "FD_q0_Axis1", "FD_q1_Axis1", "FD_q2_Axis1",
                "PD_q0_Axis1", "PD_q1_Axis1", "PD_q2_Axis1")

# Define the predictors 
predictors <- c("SoundscapeSaturation", "EntropyOfVarianceSpectrum", "AcousticComplexity", 
                "TemporalEntropy", "EventsPerSecond")



######################################################################################################
#################################### Predict nmds axis 1 #############################################


# Split the data into train (0.66) and test (0.33) by picking every 3rd row

# Every third row for the test set
test_data <- data_com[seq(3, nrow(data_com), by = 3), ]

# Remaining rows (those not in the test set) for the training set
train_data <- data_com[-seq(3, nrow(data_com), by = 3), ]

# Generate a grid of hypothetical values 
# for the prediction interval
new_data <- expand.grid(
  SoundscapeSaturation = seq(min(data_com$SoundscapeSaturation), max(data_com$SoundscapeSaturation), length.out = 10),
  EntropyOfVarianceSpectrum = seq(min(data_com$EntropyOfVarianceSpectrum), max(data_com$EntropyOfVarianceSpectrum), length.out = 10),
  AcousticComplexity = seq(min(data_com$AcousticComplexity), max(data_com$AcousticComplexity), length.out = 10),
  TemporalEntropy = seq(min(data_com$TemporalEntropy), max(data_com$TemporalEntropy), length.out = 10),
  EventsPerSecond = seq(min(data_com$EventsPerSecond), max(data_com$EventsPerSecond), length.out = 10)
)


# Initialize lists to store R-squared values, t-values, predictions and models
# for each response variable
models <- list()
r_squared_results <- list()
t_values_results <- list()
predictions_tdata_results <- list()
predictions_ndata_results <- list()



# Loop over each response variable 
for (axis in nmds_axis1) {
  
  # Create formula for the current response variable
  formula <- as.formula(paste(axis, "~", paste(predictors, collapse = " + ")))
  
  # Fit linear model to the training data
  model <- lm(formula, data = train_data)
  
  # Store the linear model 
  models[[axis]] <- model
  
  # Extract the summary of the model
  model_summary <- summary(model)
  
  # Get the coefficients (including t-values) from the model summary
  t_values <- model_summary$coefficients[, "t value"]
  
  # Store the t-values for this axis
  t_values_results[[axis]] <- t_values
  
  
  # Predict on the test data
  predictions_tdata <- predict(model, newdata = test_data)
  
  # Store predicted values
  predictions_tdata_results[[axis]] <- predictions_tdata
  
  
  # Predict on the hypothetical data
  predictions_ndata <- predict(model, newdata = new_data, interval = "prediction")
  
  # Store predicted values
  predictions_ndata_results[[axis]] <- predictions_ndata
  
  
  # Calculate R-squared for the predictions
  r_squared <- cor(predictions_tdata, test_data[[axis]])^2
  
  # Store the R-squared value
  r_squared_results[[axis]] <- r_squared
}


# Save models, results, and datasets to a single file
save(models, r_squared_results, t_values_results, 
     predictions_tdata_results, predictions_ndata_results,
     test_data, train_data, new_data, 
     file = "data/04_soundIndices_as_predictors/data_models_and_results.RData")


######################################################################################################
###################### Plot the acoustic indices as predictors for nmds Axis 1 #######################

load("data/04_soundIndices_as_predictors/data_models_and_results.RData")

# Create the plot layout
png(filename = "plots/soundIndicesAsPred.png", width = 170, height = 200, units = "mm", res = 1000)

layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol = 3, byrow = TRUE), 
       widths = c(1, 1, 1), heights = c(1, 1, 1))

par(mar = c(1, 1.8, 1.3, 1), oma = c(5, 4, 2.5, 5)) # bottom, left, top, right

# Loop through each NMDS axis for plotting
for (axis in nmds_axis1) {
  
  # Get predictions and observed values
  predicted_tdata_values <- predictions_tdata_results[[axis]]
  observed_values <- test_data[[axis]]
  predicted_ndata_values <- predictions_ndata_results[[axis]][, "fit"]
  
  # Scatterplot: Predicted vs Observed
  plot(predicted_tdata_values, observed_values, 
       #xlab = "Predicted Bird Community", 
       #ylab = "Observed Bird Community",
       pch = 1,
       col = rgb(0, 0, 0, 0.5))  # Black points with 50% transparency
  
  # Extract the lower and upper bounds of the prediction interval
  # based on the hyphothetical data
  lower_bound <- predictions_ndata_results[[axis]][, "lwr"]
  upper_bound <- predictions_ndata_results[[axis]][, "upr"]
  
  ordering <- order(predicted_ndata_values)
  predicted_ndata_values <- predicted_ndata_values[ordering]
  lower_bound <- lower_bound[ordering]
  upper_bound <- upper_bound[ordering]
  
  
  # Plot prediction interval
  polygon(c(predicted_ndata_values, rev(predicted_ndata_values)),
          c(upper_bound, rev(lower_bound)), 
          col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Fit regression line
  abline(lm(observed_values ~ predicted_tdata_values), col = "blue", lwd = 2)
  
  # Display R-squared value
  r_squared <- r_squared_results[[axis]]
  text_x <- min(predicted_tdata_values)  # X coordinate for the text
  text_y <- max(observed_values) - 0.05 * diff(range(observed_values)) # Y coordinate for the R² text
  text(text_x, text_y, 
       paste("R² =", round(r_squared, 2)), 
       pos=4, col="blue")
  
}

# Add labels for rows and columns
mtext("q = 0", side=3, line=0.5, outer=TRUE, at=0.167, cex=1.2)
mtext("q = 1", side=3, line=0.5, outer=TRUE, at=0.5, cex=1.2)
mtext("q = 2", side=3, line=0.5, outer=TRUE, at=0.833, cex=1.2)

# Left side labels (y-axis)
mtext("Observed Bird Community", side=2, line=2, outer=TRUE, at=0.83, cex=0.9)
mtext("Observed Bird Community", side=2, line=2, outer=TRUE, at=0.5, cex=0.9)
mtext("Observed Bird Community", side=2, line=2, outer=TRUE, at=0.17, cex=0.9)

# Right side labels (y-axis)
#mtext("Taxonomic", side = 4, line = 3, outer = TRUE, at = 0.83, cex = 1.2)
#mtext("Functional", side = 4, line = 3, outer = TRUE, at = 0.5, cex = 1.2)
#mtext("Phylogenetic", side = 4, line = 3, outer = TRUE, at = 0.17, cex = 1.2)

# Bottom labels (x-axis)
mtext("Predicted Bird Community", side = 1, line = 3, outer = TRUE, at = 0.17, cex = 0.9)
mtext("Predicted Bird Community", side = 1, line = 3, outer = TRUE, at = 0.5, cex = 0.9)
mtext("Predicted Bird Community", side = 1, line = 3, outer = TRUE, at = 0.83, cex = 0.9)


dev.off()



######################################################################################################
##################### Print the t-values of the acoustic indices as table ############################

load("data/04_soundIndices_as_predictors/data_models_and_results.RData")


# Create a dataframe to store t-values for each axis
t_values_df <- do.call(rbind, lapply(t_values_results, function(x) as.data.frame(t(x))))

# Transpose: Sound indices as rownames, q0-2 as colnames
results <- as.data.frame(t(t_values_df))

# Add row names as a new column
results <- cbind("Acoustic Indices" = rownames(results), results)

# Remove row names 
rownames(results) <- NULL

# Delete the first row "Intercept"
results <- results[-1, ]

# Add spaces to the acoustic indices names
results$`Acoustic Indices` <- c("Soundscape Saturation", 
                "Entropy Of Variance Spectrum", 
                "Acoustic Complexity", 
                "Temporal Entropy", 
                "Events Per Second")

# Round t-values to 2 decimal places
results[, -1] <- round(as.numeric(as.matrix(results[, -1])), 2) # Exclude the first column


results_gt <- results %>%
  gt() %>%
  tab_spanner(
    label = "Taxonomic Diversity",
    columns = c(TD_q0_Axis1, TD_q1_Axis1, TD_q2_Axis1)
  ) %>%
  tab_spanner(
    label = "Functional Diversity",
    columns = c(FD_q0_Axis1, FD_q1_Axis1, FD_q2_Axis1)
  ) %>%
  tab_spanner(
    label = "Phylogenetic Diversity",
    columns = c(PD_q0_Axis1, PD_q1_Axis1, PD_q2_Axis1)
  ) %>%
  tab_spanner(
    label = "t-values",
    columns = c(
      TD_q0_Axis1, TD_q1_Axis1, TD_q2_Axis1,
      FD_q0_Axis1, FD_q1_Axis1, FD_q2_Axis1,
      PD_q0_Axis1, PD_q1_Axis1, PD_q2_Axis1
    )
  ) %>%
  cols_label(
    TD_q0_Axis1 = "q = 0", TD_q1_Axis1 = "q = 1", TD_q2_Axis1 = "q = 2",
    FD_q0_Axis1 = "q = 0", FD_q1_Axis1 = "q = 1", FD_q2_Axis1 = "q = 2",
    PD_q0_Axis1 = "q = 0", PD_q1_Axis1 = "q = 1", PD_q2_Axis1 = "q = 2"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  ) %>%
  tab_style(
    style = cell_fill(color = "#E5E5E5"), # Grey color for every second row
    locations = cells_body(rows = seq(2, nrow(results), by = 2)) # Apply to every second row
  )

gtsave(results_gt, "plots/t-values.png")


######################################################################################################
################################## Make a heatmap of the t-values ####################################

# Convert the 'results' data frame into a matrix for heatmap plotting
results_matrix <- as.matrix(results[, -1])  # Remove the 'Acoustic Indices' column
rownames(results_matrix) <- results$`Acoustic Indices`  # Set rownames to acoustic indices

# Reverse the order of the rows
results_matrix <- results_matrix[rev(rownames(results_matrix)), ]

# Save the heatmap as a PNG file
png("plots/heatmap_t_values.png", width = 800, height = 600)

# Custom heatmap with annotations
heatmap(
  results_matrix, 
  col = colorRampPalette(c("blue", "white", "red"))(80),  # Custom color palette
  scale = "none",  # Keep raw values (no scaling)
  margins = c(8, 12),  # Adjust margins for better readability
  labRow = rownames(results_matrix),  # Acoustic Indices as row labels
  labCol = colnames(results_matrix),  # Diversity types and Hill numbers as column labels
  #main = "t-values by Diversity Type and Hill Number",
  Colv = NA, Rowv = NA
)

# Close the graphics device (this actually saves the file)
dev.off()






