rm(list=ls())

library(caret)
library(effects)

## ================================================================================================ ##
## Use abundance data / community composition "BirdCommunityNamePhylogeny_Indices.csv" from Joerg.
## There are 85 assemblages (rows) and 334 birds (columns), identified by experts.
## The table contains also 8 different sound indices and plot categories (Pasture/Cacao -> Old Growth)
## The table contains also the first and second nmds Axis for td, pd, and fd (q0, q1, q2)

## Test sound indices as predictors for the first axis of ordination of td, pd and fd  (q0, q1, q2)


## ------------ Define the data, response variable, and predictors --------------

# Load data
data_com <- read.csv("data/BirdCommunityNamePhylogeny_Indices_nmds.csv") 

# Define the response variables
nmds_axis1 <- c("TD_q0_est_Axis1", "TD_q1_est_Axis1", "TD_q2_est_Axis1",
                "FD_q0_est_Axis1", "FD_q1_est_Axis1", "FD_q2_est_Axis1",
                "PD_q0_est_Axis1", "PD_q1_est_Axis1", "PD_q2_est_Axis1")

# Define the predictors 
predictors <- c("SoundscapeSaturation", "EntropyOfVarianceSpectrum", "AcousticComplexity", 
                "TemporalEntropy", "EventsPerSecond")

#####################################################################################################################################
############################ Compute multiple R2 ####################################################################################

# Define parameters
n_splits <- 50
set.seed(123)

# Initialize lists to store results
r_squared_results <- list()

# Define the train-test split function
train_test_splits <- function(data, n_splits, seed) {
  set.seed(seed)
  splits <- replicate(n_splits, createDataPartition(data$Category10, p = 0.60, list = FALSE), simplify = FALSE)
  return(splits)
}

# Perform multiple splits and compute metrics
splits <- train_test_splits(data_com, n_splits, seed = 123)

for (nmdsAxis in nmds_axis1) {
  
  # Initialize vector to store R-squared values for this axis
  r_squared_values <- numeric(n_splits)
  
  for (i in seq_along(splits)) {
    
    # Split data
    train_indices <- splits[[i]]
    train_data <- data_com[train_indices, ]
    test_data <- data_com[-train_indices, ]
    
    # Create formula
    formula <- as.formula(paste(nmdsAxis, "~", paste(predictors, collapse = " + ")))
    
    # Fit model
    model <- lm(formula, data = train_data)
    
    # Predict and calculate R-squared
    predicted_values <- predict(model, newdata = test_data)
    r_squared_values[i] <- cor(predicted_values, test_data[[nmdsAxis]])^2
  }
  
  # Calculate summary statistics
  mean_r_squared <- mean(r_squared_values)
  sd_r_squared <- sd(r_squared_values)
  ci_low <- mean_r_squared - 1.96 * (sd_r_squared / sqrt(n_splits))
  ci_high <- mean_r_squared + 1.96 * (sd_r_squared / sqrt(n_splits))
  
  # Store results in the list
  r_squared_results[[nmdsAxis]] <- list(
    mean_r_squared = mean_r_squared,
    ci_low = ci_low,
    ci_high = ci_high
  )
}

# Save the results to a file
sink("data/soundIndices_predict_biodiversity_summaries_with_ci.txt")
cat("Summary for each NMDS axis with confidence intervals for R²:\n")
for (nmdsAxis in nmds_axis1) {
  # Retrieve stored results
  result <- r_squared_results[[nmdsAxis]]
  
  # Output results
  cat("Summary for ", nmdsAxis, ":\n", sep = "")
  cat("Mean R²:", round(result$mean_r_squared, 2), "\n")
  cat("95% CI for R²:", round(result$ci_low, 2), "-", round(result$ci_high, 2), "\n\n")
}
sink()


#########################################################################################################################
############################## Plot the predictions #####################################################################

# Order the Category10 levels to fit to the recovery gradient
desired_order <- c("A_Past", "F_Past", "A_Caca", "F_Caca", "F_PReg1", "F_CReg1", "F_PReg2", "F_CReg2", "A_Old", "F_Old")
data_com$Category10 <- factor(data_com$Category10, levels = desired_order)

# Define color mappings
colors_observed <- ifelse(startsWith(levels(data_com$Category10), "A"), rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5))
colors_predicted <- ifelse(startsWith(levels(data_com$Category10), "A"), rgb(1, 0, 1, 0.5), rgb(0, 0, 1, 0.5))


# Initialize lists to store models, predictions, and effects
models <- list()
predictions <- list()
partial_effects <- list()

# Collect the values of the individual sound indices
# Initialize an empty dataframe to store results
soundIndAsPred <- data.frame(
  SoundIndices = character(),
  Diversity = character(),
  Estimates = numeric(),
  StandardError = numeric(),
  tValue = numeric(),
  Significance = numeric(),
  stringsAsFactors = FALSE
)

# Function to convert p-value to significance code
get_significance_code <- function(p_value) {
  if (p_value <= 0.001) {
    return("***")
  } else if (p_value <= 0.01) {
    return("**")
  } else if (p_value <= 0.05) {
    return("*")
  } else if (p_value <= 0.1) {
    return(".")
  } else {
    return(" ")
  }
}

# Open a file to write the summaries
sink("data/soundIndices_predict_biodiversity_summaries.txt")

# Loop through each target column
for (nmdsAxis in nmds_axis1) {
  
  # Create a formula for the model
  formula <- as.formula(paste(nmdsAxis, "~", paste(predictors, collapse = " + ")))
  
  # Fit the linear regression model on the training data
  model <- lm(formula, data = train_data)
  
  # Store the model in the list
  models[[nmdsAxis]] <- model
  
  # Write the summary to the file
  cat("Summary for ", nmdsAxis, ":\n", sep = "")
  print(summary(model))
  cat("\n\n")
  
  # Predict the values using the test data
  predicted_values <- predict(model, newdata = test_data)
  
  # Store predicted values in the list
  predictions[[nmdsAxis]] <- predicted_values
  
  # Extract coefficients and other statistics
  coef_summary <- summary(model)$coefficients
  
  # Convert p-values to significance codes
  significance_codes <- sapply(coef_summary[, "Pr(>|t|)"], get_significance_code)
  
  # Create temporary dataframe to hold current model results
  temp_df <- data.frame(
    SoundIndices = rownames(coef_summary),
    Diversity = nmdsAxis,
    Estimates = coef_summary[, "Estimate"],
    StdError = coef_summary[, "Std. Error"],
    tValue = coef_summary[, "t value"],
    Significance = coef_summary[, "Pr(>|t|)"],
    SignificanceCode = significance_codes,
    stringsAsFactors = FALSE
  )
  
  # Append the results to the main dataframe
  soundIndAsPred <- rbind(soundIndAsPred, temp_df)
  
  
  #------------------------------------------------------------------------------------------------------------------
  
  # Plot sound indices as predictors for nmds Axis 1 (predicted vs. observed)
  png(filename = paste0("plots/SoundIndicesPredictComm/", nmdsAxis, ".png"), width = 800, height = 600)
  
  plot(predicted_values, test_data[[nmdsAxis]], xlab = "Predicted values by sound indices", ylab = nmdsAxis)
  abline(lm(test_data[[nmdsAxis]] ~ predicted_values), col = "blue")
  
  # Calculate and display R-squared value
  # a) indicates how well the model fits to the training data
  #r_squared <- summary(model)$r.squared
  
  # b) indicates how well the predictions match the test data
  r_squared <- cor(predicted_values, test_data[[nmdsAxis]])^2
  text(min(predicted_values), max(test_data[[nmdsAxis]]), paste("R² =", round(r_squared, 2)), pos = 4, col = "red")
  
  dev.off()
  
  #------------------------------------------------------------
  
  # Boxplot nmds Axis 1 and predicted values against recovery gradient (observed vs. predicted values by recovery gradient)
  png(filename = paste0("plots/RecoveryGradient/", nmdsAxis, ".png"), width = 800, height = 600)
  
  # Observed Values
  boxplot(test_data[[nmdsAxis]] ~ test_data$Category10, xlab="Recovery gradient", 
          ylab=nmdsAxis, col=rgb(1, 0, 0, 0.5))
  
  # Predicted Values
  par(new=TRUE)
  boxplot(predicted_values ~ test_data$Category10, xlab="", ylab="", 
          col=rgb(0, 0, 1, 0.5), axes=FALSE, at=1:length(unique(data_com$Category10)) + 0.2)
  
  # Add a legend
  legend("bottomright", legend=c("Observed values", "Predicted values"), fill=c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), bty="o")
  
  dev.off()
  
}

# Close the summary file
sink()

# Save the models, predictions and partial effects lists
save(models, predictions, partial_effects, file = "data/models_predictions_partialEffects.RData")

# Identify numeric columns
numeric_cols <- sapply(soundIndAsPred, is.numeric)

# Round all numeric columns to 2 decimal places
soundIndAsPred[numeric_cols] <- lapply(soundIndAsPred[numeric_cols], round, digits = 2)

# Save the soundIndAsPred dataframe as a CSV file
write.csv2(soundIndAsPred, file = "data/soundIndAsPred.csv", row.names = FALSE)

##########################################################################################################
################# Comprehensive plots ####################################################################
##########################################################################################################
# Comprehensive plot: Sound indices as predictors for nmds Axis 1

png(filename="plots/soundIndicesAsPred.png", width = 170, height = 200, units = "mm", res=1000)

layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol=3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1, 1))
par(mar = c(4, 4, 2, 2) + 0.1, oma = c(4, 4, 4, 1))

# Loop through each NMDS axis for plotting
for (nmdsAxis in nmds_axis1) {
  # Select the predictions
  predicted_values <- predictions[[nmdsAxis]]
  observed_values <- test_data[[nmdsAxis]]
  
  plot(predicted_values, observed_values, 
       xlab = "Acoustic Indices", 
       ylab = "NMDS Axis1",
       pch = 1,
       col = rgb(0, 0, 0, 0.5))  # Black points with 50% transparency
  
  # Fit a linear model
  model <- lm(observed_values ~ predicted_values)
  
  # Get the range of predicted values (from the first to last point)
  x_range <- range(predicted_values)
  
  # Predict y-values for the first and last x-values
  y_range <- predict(model, newdata = data.frame(predicted_values = x_range))
  
  # Draw a line between the first and last points
  lines(x_range, y_range, col = "blue", lwd = 2)
  
  # Calculate and display R-squared value
  r_squared <- cor(predicted_values, observed_values)^2
 
  # Retrieve stored results
  mean_r_squared <- r_squared_results[[nmdsAxis]]$mean_r_squared
  ci_low <- r_squared_results[[nmdsAxis]]$ci_low
  ci_high <- r_squared_results[[nmdsAxis]]$ci_high
  
  # Compute confidence intervals for the y-values
  pred_ci <- predict(model, newdata = data.frame(predicted_values = x_range), interval = "confidence")
  
  # Create polygon for confidence interval
  polygon(c(x_range, rev(x_range)),
          c(pred_ci[, "lwr"], rev(pred_ci[, "upr"])),
          col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Display R²
  text_x <- min(predicted_values)  # X coordinate for the text
  text_y <- max(observed_values) - 0.05 * diff(range(observed_values)) # Y coordinate for the R² text
  
  text(text_x, text_y, 
       paste("R² =", round(mean_r_squared, 2)), 
       pos=4, col="blue")
  
}

# Add labels for rows and columns
mtext("q0", side=3, line=0.5, outer=TRUE, at=0.167, cex=1.2)
mtext("q1", side=3, line=0.5, outer=TRUE, at=0.5, cex=1.2)
mtext("q2", side=3, line=0.5, outer=TRUE, at=0.833, cex=1.2)

mtext("Taxonomic Diversity", side=2, line=2, outer=TRUE, at=0.83, cex=1.2)
mtext("Functional Diversity", side=2, line=2, outer=TRUE, at=0.5, cex=1.2)
mtext("Phylogenetic Diversity", side=2, line=2, outer=TRUE, at=0.17, cex=1.2)

dev.off()

##########################################################################################################
# Comprehensive plot: Boxplot against Recovery gradient

png(filename="plots/recoveryGradient.png", width = 170, height = 200, units = "mm", res=1000)

layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol=3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1, 1))
par(mar = c(2, 2, 2, 1), oma = c(8, 4, 4, 1))

# Loop through each NMDS axis for plotting
for (nmdsAxis in nmds_axis1) {
  # Select the predictions
  predicted_values <- predictions[[nmdsAxis]]
  
  # Observed Values
  boxplot(test_data[[nmdsAxis]] ~ test_data$Category10, col=colors_observed)
  
  # Predicted Values
  par(new=TRUE)
  boxplot(predicted_values ~ test_data$Category10, 
          col=colors_predicted, axes=FALSE, at=1:length(unique(test_data$Category10)) + 0.2)
  
  # Add a legend
  legend("bottomright", legend=c("Observed (A)", "Predicted (A)", "Observed (F)", "Predicted (F)"), 
         fill=c(rgb(1, 0, 0, 0.5), rgb(1, 0, 1, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), bty="o", cex=0.7)
  
}

# Add labels for rows and columns
mtext("q0", side=3, line=0.5, outer=TRUE, at=0.167, cex=1.2)
mtext("q1", side=3, line=0.5, outer=TRUE, at=0.5, cex=1.2)
mtext("q2", side=3, line=0.5, outer=TRUE, at=0.833, cex=1.2)

mtext("Taxonomic Diversity", side=2, line=2, outer=TRUE, at=0.83, cex=1.2)
mtext("Functional Diversity", side=2, line=2, outer=TRUE, at=0.5, cex=1.2)
mtext("Phylogenetic Diversity", side=2, line=2, outer=TRUE, at=0.17, cex=1.2)

# Add information about the recovery gradient below the plots
mtext(text=paste("Recovery gradients:\n", paste(desired_order, collapse=" - ")), side=1, line=6, outer=TRUE, cex=0.8, las=1, adj=0.05)

dev.off()


##########################################################################################################
# Comprehensive dataframe: Individual sound indices as predictors

