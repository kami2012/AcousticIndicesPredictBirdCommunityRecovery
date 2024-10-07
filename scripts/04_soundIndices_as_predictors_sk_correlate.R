rm(list=ls())

library(effects)

## ================================================================================================ ##
## Use abundance data / community composition "BirdCommunityNamePhylogeny_Indices.csv" from Joerg.
## There are 85 assemblages (rows) and 334 birds (columns), identified by experts.
## The table contains also 8 different sound indices and plot categories (Pasture/Cacao -> Old Growth)
## The table contains also the first and second nmds Axis for td, pd, and fd (q0, q1, q2)

## Test sound indices as predictors for the first axis of ordination of td, pd and fd  (q0, q1, q2)


## ------------ Define the data, response variable, and predictors --------------

# 1. DATA
data_com <- read.csv("data/BirdCommunityNamePhylogeny_Indices_nmds.csv") 

## Split data_com in a training and a testing dataset
# Set seed for reproducibility
set.seed(123)

# Split the data into training (2/3) and testing (1/3) sets
sample_size <- floor(0.67 * nrow(data_com))
train_indices <- sample(seq_len(nrow(data_com)), size = sample_size)

# Create training and testing datasets
train_data <- data_com[train_indices, ]
test_data <- data_com[-train_indices, ]


# 2. RESPONSE VARIABLES
# Select the nmds axis1 (q0est_td_Axis1, q0est_td_Axis2, q1est_td_Axis1, ...) as response variables
nmds_axis1 <- c("TD_q0_est_Axis1", "TD_q1_est_Axis1", "TD_q2_est_Axis1",
                    "FD_q0_est_Axis1", "FD_q1_est_Axis1", "FD_q2_est_Axis1",
                    "PD_q0_est_Axis1", "PD_q1_est_Axis1", "PD_q2_est_Axis1")



# Order the Category10 levels to fit to the recovery gradient
desired_order <- c("A_Past", "F_Past", "A_Caca", "F_Caca", "F_PReg1", "F_CReg1", "F_PReg2", "F_CReg2", "A_Old", "F_Old")
data_com$Category10 <- factor(data_com$Category10, levels = desired_order)



# Define color mappings
colors_observed <- ifelse(startsWith(levels(data_com$Category10), "A"), rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5))
colors_predicted <- ifelse(startsWith(levels(data_com$Category10), "A"), rgb(1, 0, 1, 0.5), rgb(0, 0, 1, 0.5))


# Store models and partial effect plots in lists
models <- list()
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
  formula <- as.formula(paste(nmdsAxis, "~ AcousticDiversity + AcousticEveness + BioAcoustic + SoundscapeSaturation + EntropyOfVarianceSpectrum + AcousticComplexity + TemporalEntropy + EventsPerSecond"))
  
  # Fit the linear regression model
  model <- lm(formula, data = data_com)
  
  # Store the model in the list
  models[[nmdsAxis]] <- model
  
  # Write the summary to the file
  cat("Summary for ", nmdsAxis, ":\n", sep = "")
  print(summary(model))
  cat("\n\n")
  
  # Predict the values
  predicted_values <- predict(model)
  
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
  
  # Plot sound indices as predictors for nmds Axis 1
  png(filename = paste0("plots/SoundIndicesPredictComm/", nmdsAxis, ".png"), width = 800, height = 600)
  plot(predicted_values, data_com[[nmdsAxis]], xlab = "Predicted values by sound indices", ylab = nmdsAxis)
  abline(lm(data_com[[nmdsAxis]] ~ predicted_values), col = "blue")
  
  # Calculate and display R-squared value
  r_squared <- summary(model)$r.squared
  text(min(predicted_values), max(data_com[[nmdsAxis]]), paste("R² =", round(r_squared, 2)), pos = 4, col = "red")
  dev.off()
  
  #------------------------------------------------------------
  
  # Boxplot nmds Axis 1 and predicted values against recovery gradient
  png(filename = paste0("plots/RecoveryGradient/", nmdsAxis, ".png"), width = 800, height = 600)
  
  # Observed Values
  boxplot(data_com[[nmdsAxis]] ~ data_com$Category10, xlab="Recovery gradient", 
          ylab=nmdsAxis, col=rgb(1, 0, 0, 0.5))
  
  # Predicted Values
  par(new=TRUE)
  boxplot(predicted_values ~ data_com$Category10, xlab="", ylab="", 
          col=rgb(0, 0, 1, 0.5), axes=FALSE, at=1:length(unique(data_com$Category10)) + 0.2)
  
  # Add a legend
  legend("bottomright", legend=c("Observed values", "Predicted values"), fill=c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), bty="o")
  
  dev.off()
  
  #------------------------------------------------------------
  
  # Calculate partial effects without residuals
  partial_effect <- allEffects(model, residuals = FALSE)
  
  # Store the partial effect plot in the list
  partial_effects[[nmdsAxis]] <- partial_effect
  
  # Plot the partial effect and save as PNG
  png(filename = paste0("plots/PartialEffects/", nmdsAxis, ".png"), width = 800, height = 600)
  plot(partial_effect, main = paste("Effect plot for", nmdsAxis))
  dev.off()
  
}

# Close the summary file
sink()

# Save the models and partial effects list 
save(models, partial_effects, file = "data/models_and_partialEffects.RData")

# Identify numeric columns
numeric_cols <- sapply(soundIndAsPred, is.numeric)

# Round all numeric columns to 2 decimal places
soundIndAsPred[numeric_cols] <- lapply(soundIndAsPred[numeric_cols], round, digits = 2)

# Save the soundIndAsPred dataframe as a CSV file
write.csv2(soundIndAsPred, file = "data/soundIndAsPred.csv", row.names = FALSE)

##########################################################################################################
##########################################################################################################
# Comprehensive plot: Recovery gradient

png(filename="plots/recoveryGradient.png", width = 170, height = 200, units = "mm", res=1000)

layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol=3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1, 1))
par(mar = c(2, 2, 2, 1), oma = c(8, 4, 4, 1))

# Loop through each NMDS axis for plotting
for (nmdsAxis in nmds_axis1) {
  # Predict the values
  predicted_values <- predict(models[[nmdsAxis]])
  
  # Observed Values
  boxplot(data_com[[nmdsAxis]] ~ data_com$Category10, col=colors_observed)
  
  # Predicted Values
  par(new=TRUE)
  boxplot(predicted_values ~ data_com$Category10, 
          col=colors_predicted, axes=FALSE, at=1:length(unique(data_com$Category10)) + 0.2)
  
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
# Comprehensive plot: Sound indices as predictors for nmds Axis 1

png(filename="plots/soundIndicesAsPred.png", width = 170, height = 200, units = "mm", res=1000)

layout(matrix(c(1, 2, 3, 
                4, 5, 6, 
                7, 8, 9), ncol=3, byrow=TRUE), widths=c(1, 1, 1), heights=c(1, 1, 1))
par(mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 1))

# Loop through each NMDS axis for plotting
for (nmdsAxis in nmds_axis1) {
  # Predict the values
  predicted_values <- predict(models[[nmdsAxis]])
  
  plot(predicted_values, data_com[[nmdsAxis]], xlab = "Predicted", ylab = "Observed (NMDS Axis1)")
  abline(lm(data_com[[nmdsAxis]] ~ predicted_values), col = "blue")
  
  # Calculate and display R-squared value
  r_squared <- summary(models[[nmdsAxis]])$r.squared
  text(min(predicted_values), max(data_com[[nmdsAxis]]), paste("R² =", round(r_squared, 2)), pos = 4, col = "red")

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
# Comprehensive dataframe: Individual sound indices as predictors

