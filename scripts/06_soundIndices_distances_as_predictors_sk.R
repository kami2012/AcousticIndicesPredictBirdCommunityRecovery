rm(list=ls())

library(ecodist)

## ================================================================================================ ##
## Use biodiversity distance matrices based on expert data and sound indices distance matrices

## Check if sound indices distances predict biodiversity distances


# load biodiversity distances (taxonomic, functional, pyholgenetic / q0, q1, q2)
load("taxonomic_diversity/distances_com_exp_tax.rda")
load("functional_diversity/distances_com_exp_func.rda")
load("phylogenetic_diversity/distances_com_exp_phy.rda")


# load sound indices distances
load("data/distances_soundIndices.RData")

# Extract each sound indices distance matrix as a separate variable
distance_matrix_names <- names(distances_soundIndices)
for (name in distance_matrix_names) {
  assign(name, distances_soundIndices[[name]])
}


dis_com_exp_q0est_td <- dis_com_exp_q0est
dis_com_exp_q1est_td <- dis_com_exp_q1est
dis_com_exp_q2est_td <- dis_com_exp_q2est


# List of biodiversity matrices
biodiversity_matrices <- list(
  dis_com_exp_q0est_td, dis_com_exp_q1est_td, dis_com_exp_q2est_td,
  dis_com_exp_q0est_fd, dis_com_exp_q1est_fd, dis_com_exp_q2est_fd,
  dis_com_exp_q0est_pd, dis_com_exp_q1est_pd, dis_com_exp_q2est_pd
)

# Names for the matrices
matrix_names <- c(
  "TD_q0", "TD_q1", "TD_q2",
  "FD_q0", "FD_q1", "FD_q2",
  "PD_q0", "PD_q1", "PD_q2"
)

# Collect the values of the individual sound indices distances
# Initialize an empty dataframe to store results
soundIndDisAsPred <- data.frame(
  SoundIndices = character(),
  Diversity = character(),
  Estimates = numeric(),
  pValue = numeric(),
  SignificanceCode = character(),
  Rsquared = numeric(),
  Ftest = numeric(),
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
sink("data/soundIndicesDist_predict_biodiversityDist_summaries.txt")

# Iterate through all biodiversity matrices
for (i in 1:length(biodiversity_matrices)) {
  biodiv_matrix <- biodiversity_matrices[[i]]
  
  # Convert to distance object if needed
  biodiv_matrix <- as.dist(biodiv_matrix, upper = TRUE, diag = TRUE)
  
  # Create the formula for the model
  formula <- as.formula(paste("biodiv_matrix ~ dis_AcousticDiversity_s + dis_AcousticEveness_s + dis_BioAcoustic_s + dis_SoundscapeSaturation_s + dis_EntropyOfVarianceSpectrum_s + dis_AcousticComplexity_s + dis_TemporalEntropy_s + dis_EventsPerSecond_s"))
  
  # Fit the MRM model
  model <- MRM(formula, method = "linear")
  
  # Write the summary to the file
  cat("Summary for", matrix_names[i], "\n")
  print(model)
  cat("\n\n")
  
  # Extract coefficients and other statistics
  coef_summary <- model$coef
  
  # Ensure correct column names for estimates and p-values
  estimates <- coef_summary[, "biodiv_matrix"]
  pvalues <- coef_summary[, "pval"]
  
  # Convert p-values to significance codes
  significance_codes <- sapply(pvalues, get_significance_code)
  
  # Create temporary dataframe to hold current model results
  temp_df <- data.frame(
    SoundIndices = rownames(coef_summary),
    Diversity = matrix_names[i],
    Estimates = estimates,
    pValue = pvalues,
    SignificanceCode = significance_codes,
    Rsquared = model$r.squared["R2"], 
    Ftest = model$F.test["F"], 
    stringsAsFactors = FALSE
  )
  
  # Append the results to the main dataframe
  soundIndDisAsPred <- rbind(soundIndDisAsPred, temp_df)
}

# Close the summary file
sink()

# Save the soundIndDisAsPred dataframe as a CSV file
write.csv2(soundIndDisAsPred, file = "data/soundIndDisAsPred.csv", row.names = FALSE)

