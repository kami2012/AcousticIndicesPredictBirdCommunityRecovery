# Guidelines for Practitioners

# Table of Contents

1. [Collect Sound Recordings](#step-1-collect-sound-recordings)  
2. [Compute Acoustic Indices](#step-2-compute-acoustic-indices)  
3. [Make a Validation Dataset](#step-3-make-a-validation-dataset)  
4. [Compute Pairwise β-Diversity Indices](#step-4-compute-pairwise-β-diversity-indices)  
5. [Perform NMDS Ordination](#step-5-perform-nmds-ordination)  
6. [Train Linear Models](#step-6-train-linear-models)  
7. [Analyze All Data](#step-7-analyze-all-data)  


## Step 1: Collect Sound Recordings

**Description:** Record audio in all plots of interest.  

**Input:**  

**Analysis:** 

**Output:** [sound files]  

**Comments:**  

## Step 2: Compute Acoustic Indices

**Description:** Compute the following 5 acoustic indices for all [sound files]:  
- Soundscape Saturation
- Entropy Of Variance Spectrum
- Acoustic Complexity
- Temporal Entropy
- Events Per Second

**Input:** [sound files]  

**Analysis:** AnalysisProgram.exe 

**Output:**   

**Comments:**  

## Step 3: Make a Validation Dataset

**Description:** Select a representative subset of the [sound files] for validation. These files should be manually labeled by experts to ensure high-quality reference data.   

**Input:** [sound files]  

**Analysis:** Experts label species presence (1) /absence (0) in each file.  

**Output:** data/detections_freile_gelis_2021_2022_birds_dummy_pivot_reordered.csv [community dataset]  

![image](https://github.com/user-attachments/assets/16ebc402-7c77-482d-ae2f-aaa72cec6f84)

**Comments:**  The first column must be named "filename", the second "plot".

## Step 4: Compute Pairwise β-Diversity Indices

**Description:** Compute pairwise β-diversity indices (dissimilarity matrices), taking into account three dimensions of biodiversity
- “taxonomic diversity” (TD)
- “phylogenetic diversity” (PD)
- “functional diversity” (FD)

and focusing on 
- infrequent (q = 0)
- frequent (q = 1)
- highly frequent (q = 2)

species (based on Jaccard, Horn, and Morisita-Horn index).  

**Script:** scripts/02_betadiv_matrices_sk.R

**Input:** 

| Biodiversity Type         | Input Dataset                         | Additional Data                     |
|--------------------------|--------------------------------------|-------------------------------------|
| **Taxonomic Diversity (TD)** | `taxonomic_diversity/com_td.csv`  [community dataset]   | -                                   |
| **Phylogenetic Diversity (PD)** | `phylogenetic_diversity/com_pd.csv` [community dataset] | `bird_tree.rda` (Newick format)     |
| **Functional Diversity (FD)** | `functional_diversity/com_fd.csv` [community dataset]  | `trait_matrix.rda` (Gower distance) |


**Analysis:** iNEXTbeta3D_pair3D (an adaptation of iNEXT.beta3D [*iNEXT.beta3D* package])

```
pairwise_TD91 = iNEXTbeta3D_pair3D(com_td, div0="TD", SC = SC, datatype0 = "incidence_raw", parallel = T, cpus = 8)

pairwise_PD91 = iNEXTbeta3D_pair3D(com_pd, div0="PD", SC = SC, datatype0 = "incidence_raw", PDTree0 = tr, 
                                   parallel = T, cpus = 8)

pairwise_FD91 = iNEXTbeta3D_pair3D(com_fd, div0="FD", SC = SC, datatype0 = "incidence_raw", FDdistM0 = distM, 
                                  parallel = T, cpus = 8)
```

The resulting similarity matrices must be converted to dissimilarity matrices, e.g.:

```
dis_com_q0est_pd <- 1- pairwise_PD91[["Matrices"]][["jac_est"]]
dis_com_q1est_pd <- 1- pairwise_PD91[["Matrices"]][["hor_est"]]
dis_com_q2est_pd <- 1- pairwise_PD91[["Matrices"]][["mor_hor_est"]]
```

**Output:** 
- taxonomic_diversity/distances_com_exp_tax.rda [dissimilarity matrices]
- phylogenetic_diversity/distances_com_exp_phy.rda [dissimilarity matrices]
- functional_diversity/distances_com_exp_func.rda [dissimilarity matrices]

**Comments:**
- The sampling coverage (SC) can be calculated with DataInfobeta3D [*iNEXT.beta3D* package]
- [tree] must be in Newick format
- [traits] must be a pairwise distance matrix (Gower distance)


## Step 5: Perform NMDS Ordination

**Description:** To reduce the multidimensional complexity of the dissimilarity matrices to a two-dimensional representation (axis1 and axis2), perform an ordination with each distance matrix (TD, PD and FD) and for all orders of q (q = 0, q = 1 and q = 2). 

**Script:** scripts/03_nmds_sk.R  

**Input:** 
- taxonomic_diversity/distances_com_exp_tax.rda [dissimilarity matrices]
- phylogenetic_diversity/distances_com_exp_phy.rda [dissimilarity matrices]
- functional_diversity/distances_com_exp_func.rda [dissimilarity matrices]

**Analysis:** metaMDS [*vegan* package]

```
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

# Compute NMDS for each matrix
for (name in names(dissimilarity_matrices)) {
  set.seed(1)
  nmds <- metaMDS(dissimilarity_matrices[[name]], k = 2, maxit = 999, trymax = 100)
  nmds_results[[name]] <- nmds
}
```

**Output:** data/nmds_results.rds  

**Comments:** 
- Set a random seed for reproducibility.
- Check stress values (<0.2 is good, <0.1 is excellent).
- Adjust k (number of dimensions) if needed.


## Step 6: Train Linear Models

**Description:** Take the resulting nmds axis1 values as response variables in linear models, with the five acoustic indices (Soundscape Saturation, Entropy Of Variance Spectrum, Acoustic Complexity, Temporal Entropy and Events Per Second) as predictor variables. 

**Script:** scripts/04_soundIndices_as_predictors_sk.R  

**Input:** data/plots_categories_indices_nmds.csv [Dataset that contains the plot IDs, the plot categories, the acoustic indices and the nmds axis1 values]

![image](https://github.com/user-attachments/assets/74d9a5bb-9de4-466d-b294-a9141d61414e)

Split the data into train (0.66) and test (0.33) by picking every 3rd row:

```
# Load data
data_com <- read.csv("data/plots_categories_indices_nmds.csv")

# Sort by nmds axis1
data_com <- data_com[order(data_com$TD_q0_Axis1), ]

# Every third row for the test set
test_data <- data_com[seq(3, nrow(data_com), by = 3), ]

# Remaining rows (those not in the test set) for the training set
train_data <- data_com[-seq(3, nrow(data_com), by = 3), ]
```

**Analysis:**

```
# Define the response variables
nmds_axis1 <- c("TD_q0_Axis1", "TD_q1_Axis1", "TD_q2_Axis1",
                "FD_q0_Axis1", "FD_q1_Axis1", "FD_q2_Axis1",
                "PD_q0_Axis1", "PD_q1_Axis1", "PD_q2_Axis1")

# Define the predictors 
predictors <- c("SoundscapeSaturation", "EntropyOfVarianceSpectrum", "AcousticComplexity", 
                "TemporalEntropy", "EventsPerSecond")


# Loop over each response variable 
for (axis in nmds_axis1) {
  
  # Create formula for the current response variable
  formula <- as.formula(paste(axis, "~", paste(predictors, collapse = " + ")))
  
  # Fit linear model to the training data
  model <- lm(formula, data = train_data)

  # Predict on the test data
  predictions_tdata <- predict(model, newdata = test_data)

}

```

**Output:** data/data_models_and_results.RData  

**Comments:** Evaluate the models`performance by checking R^2 or the t-values.  

## Step 7: Analyze All Data

**Description:** If your model(s) can predict your test_data well, you can use them for analyzing your remaining sound data. 

**Input:** [sound files] 

**Analysis:** 

```
# Predict on the sound_files
predictions_sound_files <- predict(model, newdata = sound_files)

```

**Output:** 
**Comments:**

