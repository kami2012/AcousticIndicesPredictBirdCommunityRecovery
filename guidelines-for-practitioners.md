# Guidelines for Practitioners

## Step 1: Collect sound recordings

**Description:** Record audio in all plots of interest.  

**Input:**  

**Analysis:** 

**Output:** [sound files]  

**Comments:**  

## Step 2: Compute acoustic indices

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

## Step 3: Make a Validation Subset

**Description:** Select a representative subset of the [sound files] for validation. These files should be manually labeled by experts to ensure high-quality reference data.   

**Input:** [sound files]  

**Analysis:** Experts label species presence (1) /absence (0) in each file.  

**Output:** data/detections_freile_gelis_2021_2022_birds_dummy_pivot_reordered.csv [community dataset]  
![image](https://github.com/user-attachments/assets/16ebc402-7c77-482d-ae2f-aaa72cec6f84)

**Comments:**  

## Step 4: Compute pairwise β-diversity indices (distance matrices) 

**Description:** Compute pairwise β-diversity indices (distance matrices), taking into account three dimensions of biodiversity
- “taxonomic diversity” (TD)
- “phylogenetic diversity” (PD)
- “functional diversity” (FD)

and focusing on 
- infrequent (q = 0)
- frequent (q = 1)
- highly frequent (q = 2)

species (based on Jaccard, Horn, and Morisita-Horn index).  


**Input:** 

*TD*
-  taxonomic_diversity/com_td.csv [community dataset]

*PD*
- phylogenetic_diversity/com_pd.csv [community dataset]
- phylogenetic_diversity/bird_tree.rda [tree]

*FD*
- functional_diversity/com_fd.csv [community dataset]
- functional_diversity/trait_matrix.rda [traits]


**Analysis:** iNEXTbeta3D_pair3D (an adaptation of iNEXT.beta3D)

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

**Output** 
- taxonomic_diversity/distances_com_exp_tax.rda [dissimilarity matrices]
- phylogenetic_diversity/distances_com_exp_phy.rda [dissimilarity matrices]
- functional_diversity/distances_com_exp_func.rda [dissimilarity matrices]

**Comments**
- The sampling coverage (SC) can be calculated with DataInfobeta3D
- [tree] must be in Newick format
- [traits] must be a pairwise distance matrix (Gower distance)


## Step 5: Perform ordination

## Step 6: Divide the validation subset into train and test

## Step 7: Train linear models

## Step 8: Analyze all data [sound files] 





