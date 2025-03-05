# Guidelines for Practitioners

## Step 1: Collect sound recordings

*Description:* Record audio in all plots of interest.  
*Input:*  
*Analysis:*  
*Output:* [sound files]  
*Comments:*  

## Step 2: Compute acoustic indices

*Description:* Compute the following 5 acoustic indices for all [sound files]:  
- Soundscape Saturation
- Entropy Of Variance Spectrum
- Acoustic Complexity
- Temporal Entropy
- Events Per Second

*Input:* [sound files]  
*Analysis:* AnalysisProgram.exe  
*Output:*   
*Comments:*  

## Step 3: Make a Validation Subset

*Description:* Select a representative subset of the [sound files] for validation. These files should be manually labeled by experts to ensure high-quality reference data.  
*Input:* [sound files]  
*Analysis:* Experts label species presence (1) /absence (0) in each file.  
*Output:* data/detections_freile_gelis_2021_2022_birds_dummy_pivot_reordered.csv [community dataset]  
![image](https://github.com/user-attachments/assets/16ebc402-7c77-482d-ae2f-aaa72cec6f84)

*Comments:*  

## Step 4: Compute pairwise β-diversity indices (distance matrices) 

*Description:* Compute pairwise β-diversity indices (distance matrices), taking into account three dimensions of biodiversity – “taxonomic diversity” (TD), “phylogenetic diversity” (PD) and “functional diversity” (FD) – and focusing on infrequent (q = 0), frequent (q = 1), and highly frequent (q = 2) species (based on Jaccard, Horn, and Morisita-Horn index).  
*Input:* 
-  [community dataset]
-  






