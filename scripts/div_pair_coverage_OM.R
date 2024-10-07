# Pairwise similarity indices by Anne Chao
# to be loaded in analysis script by "load()"
# to provide iNEXTbeta3D_pair functions
# updated by Oliver, 15.04.2024

#install.packages("iNEXT.3D")      ## Please update iNEXT.3D to version 1.0.2 from CRAN
#install.packages("iNEXT.beta3D")  ## version 1.0.0 from CRAN
library(iNEXT.beta3D)

library(snowfall)
library(dplyr)

#' Compute observed and coverage-based standardized (estimated) Jaccard similarity index (q=0), Sorensen similarity index (q=0),
#' Horn similarity index (q=1), and Morisita-Horn similarity index (q=2) for all pairs of assemblages.
#' NOTE: The observed similarity indices typically underestimates; they are computed only for comparison.
#' We suggest researchers use the coverage-based standardized similarity indices.
#' 
#' @param com a data.frame/matrix with assemblages (rows) by species (columns).
#' @param SC a standardized sample coverage value for computing four standardized (estimated) similarity indices. 
#' By default (\code{SC = NULL}), based on all pairs of samples, SC represents the minimum among the coverage values for alpha reference samples
#' extrapolated to double the total abundance/size.  
#' The default coverage value is referred to as C_2n,alpha in Chao et al. (2023, Ecological Monographs) or Cmax_joint in this code; 
#' it represents the maximum coverage value for which one can reliably infer (dia)similarity and beta diversity. This value varies with dataset. 
#' Users can specify a common value of SC for the argument (e.g., (\code{SC=0.9})) when there are many datasets under comparison.

#' @param parallel whether to do parallel computation or not.
#' @param cpus number of cpu cores for parallel.
#' @return a data frame including observed and standardized Jaccard similarity, observed and standardized Sorensen similarity, 
#' observed and standardized Horn index, as well as observed and standardized Morisita-Horn index for all pairs of assemblages
#' under the default or user-specified sample coverage (SC).

iNEXTbeta3D_pair <- function(com, SC = NULL, parallel = F, cpus = 2){  
  
  m <- matrix(NA, ncol = dim(com)[1], nrow = dim(com)[1], byrow = T)
  
  pair <- data.frame(pos = 1:dim(com)[1]^2,
                     x = rep(rownames(com), nrow(m)), 
                     y = as.vector(sapply(rownames(com), rep, nrow(m))) )
  
  if (is.null(SC)) SC = sapply(pair$pos, function(i) DataInfobeta3D( t(com[c(pair$x[i], pair$y[i]),]) )$`SC(2n)`[4]) %>% min
  
  
  # for parallelization
  if (parallel == TRUE){
    
    loop <- function(v, com = com){
      
      sfCat(paste(v))
      
      out <- iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), q = c(0, 1, 2), level = SC, nboot = 0)[[1]]
      
      out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), q = c(0, 1, 2), base = "size", 
                            level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0)[[1]]
      beta.obs = out.obs$gamma$Gamma / out.obs$alpha$Alpha
      
      
      ## Observed Jaccard similarity q = 0 (U02) index 
      jac_obs = 1 - (1/beta.obs[1] - 1) / (1/2 - 1)
      
      ## Standardized Jaccard similarity q = 0 (U02) index
      jac_est = 1 - out$`1-U`$Dissimilarity[1]
      
      ## Observed Sorenson similarity q = 0 (C02) index
      sor_obs = 1 - (beta.obs[1] - 1) / (2 - 1)
      
      ## Standardized Sorenson similarity q = 0 (C02) index
      sor_est = 1 - out$`1-C`$Dissimilarity[1]
      
      ## Observed Horn similarity q = 1 (C12 = U12) index
      hor_obs = 1 - log(beta.obs[2]) / log(2)
      
      ## Standardized Horn similarity q = 1 (C12 = U12) index
      hor_est = 1 - out$`1-C`$Dissimilarity[2]
      
      ## Observed Morisita Horn similarity q = 2 (C22) index
      mor_hor_obs = 1 - (1/beta.obs[3] - 1) / (1/2 - 1)
      
      ## Standardized Morisita Horn similarity q = 2 (C22) index
      mor_hor_est = 1 - out$`1-C`$Dissimilarity[3]
      
      
      return(c(jac_obs = jac_obs,         jac_est = jac_est, 
               sor_obs = sor_obs,         sor_est = sor_est,
               hor_obs = hor_obs,         hor_est = hor_est, 
               mor_hor_obs = mor_hor_obs, mor_hor_est = mor_hor_est))
    }
    
    sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
    sfExport("pair"); sfExport("com")
    sfExport("loop"); sfLibrary(iNEXT.beta3D); sfLibrary(snowfall)
    
    res <- sfLapply(pair$pos, fun = loop, com = com)
    sfStop()
    
  } else {
    
    loop <- function(v, com = com){
      
      print(v)
      
      out <- iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), q = c(0, 1, 2), level = SC, nboot = 0)[[1]]
      
      out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), q = c(0, 1, 2), base = "size", 
                            level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0)[[1]]
      beta.obs = out.obs$gamma$Gamma / out.obs$alpha$Alpha
      
      
      ## Observed Jaccard similarity q = 0 (U02) index 
      jac_obs = 1 - (1/beta.obs[1] - 1) / (1/2 - 1)
      
      ## Standardized Jaccard similarity q = 0 (U02) index
      jac_est = 1 - out$`1-U`$Dissimilarity[1]
      
      ## Observed Sorenson similarity q = 0 (C02) index
      sor_obs = 1 - (beta.obs[1] - 1) / (2 - 1)
      
      ## Standardized Sorenson similarity q = 0 (C02) index
      sor_est = 1 - out$`1-C`$Dissimilarity[1]
      
      ## Observed Horn similarity q = 1 (C12 = U12) index
      hor_obs = 1 - log(beta.obs[2]) / log(2)
      
      ## Standardized Horn similarity q = 1 (C12 = U12) index
      hor_est = 1 - out$`1-C`$Dissimilarity[2]
      
      ## Observed Morisita Horn similarity q = 2 (C22) index
      mor_hor_obs = 1 - (1/beta.obs[3] - 1) / (1/2 - 1)
      
      ## Standardized Morisita Horn similarity q = 2 (C22) index
      mor_hor_est = 1 - out$`1-C`$Dissimilarity[3]
      
      
      return(c(jac_obs = jac_obs,         jac_est = jac_est, 
               sor_obs = sor_obs,         sor_est = sor_est,
               hor_obs = hor_obs,         hor_est = hor_est, 
               mor_hor_obs = mor_hor_obs, mor_hor_est = mor_hor_est))
    }
    
    res <- lapply(pair$pos, loop, com = com)
  }
  
  # finalization
  res <- do.call(rbind, res)
  res[res>1] <- 1  
  
  return(res)
  
}


###########
iNEXTbeta3D_pair2 <- function(com, SC = NULL, parallel = F, cpus = 2, div0="TD", PDTree0=NULL, FDdistM0=NULL){  
  
  m <- matrix(NA, ncol = dim(com)[1], nrow = dim(com)[1], byrow = T)
  
  pair <- data.frame(pos = 1:dim(com)[1]^2,
                     x = rep(rownames(com), nrow(m)), 
                     y = as.vector(sapply(rownames(com), rep, nrow(m))) )
  
  if (is.null(SC)) SC = sapply(pair$pos, function(i) DataInfobeta3D( t(com[c(pair$x[i], pair$y[i]),]) )$`SC(2n)`[4]) %>% min
  
  
  # for parallelization
  if (parallel == TRUE){
    
    loop <- function(v, com = com){
      
      sfCat(paste(v))
      
      out <- iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), level = SC, nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0)[[1]]
      
      out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), base = "size", 
                            level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0)[[1]]
      beta.obs = out.obs$gamma$Gamma / out.obs$alpha$Alpha
      
      
      ## Observed Jaccard similarity q = 0 (U02) index 
      jac_obs = 1 - (1/beta.obs[1] - 1) / (1/2 - 1)
      
      ## Standardized Jaccard similarity q = 0 (U02) index
      jac_est = 1 - out$`1-U`$Dissimilarity[1]
      
      ## Observed Sorenson similarity q = 0 (C02) index
      sor_obs = 1 - (beta.obs[1] - 1) / (2 - 1)
      
      ## Standardized Sorenson similarity q = 0 (C02) index
      sor_est = 1 - out$`1-C`$Dissimilarity[1]
      
      ## Observed Horn similarity q = 1 (C12 = U12) index
      hor_obs = 1 - log(beta.obs[2]) / log(2)
      
      ## Standardized Horn similarity q = 1 (C12 = U12) index
      hor_est = 1 - out$`1-C`$Dissimilarity[2]
      
      ## Observed Morisita Horn similarity q = 2 (C22) index
      mor_hor_obs = 1 - (1/beta.obs[3] - 1) / (1/2 - 1)
      
      ## Standardized Morisita Horn similarity q = 2 (C22) index
      mor_hor_est = 1 - out$`1-C`$Dissimilarity[3]
      
      
      return(c(jac_obs = jac_obs,         jac_est = jac_est, 
               sor_obs = sor_obs,         sor_est = sor_est,
               hor_obs = hor_obs,         hor_est = hor_est, 
               mor_hor_obs = mor_hor_obs, mor_hor_est = mor_hor_est))
    }
    
    sfInit(parallel = TRUE, cpus = cpus, type = "SOCK")
    sfExport("pair"); sfExport("com")
    sfExport("loop"); sfLibrary(iNEXT.beta3D); sfLibrary(snowfall)
    
    res <- sfLapply(pair$pos, fun = loop, com = com)
    sfStop()
    
  } else {
    
    loop <- function(v, com = com){
      
      print(v)
      
      out <- iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), level = SC, nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0)[[1]]
      
      out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), base = "size", 
                            level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0)[[1]]
      beta.obs = out.obs$gamma$Gamma / out.obs$alpha$Alpha
      
      
      ## Observed Jaccard similarity q = 0 (U02) index 
      jac_obs = 1 - (1/beta.obs[1] - 1) / (1/2 - 1)
      
      ## Standardized Jaccard similarity q = 0 (U02) index
      jac_est = 1 - out$`1-U`$Dissimilarity[1]
      
      ## Observed Sorenson similarity q = 0 (C02) index
      sor_obs = 1 - (beta.obs[1] - 1) / (2 - 1)
      
      ## Standardized Sorenson similarity q = 0 (C02) index
      sor_est = 1 - out$`1-C`$Dissimilarity[1]
      
      ## Observed Horn similarity q = 1 (C12 = U12) index
      hor_obs = 1 - log(beta.obs[2]) / log(2)
      
      ## Standardized Horn similarity q = 1 (C12 = U12) index
      hor_est = 1 - out$`1-C`$Dissimilarity[2]
      
      ## Observed Morisita Horn similarity q = 2 (C22) index
      mor_hor_obs = 1 - (1/beta.obs[3] - 1) / (1/2 - 1)
      
      ## Standardized Morisita Horn similarity q = 2 (C22) index
      mor_hor_est = 1 - out$`1-C`$Dissimilarity[3]
      
      
      return(c(jac_obs = jac_obs,         jac_est = jac_est, 
               sor_obs = sor_obs,         sor_est = sor_est,
               hor_obs = hor_obs,         hor_est = hor_est, 
               mor_hor_obs = mor_hor_obs, mor_hor_est = mor_hor_est))
    }
    
    res <- lapply(pair$pos, loop, com = com)
  }
  
  # finalization
  res <- do.call(rbind, res)
  res[res>1] <- 1  
  
  return(res)
  
}

