# DivPair Example, version 2
# 2024-11-04

rm(list=ls(all=TRUE))

# Pairwise similarity indices by Anne Chao
# to be loaded in analysis script by "load()"
# to provide iNEXTbeta3D_pair functions
# updated by Oliver, 15.04.2024, added phylogenetic and function diversity
# updated by Oliver and Mareike, 01.10.2024, added incidence and restructured output

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
#' OLIVER: for incidence data 2 more columns at the beginning are expected
#' column 1: sample id (column name arbitrary)
#' column 2: assemblage id (column name "plot")
#' Structure of "com" should be brought in better agreement with case "abundance" later

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

#' OLIVER: added parameters
#' div0="TD": to distinguish between diversity types, options: "TD", "PD", "FD", as in iNEXT3D
#' PDTree0=NULL: phylogenetic tree in case of div0="PD", as in iNEXT3D
#' FDdistM0=NULL: functional distance matrix for case div0="FD", as in iNEXT3D
#' datatype0 = "abundance": can be changed to "incidence_raw" as in iNEXT3D

###########
# extended function for parallelisation
iNEXTbeta3D_pair3D <- function(com, SC = NULL, parallel = F, cpus = 2, div0="TD", PDTree0=NULL, FDdistM0=NULL, datatype0 = "abundance"){  
  
  # data preprocessing for data type "abundance"
  if (datatype0 == "abundance"){
    m <- matrix(NA, ncol = dim(com)[1], nrow = dim(com)[1], byrow = T)
    
    # generate table of assemblage pairs, one pair per row
    # pos: row index of pair
    # x: assamblage 1 of specific pair
    # y: assamblage 2 of specific pair
    pair <- data.frame(pos = 1:dim(com)[1]^2,
      x = rep(rownames(com), nrow(m)), 
      y = as.vector(sapply(rownames(com), rep, nrow(m))))
  }
  
  
  # data preprocessing for data type "incidence_raw"
  if (datatype0 == "incidence_raw"){

    # expecting 2 additional columns at the beginning of dataframe "com":
    # column 1: sample id
    # column 2: named "plot" in com, assemblage id
    
    # Convert to list
    inci.raw = lapply(unique(com$plot), function(i) com %>% filter(plot == i) %>% .[,-(1:2)] %>% t)
    names(inci.raw) = unique(com$plot)
  
    # Generate table with all pairs of assemblages
    grid<-expand.grid(names(inci.raw),names(inci.raw), stringsAsFactors = F)
    # Remove doubled pairs (but keep pairs of identical assemblages)
    #grid<-grid[grid$Var1 >= grid$Var2,]
    # Combine by column, as in case of "abundance" and rename columns
    pair<-cbind(pos = 1:nrow(grid),grid)
    names(pair)[2]<-"x"
    names(pair)[3]<-"y"
  }
  
  if (datatype0 == "abundance"){
    if (is.null(SC)) SC = sapply(pair$pos, function(i) DataInfobeta3D( t(com[c(pair$x[i], pair$y[i]),]), diversity=div0, datatype=datatype0)$`SC(2n)`[4]) %>% min
  }
  
  if (datatype0 == "incidence_raw"){
    if (is.null(SC)) SC = sapply(pair$pos, function(i) DataInfobeta3D( inci.raw[c(pair$x[i], pair$y[i])], diversity=div0, datatype=datatype0)$`SC(2T)`[4]) %>% min
  }
  
  # for parallelization
  if (parallel == TRUE){
    
    loop <- function(v, com = com){
      #print(div0)
      
      sfCat(paste(v))
      
      if (datatype0 == "abundance") tt<-t(com[c(pair$x[v], pair$y[v]),])
      if (datatype0 == "incidence_raw"){
        tt<-list(inci.raw[pair$x[v]],inci.raw[pair$y[v]])
        tt<-inci.raw[c(pair$x[v], pair$y[v])]
      }
      
      out <- iNEXTbeta3D(tt, diversity=div0, q = c(0, 1, 2), level = SC, nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0, datatype=datatype0)[[1]]
      
      #out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), base = "size", 
      #                      level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0, datatype=datatype0)[[1]]
      
      if (datatype0 == "abundance"){
        out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), base = "size", 
          level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0)[[1]]
      }
      
      # OLIVER
      # Anne, please check here the "level" parameter: here "level" for incidence observed data should be specified as the number of sampling units
      if (datatype0 == "incidence_raw"){
        out.obs = iNEXTbeta3D(tt, diversity=div0, q = c(0, 1, 2), base = "size", 
                              level = ncol(tt[[1]]), nboot = 0, PDtree = PDTree0, 
                              FDdistM=FDdistM0, datatype=datatype0)[[1]]
      }
      
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
    
    # OLIVER
    # some more objects have to be exported to the parallisation environment
    sfExport("PDTree0"); sfExport("FDdistM0")
    sfExport("loop"); sfLibrary(iNEXT.beta3D); sfLibrary(snowfall)
    
    res <- sfLapply(pair$pos, fun = loop, com = com)
    sfStop()
    
  # without parallelisation  
  } else {
    
    loop <- function(v, com = com){
      
      #print(paste(v,"out of",length(pair[,1])))
      
      if (datatype0 == "abundance") tt<-t(com[c(pair$x[v], pair$y[v]),])
      if (datatype0 == "incidence_raw"){
        tt<-list(inci.raw[pair$x[v]],inci.raw[pair$y[v]])
        tt<-inci.raw[c(pair$x[v], pair$y[v])]
        #print(tt)
      }
      
      #out <- iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), level = SC, nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0, datatype=datatype0)[[1]]
      
      #out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), base = "size", 
      #                      level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0, datatype=datatype0)[[1]]
      out <- iNEXTbeta3D(tt, diversity=div0, q = c(0, 1, 2), level = SC, nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0, datatype=datatype0)[[1]]

      #out.obs = iNEXTbeta3D(tt, diversity=div0, q = c(0, 1, 2), base = "size", 
      #                      level = sum(unlist(inci.raw[c(pair$x[v], pair$y[v])])), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0, datatype=datatype0)[[1]]
      
      if (datatype0 == "abundance"){
        out.obs = iNEXTbeta3D(t(com[c(pair$x[v], pair$y[v]),]), diversity=div0, q = c(0, 1, 2), base = "size", 
                              level = sum(com[c(pair$x[v], pair$y[v]),]), nboot = 0, PDtree = PDTree0, FDdistM=FDdistM0)[[1]]
      }
      
      # OLIVER:
      # Anne, please check here the "level" parameter: here "level" for incidence observed data should be specified as the number of sampling units
      if (datatype0 == "incidence_raw"){
        out.obs = iNEXTbeta3D(tt, diversity=div0, q = c(0, 1, 2), base = "size", 
                              level = ncol(tt[[1]]), nboot = 0, PDtree = PDTree0, 
                              FDdistM=FDdistM0, datatype=datatype0)[[1]]
      }
      
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
    
    #loop(2,com)
    
    res <- lapply(pair$pos, loop, com = com)
  }
  
  # finalization
  res <- do.call(rbind, res)
  res[res>1] <- 1
  # OLIVER: correct 2nd error, too
  res[res<0] <- 0
  
  # OLIVER: output restructured to matrix form for flexibility
  if (datatype0 == "incidence_raw") n<-length(unique(com$plot))
  if (datatype0 == "abundance") n<-length(unique(rownames(com)))

  # Initialize empty list to collect distance matrices
  matlist<-list()
  # Loop through all columns of the results
  for (ii in 1:8) {
    mm<-matrix(res[,ii],ncol=n, nrow=n)
    
    # add row and column names to output matrix
    rownames(mm)<-unique(pair$x)
    colnames(mm)<-unique(pair$x)
    
    matlist[[ii]]<-mm
  }
  # Add names to matrices in list
  names(matlist)<-colnames(res)

  #return(list(res,pair,mm))
  return(list(ResultTable=res,Pairs=pair,Matrices=matlist))
  
}


