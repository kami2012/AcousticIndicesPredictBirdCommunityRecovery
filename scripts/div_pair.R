library(SpadeR)
library(snowfall)
library(vegan)

div_comp <- function(a,b, nboot = 3){
  # a <- x[to$y[v],]; b <- x[to$x[v],]; nboot <- 3
  # a <- v[1,] ; b <- v[2,]
  a <- as.numeric(a)
  b <- as.numeric(b)
  
  #jaccard and sorenson
  
  k <- which(a > 0 & b > 0)

  #SpadeR:::SimilarityPair 
  #SpadeR:::Jaccard_Sorensen_Abundance_equ
  
  s.sp <- SpadeR:::PanEstFun(a,b)
  a.sp <- SpadeR:::SpecAbunChao1(a, k = 10, conf = 0.95)[1,1] + SpadeR:::SpecAbunChao1(b, k = 10, conf = 0.95)[1,1]
  
  res <- data.frame(jac_obs = length(k)/sum(a + b > 0), 
                    jac_est = s.sp/(a.sp - s.sp), 
                    sor_obs = 2 * length(k)/(sum(a > 0) + sum(b > 0)),
                    sor_est = 2 * s.sp/a.sp)
  #js <- SpadeR:::Jaccard_Sorensen_Abundance_equ(X1 = a, X2 = b, boot = 30, datatype = "abundance")
  
  # horn
  js <- SpadeR:::Horn.Est(cbind(a,b),"abundance")
  res$hor_obs <- js[2]
  res$hor_est <- js[1]
  
  # morisita horn
  js <- SimilarityThree(X = cbind(a,b), datatype = "abundance", q = 2, method = "equal weight",nboot=0)
  res$mor_hor_obs <- js$CqN[1,1]
  res$mor_hor_est <- js$CqN[2,1]
  cat(".")
  return(res)
}

#data(dune)
#rownames(dune) <- paste0("p_", rownames(dune))

div_pair <- function(com, nboot = 3, parallel = F, cpus = 2){
  # com <- dune; nboot = 3
  m <- matrix(NA, ncol = dim(com)[1], nrow = dim(com)[1], byrow = T)
  # rownames(m) <- colnames(m) <- row.names(com)
  
  to <- data.frame(pos = 1:dim(com)[1]^2,
                   x = rep(rownames(com), nrow(m)), 
                   y = as.vector(sapply(rownames(com), rep, nrow(m))) )
  
  # for parallelization
  if(parallel == TRUE){
    
    loop <- function(v, com = com){
      # v <- 135
      sfCat(paste(v))
      v <- div_comp( com[as.character(to$x[v]),], com[as.character(to$y[v]),], nboot = nboot)
      return(v)
    }
    
    sfInit(parallel=TRUE, cpus = cpus, type = "SOCK")
    sfExport("div_comp"); sfExport("div_pair"); sfExport("to"); sfExport("com")
    sfExport("loop"); sfExport("SimilarityThree"); sfLibrary(snowfall)
    res <- sfLapply(to$pos, fun = loop, com = com)
    sfStop()
  }else{
    
    loop <- function(v, com = com){
      # v <- 2
      print(v)
      v <- div_comp( com[as.character(to$x[v]),], com[as.character(to$y[v]),], nboot = nboot)
      return(v)
    }
    
    res <- lapply(to$pos, loop, com = com)
  }
  
  # finalization
  res <- do.call(rbind, res)
  res[res>1] <- 1  
  
  return(res)
} # end div pair function

# div_pair(dune) 


#### bootstrap hack

SimilarityThree<-function (X, q, nboot = 50, datatype = "abundance", method = c("equal weight", 
                                                                                "unequal weight")) 
{
  if (datatype == "abundance") {
    p <- SpadeR:::Two_com_correct_obspi(X[, 1], X[, 2])
  }
  else {
    p <- SpadeR:::Two_com_correct_obspi_Inc(X[, 1], X[, 2])
    t <- X[1, ]
    X <- X[-1, ]
  }
  N = ncol(X)
  ni = colSums(X)
  n = sum(X)
  pool = rowSums(X)
  bX = apply(X, 2, function(x) x/sum(x))
  pool.x = rowSums(bX)/N
  if (q == 0) {
    f1 = apply(X, 2, function(x) sum(x == 1))
    f2 = apply(X, 2, function(x) sum(x == 2))
    Sobs = apply(X, 2, function(x) sum(x > 0))
    Si = Sobs + sapply(1:N, function(k) ifelse(f2[k] == 
                                                 0, f1[k] * (f1[k] - 1)/2, f1[k]^2/(2 * f2[k])))
    Sa = mean(Si)
    UqN.mle = (1/N - mean(Sobs)/sum(pool > 0))/(1/N - 1)
    CqN.mle = (N - sum(pool > 0)/mean(Sobs))/(N - 1)
    F1 = sum(pool == 1)
    F2 = sum(pool == 2)
    Sg = sum(pool > 0) + ifelse(F2 == 0, F1 * (F1 - 1)/2, 
                                F1^2/(2 * F2))
    UqN = min(1, (1/N - Sa/Sg)/(1/N - 1))
    UqN = max(0, UqN)
    CqN = min(1, (N - Sg/Sa)/(N - 1))
    CqN = max(0, CqN)
    b.UqN = numeric(nboot)
    b.UqN.mle = numeric(nboot)
    b.CqN = numeric(nboot)
    b.CqN.mle = numeric(nboot)
    
    se.U = sd(b.UqN)
    se.U.mle = sd(b.UqN.mle)
    se.C = sd(b.CqN)
    se.C.mle = sd(b.CqN.mle)
    out1 = rbind(c(UqN.mle, se.U.mle, min(1, UqN.mle + 1.96 * 
                                            se.U.mle), max(0, UqN.mle - 1.96 * se.U.mle)), c(UqN, 
                                                                                             se.U, min(1, UqN + 1.96 * se.U), max(0, UqN - 1.96 * 
                                                                                                                                    se.U)))
    out2 = rbind(c(CqN.mle, se.C.mle, min(1, CqN.mle + 1.96 * 
                                            se.C.mle), max(0, CqN.mle - 1.96 * se.C.mle)), c(CqN, 
                                                                                             se.C, min(1, CqN + 1.96 * se.C), max(0, CqN - 1.96 * 
                                                                                                                                    se.C)))
  }
  if (q == 2) {
    if (method == "equal weight") {
      a.mle = N/sum(bX^2)
      g.mle = 1/sum(pool.x^2)
      b.mle = g.mle/a.mle
      UqN.mle = (N - b.mle)/(N - 1)
      CqN.mle = (1/N - 1/b.mle)/(1/N - 1)
      Ai = sapply(1:N, function(k) sum(X[, k] * (X[, k] - 
                                                   1)/(ni[k] * (ni[k] - 1))))
      bX.1 = apply(X, 2, function(x) (x - 1)/(sum(x) - 
                                                1))
      temp = sapply(1:nrow(X), function(j) (sum(bX[j, 
                                                   ] %*% t(bX[j, ])) - sum(bX[j, ]^2)) + sum(bX[j, 
                                                                                                ] * bX.1[j, ]))
      G = 1/(sum(temp)/N^2)
      B = G/(1/mean(Ai))
      UqN = min(1, (N - B)/(N - 1))
      UqN = max(0, UqN)
      CqN = min(1, (1/N - 1/B)/(1/N - 1))
      CqN = max(0, CqN)
      b.UqN = numeric(nboot)
      b.UqN.mle = numeric(nboot)
      b.CqN = numeric(nboot)
      b.CqN.mle = numeric(nboot)
      
      se.U = sd(b.UqN)
      se.U.mle = sd(b.UqN.mle)
      se.C = sd(b.CqN)
      se.C.mle = sd(b.CqN.mle)
      out1 = rbind(c(UqN.mle, se.U.mle, min(1, UqN.mle + 
                                              1.96 * se.U.mle), max(0, UqN.mle - 1.96 * se.U.mle)), 
                   c(UqN, se.U, min(1, UqN + 1.96 * se.U), max(0, 
                                                               UqN - 1.96 * se.U)))
      out2 = rbind(c(CqN.mle, se.C.mle, min(1, CqN.mle + 
                                              1.96 * se.C.mle), max(0, CqN.mle - 1.96 * se.C.mle)), 
                   c(CqN, se.C, min(1, CqN + 1.96 * se.C), max(0, 
                                                               CqN - 1.96 * se.C)))
    }
    if (method == "unequal weight") {
      a.mle = 1/(N * sum((X/n)^2))
      g.mle = 1/sum((pool/n)^2)
      b.mle = g.mle/a.mle
      UqN.mle = (N - b.mle)/(N - 1)
      CqN.mle = (1/N - 1/b.mle)/(1/N - 1)
      A = (1/N) * (1/sum(X * (X - 1)/(n * (n - 1))))
      G = 1/sum(pool * (pool - 1)/(n * (n - 1)))
      B = G/A
      UqN = min(1, (N - B)/(N - 1))
      UqN = max(0, UqN)
      CqN = min(1, (1/N - 1/B)/(1/N - 1))
      CqN = max(0, CqN)
      b.UqN = numeric(nboot)
      b.UqN.mle = numeric(nboot)
      b.CqN = numeric(nboot)
      b.CqN.mle = numeric(nboot)
      for (i in 1:nboot) {
        if (datatype == "abundance") {
          XX = sapply(1:N, function(k) rmultinom(1, 
                                                 ni[k], p[, k]))
        }
        else {
          XX = sapply(1:N, function(k) sapply(1:nrow(p), 
                                              FUN = function(i) rbinom(1, t[1, k], p[i, 
                                                                                     1])))
        }
        pool = rowSums(XX)
        a.mle = 1/(N * sum((XX/n)^2))
        g.mle = 1/sum((pool/n)^2)
        b.mle = g.mle/a.mle
        b.UqN.mle[i] = (N - b.mle)/(N - 1)
        b.CqN.mle[i] = (1/N - 1/b.mle)/(1/N - 1)
        A = (1/N) * (1/sum(XX * (XX - 1)/(n * (n - 1))))
        G = 1/sum(pool * (pool - 1)/(n * (n - 1)))
        B = G/A
        b.UqN[i] = (N - B)/(N - 1)
        b.CqN[i] = (1/N - 1/B)/(1/N - 1)
      }
      se.U = sd(b.UqN)
      se.U.mle = sd(b.UqN.mle)
      se.C = sd(b.CqN)
      se.C.mle = sd(b.CqN.mle)
      out1 = rbind(c(UqN.mle, se.U.mle, min(1, UqN.mle + 
                                              1.96 * se.U.mle), max(0, UqN.mle - 1.96 * se.U.mle)), 
                   c(UqN, se.U, min(1, UqN + 1.96 * se.U), max(0, 
                                                               UqN - 1.96 * se.U)))
      out2 = rbind(c(CqN.mle, se.C.mle, min(1, CqN.mle + 
                                              1.96 * se.C.mle), max(0, CqN.mle - 1.96 * se.C.mle)), 
                   c(CqN, se.C, min(1, CqN + 1.96 * se.C), max(0, 
                                                               CqN - 1.96 * se.C)))
    }
  }
  out1 <- cbind(out1[, c(1, 2)], out1[, 4], out1[, 3])
  colnames(out1) = c("UqN", "se", "95%.Lower", "95%.Upper")
  rownames(out1) = c("Emperical", "Estimate")
  out2 <- cbind(out2[, c(1, 2)], out2[, 4], out2[, 3])
  colnames(out2) = c("CqN", "se", "95%.Lower", "95%.Upper")
  rownames(out2) = c("Emperical", "Estimate")
  return(list(UqN = out1, CqN = out2))
}