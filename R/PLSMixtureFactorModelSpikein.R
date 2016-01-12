PLSMixtureFactorModelSpikein <- function(Expression, Spikein, k1, k2, iniL = 3, TruncateL = 20, 
                        iter = 1000, maxK = 20, nu0 = 50, r = 0.1, s = 200, alpha = 1, alpha2 = 5,
                        kappa0 = 1, m = 0.75, g = 1, c = 1, d = 1, kappa = 2, diagH = 1, h1 = 1, h2 = 1, s1 = 1, s2 = 1, 
                        iter_to_average = 10, method = c("SpikeSlab", "IBP"), mergeCluster = TRUE, mergeN = 5, average_label = TRUE){
  
  method <- match.arg(method)
  X <- Expression 
  Y <- Spikein
  sigma <- sd(X)/2
  mu0 <- rep(0, k2)
  Sigma0 <- diag(k2)
  n <- nrow(X)
  if(nrow(Y) != n) stop("Spike-in and gene expression matrix should have the same number of samples (rows).")
  
  q <- ncol(X)
  p <- ncol(Y)
  
  if (method == "SpikeSlab") {
    res <- DirichletSpikePLSModel(Y = Y, X = X, k1 = k1, K = k2, iniL = iniL, TruncateL = TruncateL, 
                                  iter = iter, nu0 = nu0, sigma = sigma, r = r, s = s, alpha = alpha,
                                  mu0 = mu0, Sigma0 = Sigma0, kappa0 = kappa0, m = m, g = g, 
                                  c = c, d = d, diagH = diagH, h1 = h1, h2 = h2, s1 = s1, s2 = s2, iter_to_average = iter_to_average)
  } else if (method == "IBP") {
    res <- DirichletIBPPLSModel(Y = Y, X = X, k1 = k1, K = k2, iniL = iniL, TruncateL = TruncateL, 
                                iter = iter, maxK = maxK, nu0 = nu0, sigma = sigma, r = r, s = s, alpha = alpha, 
                                alpha2 = alpha2, mu0 = mu0, Sigma0 = Sigma0, kappa0 = kappa0, m = m, g = g, 
                                c = c, d = d, kappa = kappa, diagH = diagH, h1 = h1, h2 = h2, s1 = s1, s2 = s2, iter_to_average = iter_to_average)
  }
  
  Q <- res$Q
  
  if(average_label == TRUE){
    ClusterAssign <- apply(res$AveC, 2, which.max)
    ClusterList <- unique(ClusterAssign)
  } else {
    ClusterAssign <- res$C + 1
    ClusterList <- unique(ClusterAssign) 
  }
  
  if (method == "IBP"){
    activeK <- ncol(Q)
    MuList <- res$Mu[ClusterList, 1:activeK]
    SigmaList <- res$Sigma[1:activeK, 1:activeK, ClusterList]  
  } else if (method == "SpikeSlab") {
    MuList <- res$Mu[ClusterList, ]
    SigmaList <- res$Sigma[, , ClusterList]
  }
  
  ClusterLabel <- NULL
  for(i in 1:n){
    ClusterLabel <- c(ClusterLabel, which(ClusterList == ClusterAssign[i]))
  }
  
  BeforeMergeClusterLabel <- ClusterLabel
    
  l <- length(ClusterList)  
  
  if(mergeCluster == TRUE){ 
    mergelist <- which(table(ClusterLabel) <= mergeN)
    if(length(mergelist) > 0){
      currentLabel <- ClusterLabel
      currentMu <- MuList[-mergelist, ]
      currentCluster <- c(1:l)[-mergelist]
      for(j in 1:length(mergelist)){  
        distance <- apply(currentMu, 1, function(x){
          sum((x - MuList[mergelist[j], ])^2) 
        })  
        currentLabel[currentLabel == mergelist[j]] <- currentCluster[which.min(distance)]
      }
      ClusterLabel <- NULL
      for(i in 1:n){
        ClusterLabel <- c(ClusterLabel, which(currentCluster == currentLabel[i]))
      }
      MuList <- currentMu
      SigmaList <- SigmaList[, , -mergelist]
    }    
  }
  
  lambda <- res$lambdaU
  denoised <- Q%*%lambda
  
  E_x <- X - denoised - res$Z%*%res$lambdaX
  E_y <- Y - res$Z%*%res$lambdaY 
  
  psi_y <- res$psi[1:p]
  psi_x <- res$psi[-(1:p)]
  
  Likelihood <- 0
  for(i in 1:n){
    Likelihood <- Likelihood + dmvnrmRowArma(Q[i, ], MuList[ClusterLabel[i], ], SigmaList[, , ClusterLabel[i]], TRUE) 
  }
  for(j in 1:p){
    Likelihood <- Likelihood + sum(dnorm(c(E_y[, j]), 0, psi_y[j], TRUE))
  }
  for(j in 1:q){
    Likelihood <- Likelihood + sum(dnorm(c(E_x[, j]), 0, psi_x[j], TRUE))
  }
  Likelihood <- Likelihood + sum(dnorm(c(lambda), 0, res$sigma, TRUE))
  
  Ex_Likelihood <- 0
  for(j in 1:q){
    Ex_Likelihood <- Ex_Likelihood + sum(dnorm(c(E_x[, j]), 0, psi_x[j], TRUE))
  }
  
  adjustRes <- .identifiablityAdjust(res$Q, res$lambdaU)
  Q <- adjustRes[[1]]
  Lambda <- adjustRes[[2]]
  report <- list(Lambda, Q, ClusterLabel, denoised, Likelihood, MuList, Ex_Likelihood, BeforeMergeClusterLabel)
  names(report) <- c("SparseFactor", "Loading", "Cluster", "Denoised", "Likelihood", "Mu", "Ex_Likelihood", "ClusterBeforeMerging")
  class(report) <- "CitrusReport"
  return(report)
}
