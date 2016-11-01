MixtureFactorModelWithoutSpikein <- function(Expression, K, iniL = 3, ModelZero = FALSE, kappa_int = 1, kappagrid = seq(0.1, 2, by = 0.2),
                    TruncateL = 20, iter = 1000, maxK = 20, nu0 = 50, r = 0.1, s = 200, alpha = 1, alpha2 = 5,
                    kappa0 = 1, m = 0.75, g = 1, h = 1, c = 1, d = 1, kappa_ibp = 2, method = c("SpikeSlab", "IBP"), 
                    s1 = 1, s2 = 1, iter_to_average = 10, mergeCluster = TRUE, mergeN = 5, average_label = TRUE){
  
  method <- match.arg(method)
  
  X <- apply(Expression, 2, function(z){
    z - mean(z)
  }) 
  
  #tau <- apply(Expression, 2, mean)
  
  tau <- rep(0, ncol(Expression))
  
  sigma <- sd(X)/2
  mu0 <- rep(0, K)
  Sigma0 <- diag(K)
  n <- nrow(X)
  q <- ncol(X)
  
  if(ModelZero == TRUE){
    H <- matrix(0, ncol = ncol(X), nrow = nrow(X))
    H[Expression == 0] <- 1
    if(method == "SpikeSlab"){
      res <- DirichletSpikeModelZero(X = X, H = H, tau = tau, kappa_int = kappa_int, kappagrid = kappagrid, 
                                 K = K, iniL = iniL, TruncateL = TruncateL, iter = iter, 
                                 nu0 = nu0, sigma = sigma, r = r, s = s, alpha = alpha, mu0 = mu0, 
                                 Sigma0 = Sigma0, kappa0 = kappa0, m = m, g = g, h = h, c = c, d = d, 
                                 s1 = s1, s2 = s2, iter_to_average = iter_to_average)
      kappa <- res$kappa
    } else if (method == "IBP") {
      res <- DirichletIBPModelZero(X = X, H = H, tau = tau, kappa_int = kappa_int, kappagrid = kappagrid,
                               K = K, iniL = iniL, TruncateL = TruncateL, iter = iter, maxK = maxK, 
                               nu0 = nu0, sigma = sigma, r = r, s = s, alpha = alpha, alpha2 = alpha2, 
                               mu0 = mu0, Sigma0 = Sigma0, kappa0 = kappa0, m = m, g = g, h = h, 
                               c = c, d = d, kappa_ibp = kappa_ibp, s1 = s1, s2 = s2, iter_to_average = iter_to_average)
      kappa <- res$kappa
    }
  } else {
    if (method == "SpikeSlab") {
      res <- DirichletSpikeModel(X = X, K = K, iniL = iniL, TruncateL = TruncateL, iter = iter, 
                               nu0 = nu0, sigma = sigma, r = r, s = s, alpha = alpha, mu0 = mu0, 
                               Sigma0 = Sigma0, kappa0 = kappa0, m = m, g = g, h = h, c = c, d = d, 
                               s1 = s1, s2 = s2, iter_to_average = iter_to_average)
    
    } else if (method == "IBP") {
      res <- DirichletIBPModel(X = X, K = K, iniL = iniL, TruncateL = TruncateL, iter = iter, maxK = maxK, 
                             nu0 = nu0, sigma = sigma, r = r, s = s, alpha = alpha, alpha2 = alpha2, 
                             mu0 = mu0, Sigma0 = Sigma0, kappa0 = kappa0, m = m, g = g, h = h, 
                             c = c, d = d, kappa = kappa_ibp, s1 = s1, s2 = s2, iter_to_average = iter_to_average)
    }
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
  } else {
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
      if(!is.vector(currentMu)){
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
      }
      MuList <- currentMu
      SigmaList <- SigmaList[, , -mergelist]
    }    
  }
  
  lambda <- res$lambda
  
  denoised <- Q%*%lambda
  E_x <- X - denoised 
  
  Likelihood <- 0
  
  for(i in 1:n){
    test <- try(Likelihood <- Likelihood + dmvnrmRowArma(Q[i, ], MuList[ClusterLabel[i], ], SigmaList[, , ClusterLabel[i]], TRUE))
    if(class(test) == "try-error") {
      print(MuList)
      print(i)
      print(ClusterLabel)
    }
  }
  psi_x <- res$psi
  for(j in 1:q){
    Likelihood <- Likelihood + sum(dnorm(c(E_x[, j]), 0, psi_x[j], TRUE))
  }
  Likelihood <- Likelihood + sum(dnorm(c(lambda), 0, res$sigma, TRUE))
  
  adjustRes <- .identifiablityAdjust(res$Q, res$lambda)
  Q <- adjustRes[[1]]
  Lambda <- adjustRes[[2]]
  if(ModelZero == TRUE){
    report <- list(Lambda, Q, ClusterLabel, denoised, Likelihood, MuList, BeforeMergeClusterLabel, kappa)
    names(report) <- c("SparseFactor", "Loading", "Cluster", "Denoised", "Likelihood", "Mu", "ClusterBeforeMerging", "DropoutRate")
  } else {
    report <- list(Lambda, Q, ClusterLabel, denoised, Likelihood, MuList, BeforeMergeClusterLabel)
    names(report) <- c("SparseFactor", "Loading", "Cluster", "Denoised", "Likelihood", "Mu", "ClusterBeforeMerging")
  }

  class(report) <- "CitrusReport"
  
  return(report)
}

