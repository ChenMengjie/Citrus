#' This is function to do factor analysis.
#' @param Expression A matrix..

FactorAnalysisWithSpikein <- function(Expression, Spikein, k1, k2, iter = 500, alpha = 10, 
                            kappa = 2, diagH = 1, g = 1, h1 = 1, h2 = 1, c = 1, d = 1, limit = 25,
                            method = c("EM", "EMSparse", "EMSparseTraining", "EMSparseNfold", "IBP", "PCA"), 
                            penalty = seq(0.01, 10, length.out = 20), tol = 0.1, givenpenalty = NULL, kfold = NULL,
                            Chunk = TRUE, chunk.size = 1000) {
  method <- match.arg(method)
  X <- Expression 
  Y <- Spikein
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  
  if (method == "EM") {
    if(Chunk == TRUE & q > chunk.size){
      chunk <- ceiling(q/chunk.size)
      res <- PLSfactorEMchunk(X, Y, k1, k2, iter, chunk)
    } else {
      res <- PLSfactorEM(X, Y, k1, k2, iter)
    }
    lambda <- res$lambda
    lambdaU <- lambda[(k1+1):(k1+k2), (p+1):(p+q)]
    lambdaX <- lambda[1:k1, (p+1):(p+q)]
    lambdaY <- lambda[1:k1, 1:p]
    Z <- as.matrix(res$Z[, 1:k1])
    Q <- as.matrix(res$Z[, (k1+1):(k1+k2)])
    denoised <- Q%*%lambdaU
    
  } else if (method == "EMSparseTraining"){
    
    res <- PLSfactorEMpenaltyTrain(Y, X, k1, k2, iter, penalty, tol)
    lambdaU <- res$lambdaU
    lambdaX <- res$lambdaX
    lambdaY <- res$lambdaY
    Z <- as.matrix(res$Z[, 1:k1])
    Q <- as.matrix(res$Z[, (k1+1):(k1+k2)])
    denoised <- Q%*%lambdaU
    
  } else if (method == "EMSparse"){
    
    if(is.null(givenpenalty)) stop("Method \"EMSparse\" requires the specification of penalty for sparisty.")
    res <- PLSfactorEMpenaltyGivenPen(Y, X, k1, k2, iter, givenpenalty, tol)
    lambdaU <- res$lambdaU
    lambdaX <- res$lambdaX
    lambdaY <- res$lambdaY
    Z <- as.matrix(res$Z[, 1:k1])
    Q <- as.matrix(res$Z[, (k1+1):(k1+k2)])
    denoised <- Q%*%lambdaU
    
  } else if (method == "EMSparseNfold"){
   
    if(is.null(kfold)) stop("Method \"EMSparseNfold\" requires the specification of k-fold cross validation.")
    res <- PLSfactorEMpenaltyNfoldCV(Y, X, k1, k2, iter, kfold, penalty, tol)
    lambdaU <- res$lambdaU
    lambdaX <- res$lambdaX
    lambdaY <- res$lambdaY
    Z <- as.matrix(res$Z[, 1:k1])
    Q <- as.matrix(res$Z[, (k1+1):(k1+k2)])
    denoised <- Q%*%lambdaU
    
  } else if (method == "IBP") {
    
    sigma <- sd(X)/2
    res <- PLSIBPfactormodel(Y, X, k1, k2, iter, sigma = sigma, alpha = alpha, kappa = kappa, diagH = diagH, 
                             g = g, h1 = h1, h2 = h2, c = c, d = d, limit = limit)
    lambdaU <- res$lambdaU
    lambdaX <- res$lambdaX
    lambdaY <- res$lambdaY
    Z <- as.matrix(as.matrix(res$Z))
    Q <- as.matrix(as.matrix(res$Q))
    denoised <- Q%*%lambdaU
    
  } else if (method == "PCA") {
    
    res <- NaivePCA(X, Y, k1, k2)
    lambda <- res$lambda
    lambdaU <- lambda[(k1+1):(k1+k2), (p+1):(p+q)]
    lambdaX <- lambda[1:k1, (p+1):(p+q)]
    lambdaY <- lambda[1:k1, 1:p]
    Z <- as.matrix(res$Z[, 1:k1])
    Q <- as.matrix(res$Z[, (k1+1):(k1+k2)])
    denoised <- Q%*%lambdaU
    
  }
  
  if(method != "PCA"){    
    E_x <- X - denoised - Z%*%lambdaX
    psi_x <- res$psi[-(1:p)]
    psi_x[psi_x < 0.000001] <- 0.000001 
    Ex_Likelihood <- 0
    for(j in 1:q){
      #if(is.na(sum(dnorm(c(E_x[, j]), 0, psi_x[j], TRUE)))) print(j)
      Ex_Likelihood <- Ex_Likelihood + sum(dnorm(c(E_x[, j]), 0, psi_x[j], TRUE))
    }
    report <- list(lambdaU, Q, denoised, Ex_Likelihood, Z, lambdaX, lambdaY, psi_x)
    names(report) <- c("Factor", "Loading", "Denoised", "Ex_Likelihood", "Z", "lambdaX", "lambdaY", "psi_x")
    
  } else {
    
    report <- list(lambdaU, Q, denoised, Z, lambdaX, lambdaY)
    names(report) <- c("Factor", "Loading", "Denoised", "Z", "lambdaX", "lambdaY")
    
  }
  return(report)
}