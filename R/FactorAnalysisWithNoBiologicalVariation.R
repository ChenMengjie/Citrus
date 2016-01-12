FactorAnalysisWithNoBiologicalVariation <- function(Expression, Spikein, k, iter = 500){
  
  X <- Expression 
  Y <- Spikein
  
  n <- nrow(X)
  q <- ncol(X)
  p <- ncol(Y)
  res <- factorEM(X, Y, k, iter)
  
  lambda <- as.matrix(res$lambda)
  lambdaX <- lambda[, (p+1):(p+q)]
  Z <- as.matrix(res$Z)
  E_x <- X - Z%*%lambdaX
  
  psi_x <- res$psi[-(1:p)]
  Ex_Likelihood <- 0
  for(j in 1:q){
    Ex_Likelihood <- Ex_Likelihood + sum(dnorm(c(E_x[, j]), 0, psi_x[j], TRUE))
  }
  report <- list(E_x, Ex_Likelihood, Z, lambdaX, psi_x)
  names(report) <- c("Residuals", "Ex_Likelihood", "Z", "lambdaX", "psi_x")
  return(report)
}