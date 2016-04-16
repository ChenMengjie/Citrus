#' @title A function to remove unwanted variation from target genes given control genes.
#' @description To infer latent confounding factors from scRNAseq studies and remove unwanted variation, 
#' we develop a novel statistical method, which we refer to as scPLS. 
#' scPLS is based on the partial least squares regression models and incorporates both control and target genes to infer hidden confounding effects.
#' In addition, our method can model other systematic biological variation and heterogeneity, 
#' which are often observed in the target genes. By incorporating such systematic heterogeneity, 
#' we can further improve the estimation of the confounding factors and the removal of unwanted variation. 
#' To make our method widely applicable, we also develop a novel efficient estimation algorithm that is scalable 
#' to hundreds of cells and thousands of genes. 
#' @param Target: A $n \times q$ matrix for $q$ target genes from $n$ samples.
#' @param Control: A $n \times p$ matrix for $p$ control genes from $n$ samples. To remove cell cycle effect, input the matrix 
#' of cell cycle gene expression. To remove technical effect, input the matrix of control gene expression.  
#' @param k1: The number of technical factors.
#' @param k2: The number of structured biological factors.
#' @param iter: The number of iterations for the EM algorithm. The default is 500.
#' @param alpha: Hyperparameter if using the IBP prior, other hyperparameters including kappa, diagH, g, h1, h2, c, d.
#' @param limit: the maximal number of factors if using the IBP prior. 
#' @param method: The method used to infer the latent factors. The default is ``EM" algorithm.
#' Other choices include ``EMSparse" algorithm: penalty on sparsity of the factor matrix is specified. 
#' ``EMSparseTraining" algorithm: penalty on sparsity is learned from training samples.
#' ``EMparseNfold" algorithm: penalty on sparsity is learned from N fold cross-validation.
#' ``IBP" algorithm: the sparse factor matrix modeled by an IBP prior.
#' ``PCA" algorithm: The initializer for the EM algorithm. The latent factors are estimated from a Singular Value Decomposition.
#' @param penalty: A sequence of penality will be tested in training or cross-validation.    
#' @param tol: Tolerance for the convergence.
#' @param givenpenalty: Specified penalty level on sparsity of the factor matrix. The default is NULL. 
#' @param kfold: The fold number for cross-validation.
#' @param Chunk: Whether to use EM-in-chunks algorithm. The default is TRUE.
#' @param chunk.size: The chunk size (number of genes) for EM-in-chunks algorithm. The default is 1000.        


scPLS <- function(Target, Control, k1, k2, iter = 500, alpha = 10, 
                            kappa = 2, diagH = 1, g = 1, h1 = 1, h2 = 1, c = 1, d = 1, limit = 25,
                            method = c("EM", "EMSparse", "EMSparseTraining", "EMSparseNfold", "IBP", "PCA"), 
                            penalty = seq(0.01, 10, length.out = 20), tol = 0.1, givenpenalty = NULL, kfold = NULL,
                            Chunk = TRUE, chunk.size = 1000, center = TRUE) {
  
  Targetmean <- apply(Target, 2, mean)
  
  method <- match.arg(method)
  
  if(center == TRUE){
    X <- apply(Target, 2, function(z){
      z - mean(z)
    }) 
    Y <- apply(Control, 2, function(z){
      z - mean(z)
    }) 
  } else {
    X <- Target
    Y <- Control
  }
  
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
    report <- list(lambdaU, Q, Ex_Likelihood, Z, lambdaX, lambdaY, psi_x)
    names(report) <- c("Factor", "Loading", "Ex_Likelihood", "Z", "lambdaX", "lambdaY", "psi_x")
    
  } else {
  
    report <- list(lambdaU, Q, Z, lambdaX, lambdaY)
    names(report) <- c("Factor", "Loading",  "Z", "lambdaX", "lambdaY")
    
  }
  
  adjusted <- X - report$Z%*%report$lambdaX
  report$Method <- method
  if(center == TRUE){
    Adjusted <- apply(adjusted, 1, function(x){ x + Targetmean})
    report$Adjusted <- t(Adjusted)
  } else {
    report$Adjusted <- adjusted
  }
  
  if(is.vector(report$lambdaX)){
    XX_t <- t(matrix(report$lambdaX, nrow=1))%*%matrix(report$lambdaX, nrow=1)
  } else {
    XX_t <- t(report$lambdaX)%*%report$lambdaX
  }
  if(is.vector(report$Factor)){
    UU_t <- t(matrix(report$Factor, nrow=1))%*%matrix(report$Factor, nrow=1)
  } else {
    UU_t <- t(report$Factor)%*%report$Factor
  }
  
  Xcomponent <- diag(XX_t)
  Ucomponent <- diag(UU_t)
  varAll <- apply(Target, 2, var)
  varCon <- apply(report$Z%*%report$lambdaX, 2, var)
  if(method != "PCA"){
    varModel <- Xcomponent + Ucomponent + report$psi_x
    report$VarianceSummary <- data.frame(Sample = round(varAll, 3), 
                                       SampleConfounding = round(varCon, 3), 
                                       Model= round(varModel, 3), 
                                       ModelConfounding = round(Xcomponent, 3))
  } else {
    report$VarianceSummary <- data.frame(Sample = round(varAll, 3), 
                                         SampleConfounding = round(varCon, 3), 
                                         ModelConfounding = round(Xcomponent, 3))
  }
  return(report)
}