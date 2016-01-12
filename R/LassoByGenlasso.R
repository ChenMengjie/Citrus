#' Function to call lasso from genelasso package.
#' @param Y vector for response
#'  
LassoByGenlasso <- function(Y, X, D, lam){
  require(genlasso)
  res <- coef(genlasso(Y, X, D), lambda = lam)$beta
  return(res)
} 
