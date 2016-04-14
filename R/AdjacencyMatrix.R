AdjacencyMatrix <- function(x){
  
  pp <- unique(x)
  n <- length(x)
  all <- matrix(0, ncol = n, nrow = n)
  for(i in 1:length(pp)){
    IDs <- which(x == pp[i])	
    all[IDs, IDs] <- 1
  }
  diag(all) <- 0
  return(all)
  
}