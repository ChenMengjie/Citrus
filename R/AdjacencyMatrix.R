AdjacencyMatrix <- function(x){
  pp <- unique(x)
  n <- length(x)
  all <- matrix(0, ncol = n, nrow = n)
  for(i in 1:length(pp)){
    IDs <- which(x %in% pp[i])	
    all[IDs, IDs] <- 1
  }
  diag(all) <- 1
  return(all)
  
}