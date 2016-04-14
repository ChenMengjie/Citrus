RandIndex <- function(x, y){
  
  n <- length(x)
  adjx <- AdjacencyMatrix(x)
  adjy <- AdjacencyMatrix(y)
  diff <- adjx - adjy	
  disgree <- length(diff[diff==1])
  RD <- 1 - disgree/(n*(n-1))
  return(round(RD, 2))
  
}
