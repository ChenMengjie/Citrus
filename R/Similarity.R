Similarity <- function(dataMat, cluster, ...){
  if(nrow(dataMat) != length(cluster))
    stop("The number of rows of dataMat should be equal to the length of cluster.")

  res <- EvaluateSimilarity(dataMat, cluster)
  names(res) <- c("WithinCluster", "BetweenCluster")
  return(res)
}
