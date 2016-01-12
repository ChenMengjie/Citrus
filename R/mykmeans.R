mykmeans <- function(x, k){
  return(kmeans(x, k)$cluster)
}