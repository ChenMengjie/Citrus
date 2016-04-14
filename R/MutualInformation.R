MutualInformation <- function(x, y){
  
  k <- length(unique(x))	
  l <- length(unique(y))
  n <- length(x)
  
  cluster1 <- unique(x)
  cluster2 <- unique(y)
  
  MI <- 0
  for(i in 1:k){
    for(j in 1:l){
      set1 <- which(x == cluster1[i])
      set2 <- which(y == cluster2[j])
      set12 <- intersect(set1, set2)
      tij <- length(set12)
      if(tij > 0){
        MI <- MI + tij*log(tij*n/(length(set1)*length(set2)))/log(k*l)		
      }	
    }
  }
  
  return(round(MI/n, 2))
  
}

