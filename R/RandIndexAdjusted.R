RandIndexAdjusted <- function(x, y){
  
  k <- length(unique(x))	
  l <- length(unique(y))
  n <- length(x)
  
  cluster1 <- unique(x)
  cluster2 <- unique(y)
  
  tij_sum <- 0
  
  for(i in 1:k){
    for(j in 1:l){
      set1 <- which(x == cluster1[i])
      set2 <- which(y == cluster2[j])
      set12 <- intersect(set1, set2)
      tij <- length(set12)
      if( tij > 1){
        tij_sum <- tij_sum + tij*(tij-1)/2
      }				
    }
  }
  
  ti_sum <- 0
  for(i in 1:k){
    set1 <- which(x == cluster1[i])
    ti <- length(set1)
    if(ti > 1){
      ti_sum <- ti_sum + ti*(ti-1)/2
    }		
  }
  
  tj_sum <- 0
  for(j in 1:l){
    set2 <- which(y == cluster2[j])
    tj <- length(set2)
    if(tj > 1){
      tj_sum <- tj_sum + tj*(tj-1)/2
    }		
  }
  
  n2 <- n*(n-1)/2
  
  expected <- ti_sum*tj_sum/n2
  
  RI <- (tij_sum - expected)/((ti_sum + tj_sum)/2 - expected) 
  
  RI <- ifelse(RI < 0, 0, RI)
  
  return( round(RI, 2))
  
} 

