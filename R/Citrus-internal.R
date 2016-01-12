.identifiablityAdjust <- function(Q, lambda){
  k <- nrow(lambda)
  lambda_adjust <- lambda
  Q_adjust <- Q
  for(i in 1:k){
    aa <- lambda[i, ]
    aa.sd <- sd(aa[aa != 0])
    lambda_adjust[i, ] <- lambda[i, ]/aa.sd
    Q_adjust[, i] <- Q[, i]*aa.sd
  }
  return(list(Q_adjust, lambda_adjust))
}

.reorderClusterByMu <- function(C, Mu){
  
  current.list <- unique(C)
  current.Mu <- Mu[current.list, ]
  rownames(current.Mu) <- current.list
  iter <- length(current.list)
  current.cluster <- current.list[which.max(table(C))]
  reordered <- current.cluster 
  
  while(iter > 2){
    mu <- current.Mu[current.list == current.cluster, ]
    flag <- current.list != current.cluster
    current.Mu <- current.Mu[flag, ]
    current.list <- current.list[flag] 
    distance <- apply(current.Mu, 1, function(x){
      sum((x - mu)^2) 
    })	
    current.cluster <- current.list[which.min(distance)]
    reordered <- c(reordered, current.cluster)
    iter <- iter - 1
  }
  reordered <- c(reordered, current.list[current.list != current.cluster])
  return(reordered)	
}

.reorderAmat <- function(C, reordered, Z){
  len <- length(reordered)
  res <- NULL
  for(i in 1:len){
    subZ <- Z[which(C == reordered[i]), ]
    res <- rbind(res, subZ)
  }
  return(res)
}

.reorderAvec <- function(C, reordered){
  len <- length(reordered)
  res <- NULL
  for(i in 1:len){
    subZ <- which(C == reordered[i])
    res <- c(res, subZ)
  }
  return(res)
}
