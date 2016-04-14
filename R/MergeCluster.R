MergeCluster <- function(res, mergeN = 5){ 
  if(class(res) != "CitrusReport")
    stop("Input must be a CitrusReport!")
  
  n <- nrow(res$Loading)
  
  ClusterAssign <- res$Cluster
  ClusterList <- unique(ClusterAssign)

  MuList <- res$Mu
  ClusterLabel <- NULL
  for(i in 1:n){
    ClusterLabel <- c(ClusterLabel, which(ClusterList == ClusterAssign[i]))
  }
  
  BeforeMergeClusterLabel <- ClusterLabel
  l <- length(ClusterList)  

  mergelist <- which(table(ClusterLabel) <= mergeN)
  if(length(mergelist) > 0){
    currentLabel <- ClusterLabel
    currentMu <- MuList[-mergelist, ]
    if(!is.vector(currentMu)){
      currentCluster <- c(1:l)[-mergelist]
      for(j in 1:length(mergelist)){  
        distance <- apply(currentMu, 1, function(x){
          sum((x - MuList[mergelist[j], ])^2) 
        })  
        currentLabel[currentLabel == mergelist[j]] <- currentCluster[which.min(distance)]
      }
      ClusterLabel <- NULL
      for(i in 1:n){
        ClusterLabel <- c(ClusterLabel, which(currentCluster == currentLabel[i]))
      }
    }
    MuList <- currentMu
  }    
  
  res$Cluster <- ClusterLabel
  res$Mu <- MuList
  res$BeforeMergeClusterLabel <- BeforeMergeClusterLabel

  return(res)
}