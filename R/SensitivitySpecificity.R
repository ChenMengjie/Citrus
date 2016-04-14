SensitivitySpecificity <- function(x, y){
  
  n <- length(x)
  adjx <- AdjacencyMatrix(x)
  adjy <- AdjacencyMatrix(y)
  
  Aequal <- length(adjx[adjx==1])
  gequal <- length(adjx[adjx==1 & adjy==1])
  Sensitivity <- round(gequal/Aequal, 2)
  
  Anotequal <- length(adjx[adjx==0]) - n
  gnotequal <- length(adjx[adjx==0 & adjy==0]) - n
  Specificity <- round(gnotequal/Anotequal, 2)
  
  gequalprime <- length(adjy[adjy==1])
  Aequalprime <- length(adjy[adjy==1 & adjx==1])
  Positive <- round(Aequalprime/gequalprime, 2)
  
  gnotequalprime <- length(adjy[adjy==0])
  Anotequalprime <- length(adjy[adjy==0 & adjx==0])
  Negative <- round(Anotequalprime/gnotequalprime, 2)
  
  return(c(Sensitivity, Specificity, Positive, Negative))
  
}

