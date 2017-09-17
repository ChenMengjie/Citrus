clusterSummary <- function(clustering1, clustering2){
  all <- NULL
  all <- c(all, RandIndex(clustering1, clustering2))
  all <- c(all, SensitivitySpecificity(clustering1, clustering2))
  all <- c(all, MutualInformation(clustering1, clustering2))
  all <- c(all, RandIndexAdjusted(clustering1, clustering2))
  names(all) <- c("RandIndex", "Sensitivity", "Specificity", "PPV", "NPV", "Mutual", "RandIndexAdj")
  return(all)
}