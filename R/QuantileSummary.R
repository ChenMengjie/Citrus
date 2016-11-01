
QuantileSummary <- function(VarTable, quantiles = seq(0.1, 1, by = 0.1), rankingby){

per.quantile <- NULL
for(i in 1:length(quantiles)){
  per.quantile <- c(per.quantile, quantile(rankingby, quantiles[i]))
}

gene.quantile <- NULL
for(i in 1:length(quantiles)){
  if(i == 1){
    gene.quantile <- rbind(gene.quantile, apply(VarTable[rankingby < per.quantile[i], ], 2, mean))
  } else if(i == length(quantiles)){
    gene.quantile <- rbind(gene.quantile, apply(VarTable[rankingby >= per.quantile[i-1], ], 2, mean))		
  } else {
    gene.quantile <- rbind(gene.quantile, apply(VarTable[rankingby >= per.quantile[i-1] & rankingby < per.quantile[i], ], 2, mean))
  }
}

data <- data.frame(gene.quantile, quantiles*100)
colnames(data)[ncol(data)] <- "Quantile"
return(data)
}

