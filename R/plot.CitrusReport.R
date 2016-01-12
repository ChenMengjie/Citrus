#' A function to make plot for ...
plot.CitrusReport <- function(X, SampleColVector = NULL, save = FALSE, name = NULL, ...){
  colorlist <- c("red", "violet", "blue", "gold", "slateblue1", "limegreen",  "skyblue1", "khaki", 
                 "pink", "chocolate", "salmon",  "maroon", "cyan", "lavender", 
                 "purple", "yellow", "turquoise4")
  l <- length(unique(X$Cluster))
  if(l > length(colorlist)){
    colorlist <- rainbow(l)
  }

  require(gplots)
  palette.gr.marray <- colorRampPalette(c("blue", "white", "red"))(56)
  if(class(X) != "CitrusReport"){
    stop("Input must be a CitrusReport!")
  } else {
    Cluster <- X$Cluster
    Mu <- X$Mu
    reordered <- .reorderClusterByMu(Cluster, Mu)
    res <- .reorderAmat(Cluster, reordered, X$Denoised)
    newOrder <- .reorderAvec(Cluster, reordered)  
    
    if(!is.null(SampleColVector)){
      if(is.null(nrow(SampleColVector))){
        if(length(SampleColVector) != nrow(res)){
          stop("The length of SampleColVector is not equal to the number of samples in the data!")
        } else {  
          colorMat <- cbind(colorlist[Cluster[newOrder]], SampleColVector[newOrder])
        }
      } else {
        if(nrow(SampleColVector) != nrow(res)){
          stop("The number of rows of SampleColVector is not equal to the number of samples in the data!")
        } else {
          colorMat <- cbind(colorlist[Cluster[newOrder]], SampleColVector[newOrder, ])
        }
      }
    } else {
      colorMat <- matrix(colorlist[Cluster[newOrder]], ncol = 1)
    }
    if(save == TRUE){
      if(is.null(name)){
        filename <-  paste0(getwd(), "/CitrusClusteringHeatmap.png")
      } else {
        filename <- name
      }
      png(filename, width = 10, height = 10, unit = "in", res = 400)  
      heatmap.3(as.matrix(t(res)), trace = "none", col = palette.gr.marray, margins = c(5, 5), symbreaks = T, 
                labRow = NA, labCol = NA, dendrogram = "none", ColSideColors = colorMat, Colv = FALSE, ...)
      dev.off()
    } else {
    heatmap.3(as.matrix(t(res)), trace = "none", col = palette.gr.marray, margins = c(5, 5), symbreaks = T, 
              labRow = NA, labCol = NA, dendrogram = "none", ColSideColors = colorMat, Colv = FALSE, ...)
    }
  }
}
