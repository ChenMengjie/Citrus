

plotBlock <- function(Z, sep.height = 5, sep.width = 5, Color = TRUE, supportOnly = TRUE){
  require(gplots)
  require(caret)
  
  N <- nrow(Z)
  P <- ncol(Z)
  csep <- seq(0, sep.height*(P-1), by = sep.height)
  
  if(supportOnly == TRUE){
    sepcolor <- apply(Z, 2, function(x){
      y <- rep("white", P)
      y[x != 0] <- "red"
      return(y)
    })
  } else if(Color == TRUE){
    max.value <- max(abs(Z))
    aa <- floor(Z*128/max.value)
    palette.gr.marray <- colorRampPalette(c("blue", "white", "red"))(256)
    colors_pos <- palette.gr.marray[129:256]
    colors_neg <- palette.gr.marray[1:128]
    sepcolor <- apply(aa, 2, function(x){
      y <- rep("white", P)
      y[x > 0] <- colors_pos[x[x > 0]]
      bb <- x[x < 0]
      y[x < 0] <- colors_neg[128 -abs(bb) +1]
      return(y)
    })
  } else {
    sepcolor <- apply(Z, 2, function(x){
      y <- rep("white", P)
      y[x > 0] <- "red"
      y[x < 0] <- "blue"
      return(y)
    })
  } 
  
  plot(0, 0, ylim=c(-1, (N+1)*sep.height), xlim=c(0, (P+1)*sep.width), type="n",
       axes=F, ylab="", xlab="")
  for(i in 1:N){ 
    rect(xleft = csep, ybottom = rep((N-i+1)*sep.height + sep.height/2, length(csep)),
         xright = csep + sep.width, 
         ytop = rep((N-i)*sep.height + sep.height/2, length(csep)), lty = 1, lwd = 1, 
         col = sepcolor[i, ], border = NA)
  }
  rect(xleft = csep[1], ybottom = N*sep.height + sep.height/2,
       xright = csep[P] + sep.width, 
       ytop =  sep.height/2, lty = 1, lwd = 1,  col = NA, border = "black")
  
}  

