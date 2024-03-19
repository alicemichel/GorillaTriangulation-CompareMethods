
distGorMic <- function(gorilla, mic) {
  xGor = gorilla[,1]
  yGor = gorilla[,2]
  tmp <- vector()
  for (i in 1:nrow(mic)){
    x2 = mic[i,1]
    y2 = mic[i,2]
    tmp <- c(tmp,sqrt( (xGor-x2)^2 + (yGor-y2)^2))
  }
  return(tmp)
}


pamNeighbors <- function(pam, xy, n=10, self=TRUE){
  pam.pt <- xy[pam,]
  xy$d <- distGorMic(pam.pt, xy)
  s <- ifelse(self==TRUE, 1, 2)
  pamNeighborIDs = rownames(xy[order(xy$d),][s:(n+1),])
  return(pamNeighborIDs)
}


