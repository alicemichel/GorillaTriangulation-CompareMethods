library(rootSolve)
library(adehabitatHR)

hyperbolaComput <- function(p1, p2, dt, v=340, xvals= 1:999, yrange=c(-1000, 1000)){
  x1 <- p1[1]
  y1 <- p1[2]
  x2 <- p2[1]
  y2 <- p2[2]
  res <- data.frame(x=numeric(0), y=numeric(0))
  for (x in xvals){
    fun <- function (y) {
      sqrt((x-x1)^2+(y-y1)^2)-sqrt((x-x2)^2+(y-y2)^2)-v*dt
    }
    solution <- uniroot.all(fun, yrange)
    if (length(solution)==0) solution <- NA
    #print(solution)
    res <- rbind(res, data.frame(x=x, y= solution))
  }
  return(res)
}

#########################
#########################

goriLoc <- function(lags1ind, xy, xjump = 10, xextr = 2000, temperature, main=date){
  
  temperature = 25
  v = 331.5 * ((temperature + 273.15)/273.15)^0.5
  
  pams.xy <- xy[rownames(xy) %in% unique(lags1ind$PAM),]
  pams.xy <- pams.xy[complete.cases(pams.xy),]
  pams.xy$PAM <- rownames(pams.xy)
  
  x <- seq(min(pams.xy$lon) - xextr, max(pams.xy$lon) + xextr, xjump)
  y <- seq(min(pams.xy$lat) - xextr, max(pams.xy$lat) + xextr, xjump)
  
  xvals = x
  yrange = c(y[1],y[length(y)])
  
  pamCoords <- merge(pams.xy, lags1ind[,2:3], by="PAM")
  names(pamCoords) <- c("PAM","x","y","timeStamp")
  
  ans <- list()
  for (i in 1:(nrow(pamCoords)-1)){
    for (j in (i+1):nrow(pamCoords)){
      res <- hyperbolaComput(p1= as.numeric(pamCoords[i, c("x", "y")]), p2= as.numeric(pamCoords[j, c("x", "y")]), dt=pamCoords[i, c("timeStamp")]-pamCoords[j, c("timeStamp")], v=v, xvals=xvals, yrange=yrange)
      ans[[length(ans)+1]] <- res			
    }
  }
  ## calculate intersection between any 2 hyperbola
  pamcombos <- data.frame(combn(1:nrow(pamCoords), 3))
  intersecsdf <- data.frame(x=numeric(0), y=numeric(0))
  for (i in 1:ncol(pamcombos)){
    tryCatch(
      expr = {
        pamUsed <- pamcombos[,i]
        
        x1 <- pamCoords[pamUsed[1], "x"]
        x2 <- pamCoords[pamUsed[2], "x"]
        x3 <- pamCoords[pamUsed[3], "x"]
        y1 <- pamCoords[pamUsed[1], "y"]
        y2 <- pamCoords[pamUsed[2], "y"]
        y3 <- pamCoords[pamUsed[3], "y"]
        dt1 <- pamCoords[pamUsed[1], "timeStamp"]- pamCoords[pamUsed[2], "timeStamp"]
        dt2 <- pamCoords[pamUsed[1], "timeStamp"]- pamCoords[pamUsed[3], "timeStamp"]
        equations <- function(x) c(sqrt((x[1]-x1)^2+(x[2]-y1)^2)-sqrt((x[1]-x2)^2+(x[2]-y2)^2)-v*dt1, sqrt((x[1]-x1)^2+(x[2]-y1)^2)-sqrt((x[1]-x3)^2+(x[2]-y3)^2)-v*dt2)
        
        ans1 <- na.omit(ans[[1]])
        ans2 <- na.omit(ans[[2]])
        matDist <- as.matrix(dist(rbind(ans1, ans2)))[(nrow(ans1)+1):(nrow(ans2)+nrow(ans1)), 1:nrow(ans1)]
        indexMin <- which(matDist==min(matDist), arr.ind=TRUE)
        startingVals <- as.numeric(ans1[indexMin[2],])
        
        intersec <- multiroot(f = equations, start = startingVals)$root
        pt <- data.frame(t(intersec))
        names(pt) <- c("x","y")
        pt$pams <- paste0(pamCoords$PAM[pamUsed], collapse = "")
        intersecsdf=rbind(intersecsdf,pt)
      },
      
      error = function(e){
        message(conditionMessage(e))
        NULL
      }
    )
  }
  
  intersecs = intersecsdf
  
  intersecs$x <- ifelse(intersecs$x>max(x) | intersecs$x<min(x) | intersecs$y>max(y) | intersecs$y<min(y),
                        NA,
                        intersecs$x)
  intersecs <- intersecs[complete.cases(intersecs),]
  
  # Format spatial for KDE:
  coordinates(intersecs) <- c("x", "y")
  kud <- kernelUD(intersecs[-1], h="href")
  kud
  kdr <- raster::raster(kud)
  plot(kdr, legend=FALSE, col=cm.colors(8, alpha=0.6), las=1, axes=T, xlim=range(xvals), ylim=yrange, main=paste(main, "- Individual", unique(lags1ind$IndID)))
  points(xy, pch = 20)
  text(xy, labels=rownames(xy), pos=1, cex=.75)
  points(pams.xy, col="blue", cex=2)
  points(intersecsdf, col="red", pch=20)
  idx = raster::which.max(kdr)
  pos = raster::xyFromCell(kdr,idx)
  points(pos, col="blue", pch=8, cex=2, lwd=3)
  legend("bottom", inset = 0.05, pch=8, legend=paste(round(pos,1), collapse=","), box.lwd = 0.5, col="blue")
  
  return(list(intersection = intersecsdf, hyperbola = ans, optimum = pos))
}



