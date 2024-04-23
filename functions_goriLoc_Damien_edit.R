# load libraries and functions
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

goriLoc <- function(lagdf.IndID.PAM.lag, mic.coords.rownames.pams, xjump = 10, xextr = 500, temperature, main=date){
  
  v = 331.5 * ((temperature + 273.15)/273.15)^0.5
  
  pams.mic.coords.rownames.pams <- mic.coords.rownames.pams[rownames(mic.coords.rownames.pams) %in% unique(lagdf.IndID.PAM.lag$PAM),]
  pams.mic.coords.rownames.pams <- pams.mic.coords.rownames.pams[complete.cases(pams.mic.coords.rownames.pams),]
  pams.mic.coords.rownames.pams$PAM <- rownames(pams.mic.coords.rownames.pams)
  
  x <- seq(min(pams.mic.coords.rownames.pams$lon) - xextr, max(pams.mic.coords.rownames.pams$lon) + xextr, xjump)
  y <- seq(min(pams.mic.coords.rownames.pams$lat) - xextr, max(pams.mic.coords.rownames.pams$lat) + xextr, xjump)
  
  xvals = x
  yrange = c(y[1],y[length(y)])
  
  pamCoords <- merge(pams.mic.coords.rownames.pams, lagdf.IndID.PAM.lag[,2:3], by="PAM")
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
  
  # set up plot:
  plot(x = range(xvals), y = yrange, type="n", main=paste(main, "- Individual", unique(lagdf.IndID.PAM.lag$IndID)), 
       las=1, xlab=NA, ylab=NA)
  
  # Format spatial for KDE:
  coordinates(intersecs) <- c("x", "y")
  pos <- NA
  if (length(intersecs)>=5){
    kud <- kernelUD(intersecs[0], h="href")
    kdr <- raster::raster(kud)
    alPal <- colorRampPalette(c('white','grey80'))
    plot(kdr, legend=FALSE, col=alPal(8), add=T)
    idx = raster::which.max(kdr)
    pos = raster::xyFromCell(kdr,idx)
    legend("bottom", inset = 0.05, pch=8, legend=paste(round(pos,1), collapse=","), box.lwd = 0.5, col="blue")
  }
  points(mic.coords.rownames.pams, pch = 20)
  text(mic.coords.rownames.pams, labels=rownames(mic.coords.rownames.pams), pos=1, cex=.75)
  points(pams.mic.coords.rownames.pams, col="blue", cex=2)
  
  keyPAM <- lagdf.IndID.PAM.lag[!is.na(lagdf.IndID.PAM.lag$lag) & lagdf.IndID.PAM.lag$lag==0,]$PAM
  intersecsdf$pam1x <- intersecsdf$pam1y <- intersecsdf$pam2x <- intersecsdf$pam2y <- intersecsdf$pam3x <- intersecsdf$pam3y <- NA
  for (i in 1:nrow(intersecsdf)){
    intersecsdf[i,]$pam1x <- pamCoords[pamCoords$PAM == substr(intersecsdf[i,]$pams, 1, 1),]$x
    intersecsdf[i,]$pam1y <- pamCoords[pamCoords$PAM == substr(intersecsdf[i,]$pams, 1, 1),]$y
    intersecsdf[i,]$pam2x <- pamCoords[pamCoords$PAM == substr(intersecsdf[i,]$pams, 2, 2),]$x
    intersecsdf[i,]$pam2y <- pamCoords[pamCoords$PAM == substr(intersecsdf[i,]$pams, 2, 2),]$y
    intersecsdf[i,]$pam3x <- pamCoords[pamCoords$PAM == substr(intersecsdf[i,]$pams, 3, 3),]$x
    intersecsdf[i,]$pam3y <- pamCoords[pamCoords$PAM == substr(intersecsdf[i,]$pams, 3, 3),]$y
  }
  intersecsdf$pamsCntrX <- (intersecsdf$pam1x + intersecsdf$pam2x + intersecsdf$pam3x)/3
  intersecsdf$pamsCntrY <- (intersecsdf$pam1y + intersecsdf$pam2y + intersecsdf$pam3y)/3
  
  intersecsdf$dist2key <- distGorMic(pamCoords[pamCoords$PAM==keyPAM,2:3], data.frame(x=intersecsdf$pamsCntrX,y=intersecsdf$pamsCntrY))
  alPal2 <- colorRampPalette(c('red','grey80'))
  intersecsdf$color <- alPal2(10)[as.numeric(cut(intersecsdf$dist2key,breaks = 10))]
  
  points(intersecsdf[,c("x","y")], pch=15, cex=1.9, col=intersecsdf$color)
  text(intersecsdf[,c("x","y")], labels=intersecsdf$pams, cex=0.4, col="white")
  
  points(pos, col="blue", pch=8, cex=2, lwd=3)
  
  
  return(list(intersection = intersecsdf, hyperbola = ans, optimum = pos))
}



