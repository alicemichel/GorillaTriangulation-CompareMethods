

cutNbuff <- function(start, cut, buffer=1){
  
  newst <- rep(NA, length(start))
  for (i in 1:length(start)){
    newst[i] <- start[i] + cut[i] - buffer
  }
  return(newst)
}


plotConseq <- function(dets.ordered, start0=TRUE, colors=NA, pal=palette()){
  
  opar <- par() 
  
  dets <- dets.ordered
  
  start.date <- substr(dets$sound.files, start = 4, stop = 11)
  inds=sort(unique(dets$ind))
  
  # plot the consecutive chest beats of each individual over one another
  dets$fullTimeFromClipStart <- as.POSIXct(start.date, format = "%Y%m%d", tz = "GMT") + fullTimeFromClipStart(dets$sound.files,dets$startClip.GPS)
  
  dets <- dets[order(dets$fullTimeFromClipStart),]
  
  tb <- data.frame(table(dets$ind, dets$fullTimeFromClipStart))
  ind1=as.character(tb[1,]$Var1)
  tb1 <- tb[tb$Var1==ind1 & tb$Freq==1,]
  
  par(mfrow=c(2,1), mar=c(0,5,2,1), oma=c(4,0,0,0))
  
  plot(dets[dets$ind==ind1,]$fullTimeFromClipStart, 1:nrow(dets[dets$ind==ind1,]), bty="l", pch=20, col=pal[1], las=1, xaxt="n", xlab=NA, ylab="Cumulative # chest beats\nper individual", main=paste("Individuals chest beating on", start.date[1]), xlim=c(min(dets$fullTimeFromClipStart)-1, max(dets$fullTimeFromClipStart)+1), ylim=c(0,max(table(dets$ind))))
  lines(dets[dets$ind==ind1,]$fullTimeFromClipStart, 1:nrow(dets[dets$ind==ind1,]), col=pal[1])
  text(dets[dets$ind==ind1,]$fullTimeFromClipStart[1], 0, labels=ind1, cex=0.5, col=pal[1])
  
  
  dets$colvec <- ifelse(dets$ind==ind1,1,NA)
  
  otherinds <- inds[inds!=ind1]
  
  #count number of ind1 before the next ind for every next ind
  for (other in otherinds){
    
    dets[which(dets$ind %in% other),]$colvec <- which(otherinds %in% other)+1
    
    if (length(colors)>1){
      dets$colvec <- colors
    }
    
    tbi <- tb[tb$Var1==other & tb$Freq==1,]
    first1 <- as.POSIXct(tbi[1,]$Var2, "%Y-%m-%d %H:%M:%S")
    tb1before <- tb1[as.POSIXct(tb1$Var2, "%Y-%m-%d %H:%M:%S")<first1,]
    cntb4 <- nrow(tb1before)
    
    if (start0==TRUE){
      cntb4=1
    }
    
    points(dets[dets$ind==other,]$fullTimeFromClipStart, cntb4:(nrow(dets[dets$ind==other,])+cntb4-1), pch=20, col=pal[dets[dets$ind==other,]$colvec])
    text(dets[dets$ind==other,]$fullTimeFromClipStart[1], 0, labels=other, col=pal[dets[dets$ind==other,]$colvec], cex=0.5)
    lines(dets[dets$ind==other,]$fullTimeFromClipStart, cntb4:(nrow(dets[dets$ind==other,])+cntb4-1), lwd=1, col=pal[dets[dets$ind==other,]$colvec])
  }
  legend("topleft", cex = 1.1, inset = 0.05, title.adj = 0, legend=sort(unique(dets$indMoved)), col=pal[aggregate(colvec ~ indMoved, data=dets, FUN=unique)$colvec ], pch=15, bty="n", title = "Individuals")
  
  
  plot(dets$fullTimeFromClipStart, 1:nrow(dets), col = pal[dets$colvec], bty="l", pch=15, cex=1.5, las=1, xlab=NA, ylab="Cumulative # chest beats\noverall", main=NA, cex.axis=1)
  mtext(side = 1, line = 3, adj = 0.5, text="Time")
  
  par(opar)      
  
}

