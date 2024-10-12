transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

timeline <- function(dets, pal = RColorBrewer::brewer.pal(5, "Set2"), pdf=TRUE){
  
  start.date = substr(dets$sound.files, start = 4, stop = 11)[1]
  
  full.starts.all <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = dets$sound.files, clip.start = dets$startClip.GPS)
  full.starts.all <- as.POSIXct(ifelse(full.starts.all>"2023-01-17 19:00:00 PST", full.starts.all, full.starts.all+86400))
  full.starts.all.df <- data.frame("ind"=dets$indMoved, "start"=full.starts.all)

  num_colors <- nlevels(as.factor(dets$indMoved))
  colors <- pal[as.factor(dets$indMoved)]
  
  cb1 <- min(full.starts.all)
  cblast <- max(full.starts.all)
  
  par(mfrow=c(2,1))
  if (pdf==TRUE){
    pdf(paste0(start.date,"timeline", paste(unique(dets$indMoved), collapse=""),".pdf"), width=12, height=3)
  }
  plot(full.starts.all.df$start, rep(1, nrow(full.starts.all.df)), axes=F, pch="|", xlab=NA, ylab=NA, cex=3, main=as.Date(full.starts.all[[1]][1]-86400), xlim=c(cb1-5, cblast+5),
       col=colors)
  #plot(x=ind1$full.start, y=rep(1, length(ind1$full.start)), axes=F, pch="|", xlab=NA, ylab=NA, cex=3, main=as.Date(ind1$full.start[1]-86400), xlim=c(cb1-5, cblast+5))
  axis.POSIXct(side = 1, tick=F, pos=0.8)
  abline(h=1)
  #points(x=ind2$full.start, y=rep(1, length(ind2$full.start)), pch="|", col="red", cex=3)
  mtext(text = "Chest beats from individuals", cex=1)
  legend("top", inset=-0.2, legend = unique(dets$indMoved), text.col = unique(colors), bty="n", horiz=TRUE, cex = 1, xpd=T)
  #mtext(text = paste(unique(substr(names(namereps), 1, 1)), collapse=", "), line = -1)
  #mtext(text = "and", line = -1)
  #mtext(text = paste("          ", unique(ind1$ind)), line = -1)
  if (pdf==TRUE){dev.off()}
  # Get overall first and last CB of all individuals that night:
}
  
  # slide CBs of individual who lasts less long within this time range:
  

# ind1 <- night[[1]]$detections
# ind2 <- night[[2]]$detections
# ind1$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind1$sound.files, clip.start = ind1$startClip.GPS)
# ind2$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind2$sound.files, clip.start = ind2$startClip.GPS)
# 

gapsSim <- function(dets,short,pal = RColorBrewer::brewer.pal(5, "Set2"), pdf=TRUE){
  
  num_colors <- nlevels(as.factor(dets$indMoved))
  dets$colors <- pal[as.factor(dets$indMoved)]
  
  start.date = substr(dets$sound.files, start = 4, stop = 11)[1]
  
  full.starts.all <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = dets$sound.files, clip.start = dets$startClip.GPS)
  full.starts.all <- as.POSIXct(ifelse(full.starts.all>"2023-01-17 19:00:00 PST", full.starts.all, full.starts.all+86400))
  full.starts.all.df <- data.frame("ind"=dets$indMoved, "start"=full.starts.all)
  
  # for all individuals, the chest beat period that night:
  cb1 <- min(full.starts.all)
  cblast <- max(full.starts.all)
  
  ind2 <- dets[dets$indMoved==short,]
  ind2$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind2$sound.files, clip.start = ind2$startClip.GPS)
  colori <- ind2$colors[1]
  
  others <- unique(dets$indMoved)[unique(dets$indMoved)!=short]
  
  num <- data.frame()
  for (long in others){
    ind1 <- dets[dets$indMoved==long,]
    ind1$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind1$sound.files, clip.start = ind1$startClip.GPS)
  
    minIIgaps <- min(abs(ind2[1,]$full.start - ind1$full.start))
    indTgaps <- NA
    indTfrom1 <- 0
    for (i in 2:nrow(ind2)){
      minIIgaps <- c(minIIgaps, (ind2[i,]$full.start - ind1$full.start)[abs(ind2[i,]$full.start - ind1$full.start)==min(abs(ind2[i,]$full.start - ind1$full.start))])
      indTgaps <- c(indTgaps, ind2[i,]$full.start - ind2[i-1,]$full.start)
      indTfrom1 <- c(indTfrom1, ind2[i,]$full.start - ind2[1,]$full.start)
    }
    
    # move by 1 second over the whole thing, starting from 0, ending when the first one is at the last end time and the rest loop around
    seclength <- round(as.numeric(cblast - cb1) * 60) #seconds
    starts <- list()
    gaps <- data.frame()
    for (i in 1:seclength){ #no original! plot separately (in red)
      tmp <- cb1 + 1*i + indTfrom1*60
      tmp2 <- c()
      for (j in 1:length(tmp)){
        if (tmp[j]>cblast){
          tmp[j] <- tmp[j]-cblast + cb1 #amount over, add to cb1 to get new start
        }
        tmp2 <- c(tmp2, tmp[j] - ind1$full.start[abs(tmp[j] - ind1$full.start)==min(abs(tmp[j] - ind1$full.start))])
      }
      starts[[i]] <- tmp
      gaps <- rbind(gaps, tmp2[-1]) #the first always = the last for ind1 since he's the longest
    }
    names(gaps) <- paste0("CB",2:(length(gaps)+1))
    abs.gaps.sim <- abs(as.vector(unlist(gaps)))
    
    #simulate log-normal distribution with that mean and sd:
    par(mfrow=c(1,1))
    plot(density(lnorm.dens <- rlnorm(n=3000,
                                      meanlog=mean(log(abs.gaps.sim)),
                                      sdlog=sd(log(abs.gaps.sim)))), bty="l", main="Log normal",las=1, xlim=c(0,50))
    #simulate Gamma distribution with that shape and error:
    gamma <- MASS::fitdistr(abs.gaps.sim, "Gamma")
    plot(density(gamma.dens <- rgamma(n=3000,
                                      shape=gamma$estimate[1],
                                      rate=gamma$estimate[2])), bty="l", main="Gamma",las=1, xlim=c(0,50))
    #calculate confidence intervals for each:
    ci95gamma <- quantile(gamma.dens, probs = c(0.05,.95))
    ci95lnorm <- quantile(lnorm.dens, probs = c(0.05,.95))
    
    #plot the simulated data:
    title=paste(unique(ind2$indMoved), "shifted over", unique(ind1$indMoved)); setx="Minimum time between CBs of different individuals (minutes)"
    hist(abs.gaps.sim, bty="l", las=1, breaks=20, #breaks = quantile(abs.gaps.sim, probs = c(0,0.01,seq(0.1,0.9,0.1),1)),
         xlab=setx, main=title)
    #or plot the gamma distribution:
    if (pdf==TRUE){
      pdf(paste0(start.date,"distributionSimGaps", unique(ind1$indMoved), unique(ind2$indMoved),".pdf"), width=8, height = 6)
    }
    plot(density(gamma.dens <- rgamma(n=3000,shape=gamma$estimate[1],rate=gamma$estimate[2]), from = 0), yaxs="i", xaxs="i",
         bty="l", xlab=setx, main=title,las=1, zero.line=FALSE, ylab="Gamma density of simulated values")
    sign=ifelse(ci95gamma[1]>min(abs.gaps.sim),-0.001,0.001)
    #shade lower 99% CI gamma:
    abline(v=c(ci95gamma[1]-seq(ci95gamma[1],min(abs.gaps.sim),sign)), col=transp("lightblue",0.01),lwd=2)
    #shade lower 99% CI log-normal:
    # abline(v=c(ci95lnorm[1]-seq(ci95lnorm[1],min(abs.gaps.sim),sign)), col=transp("lightgreen",0.006),lwd=2)
    #add lines for obs data:
    abs.min.gaps.obs <- abs(as.numeric(minIIgaps)/60)
    
    points(x=abs.min.gaps.obs, y=rep(0, length(minIIgaps)), col=colori, pch=20, xpd=T)
    
    # what fraction are significantly replies:
    length(which(abs.min.gaps.obs<ci95lnorm[1]))/length(abs.min.gaps.obs)
    length(which(abs.min.gaps.obs<ci95gamma[1]))/length(abs.min.gaps.obs)
    
    text(ci95gamma[1]+2, 0.005, labels=paste0("n=", length(which(abs.min.gaps.obs<ci95gamma[1])), "/",length(abs.min.gaps.obs)," of ", unique(ind2$indMoved), "'s chest beats are closer to ", unique(ind1$indMoved), "'s (n=", nrow(ind1), ") than \n<5% of the simulated inter-individual time gaps"), col=colori, adj=0)
    if (pdf==TRUE){dev.off()}
    
    
    which(abs.min.gaps.obs>=ci95gamma[1])
    
    tmp <- c(short, long, length(which(abs.min.gaps.obs<ci95gamma[1])), length(which(abs.min.gaps.obs<ci95lnorm[1])), nrow(ind2), paste0(which(abs.min.gaps.obs<ci95gamma[1]), collapse="_"), paste0(which(abs.min.gaps.obs<ci95lnorm[1]), collapse="_"))
    num <- rbind(num, tmp)
  }
  names(num) <- c("FocalInd", "BkgdInd", "num_<5%_gamma", "num_<5%_lnorm", "num_cbs_total", "ids_of_closer_gamma", "ids_of_closer_lnorm")
  
  
  num$reply_to_any_gamma <- length(unique(as.numeric(unlist(strsplit(num$ids_of_closer_gamma, "_")))))
  num$reply_to_any_lnorm <- length(unique(as.numeric(unlist(strsplit(num$ids_of_closer_lnorm, "_")))))
  
  return(num)
}




simGapsListCBs <- function(night1,short,long,sign="neg"){
  
  night <- night1[c(names(night1)[long], names(night1)[short], (names(night1)[!(names(night1) %in% c(names(night1)[short], names(night1)[long]))]))]
  
  sign=ifelse(sign=="neg",-0.001,0.001)
  
  start.date = substr(night[[1]]$detections$sound.files, start = 4, stop = 11)[1]
  
  full.starts.all <- lapply(night, FUN = function(x) as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = x$detections$sound.files, clip.start = x$detections$startClip.GPS))
  full.starts.all <- lapply(full.starts.all, function(x) as.POSIXct(ifelse(x>"2023-01-17 19:00:00 PST", x, x+86400)))
  
  cb1 <- min(do.call(c, lapply(full.starts.all, min)))
  cblast <- max(do.call(c, lapply(full.starts.all, max)))
  
  ind1 <- night[[1]]$detections
  ind2 <- night[[2]]$detections
  ind1$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind1$sound.files, clip.start = ind1$startClip.GPS)
  ind2$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind2$sound.files, clip.start = ind2$startClip.GPS)
  
  # Plot the chest beat timeline colored by individual:
  
  par(mfrow=c(4,1))
  #pdf(paste0(start.date,"timeline", paste(names(night), collapse=""),".pdf"), width=10, height=3)
  namereps=do.call(c, lapply(full.starts.all, function(x) rep(1, length(x))))
  cols <- c("red","blue","turquoise","goldenrod","violet")[as.numeric(as.factor(substr(names(namereps), 1, 1)))]
  #RColorBrewer::brewer.pal(n = 5, name = "Set1")[as.numeric(as.factor(substr(names(namereps), 1, 1)))]
  #nationalparkcolors::park_palette("Everglades",5)[as.numeric(as.factor(substr(names(namereps), 1, 1)))]
  plot(do.call(c,full.starts.all), namereps, axes=F, pch="|", xlab=NA, ylab=NA, cex=3, main=as.Date(full.starts.all[[1]][1]-86400), xlim=c(cb1-5, cblast+5),
       col=cols)
  #plot(x=ind1$full.start, y=rep(1, length(ind1$full.start)), axes=F, pch="|", xlab=NA, ylab=NA, cex=3, main=as.Date(ind1$full.start[1]-86400), xlim=c(cb1-5, cblast+5))
  axis.POSIXct(side = 1, tick=F, pos=0.8)
  abline(h=1)
  #points(x=ind2$full.start, y=rep(1, length(ind2$full.start)), pch="|", col="red", cex=3)
  mtext(text = "Chest beats from individuals", cex=1)
  legend("top", inset=-0.1, legend = unique(substr(names(namereps), 1, 1)), text.col = unique(cols), bty="n", horiz=TRUE, cex = 1.5)
  #mtext(text = paste(unique(substr(names(namereps), 1, 1)), collapse=", "), line = -1)
  #mtext(text = "and", line = -1)
  #mtext(text = paste("          ", unique(ind1$ind)), line = -1)
  #dev.off()
  # Get overall first and last CB of all individuals that night:
  
  
  # slide CBs of individual who lasts less long within this time range:
  
  # for each individual:
  cb1 <- min(c(min(ind1$full.start), min(ind2$full.start)))
  cblast <- max(c(max(ind1$full.start), max(ind2$full.start)))
  
  minIIgaps <- min(abs(ind2[1,]$full.start - ind1$full.start))
  indTgaps <- NA
  indTfrom1 <- 0
  for (i in 2:nrow(ind2)){
    minIIgaps <- c(minIIgaps, (ind2[i,]$full.start - ind1$full.start)[abs(ind2[i,]$full.start - ind1$full.start)==min(abs(ind2[i,]$full.start - ind1$full.start))])
    indTgaps <- c(indTgaps, ind2[i,]$full.start - ind2[i-1,]$full.start)
    indTfrom1 <- c(indTfrom1, ind2[i,]$full.start - ind2[1,]$full.start)
  }
  (mean.obs.gaps <- mean(abs.min.gaps.obs))
  
  # move by 1 second over the whole thing, starting from 0, ending when the first one is at the last end time and the rest loop around
  seclength <- round(as.numeric(cblast - cb1) * 60) #seconds
  starts <- list()
  gaps <- data.frame()
  for (i in 1:seclength){ #no original! plot separately (in red)
    tmp <- cb1 + 1*i + indTfrom1*60
    tmp2 <- c()
    for (j in 1:length(tmp)){
      if (tmp[j]>cblast){
        tmp[j] <- tmp[j]-cblast + cb1 #amount over, add to cb1 to get new start
      }
      tmp2 <- c(tmp2, tmp[j] - ind1$full.start[abs(tmp[j] - ind1$full.start)==min(abs(tmp[j] - ind1$full.start))])
    }
    starts[[i]] <- tmp
    gaps <- rbind(gaps, tmp2[-1]) #the first always = the last for ind1 since he's the longest
  }
  names(gaps) <- paste0("CB",2:(length(gaps)+1))
  abs.gaps.sim <- abs(as.vector(unlist(gaps)))
  mean.sim.gaps <- sort(rowMeans(abs(gaps)))
  
  par(mfrow=c(1,1))
  hist(mean.sim.gaps, xlim=c(max(0, min(mean.obs.gaps, mean.sim.gaps)), max(c(mean.obs.gaps, mean.sim.gaps))), main=paste(unique(ind2$ind), "and", unique(ind1$ind)), xlab="Simulated inter-individual nearest chest beat time gaps (secs)", color="grey95", border="grey70")
  abline(v=mean.obs.gaps, lty=4, col="red", lwd=2)
  
  pval <- rank(c(mean.obs.gaps, mean.sim.gaps))[1]/length(unique(rank(c(mean.obs.gaps, mean.sim.gaps))))
  
  tmp <- c(short, long, nrow(ind2), nrow(ind1), mean.obs.gaps, seclength, pval)
  num <- rbind(num, tmp)
  names(num) <- c("FocalInd", "BkgdInd", "num_cbs_focal", "num_cbs_bkgd", "mean_gap", "num_sims", "pval")
  return(num)
}



#   
#   
#   #simulate log-normal distribution with that mean and sd:
#   par(mfrow=c(1,1))
#   plot(density(lnorm.dens <- rlnorm(n=3000,
#                                     meanlog=mean(log(abs.gaps.sim)),
#                                     sdlog=sd(log(abs.gaps.sim)))), bty="l", main="Log normal",las=1, xlim=c(0,50))
#   #simulate Gamma distribution with that shape and error:
#   gamma <- MASS::fitdistr(abs.gaps.sim, "Gamma")
#   plot(density(gamma.dens <- rgamma(n=3000,
#                                     shape=gamma$estimate[1],
#                                     rate=gamma$estimate[2])), bty="l", main="Gamma",las=1, xlim=c(0,50))
#   #calculate confidence intervals for each:
#   ci95gamma <- quantile(gamma.dens, probs = c(0.05,.95))
#   ci95lnorm <- quantile(lnorm.dens, probs = c(0.05,.95))
#   
#   #plot the simulated data:
#   title=paste(unique(ind2$ind), "shifted over", unique(ind1$ind)); setx="Minimum time between CBs of different individuals (minutes)"
#   hist(abs.gaps.sim, bty="l", las=1, breaks=20, #breaks = quantile(abs.gaps.sim, probs = c(0,0.01,seq(0.1,0.9,0.1),1)),
#        xlab=setx, main=title)
#   #or plot the gamma distribution:
#   #pdf(paste0(start.date,"distributionSimGaps", unique(ind1$ind), unique(ind2$ind),".pdf"), width=8, height = 6)
#   plot(density(gamma.dens <- rgamma(n=3000,shape=gamma$estimate[1],rate=gamma$estimate[2]), from = 0), yaxs="i", xaxs="i",
#        bty="l", xlab=setx, main=title,las=1, xlim=c(0,50), zero.line=FALSE, ylab="Gamma density of simulated values")
#   #shade lower 99% CI gamma:
#   abline(v=c(ci95gamma[1]-seq(ci95gamma[1],min(abs.gaps.sim),sign)), col=transp("lightblue",0.01),lwd=2)
#   #shade lower 99% CI log-normal:
#   # abline(v=c(ci95lnorm[1]-seq(ci95lnorm[1],min(abs.gaps.sim),-0.001)), col=transp("lightgreen",0.006),lwd=2)
#   #add lines for obs data:
#   abs.min.gaps.obs <- abs(as.numeric(minIIgaps)/60)
#   points(x=abs.min.gaps.obs, y=rep(0, length(minIIgaps)), col="red", pch=20, xpd=T)
#   
#   # what fraction are significantly replies:
#   length(which(abs.min.gaps.obs<ci95lnorm[1]))/length(abs.min.gaps.obs)
#   length(which(abs.min.gaps.obs<ci95gamma[1]))/length(abs.min.gaps.obs)
#   
#   text(ci95gamma[1]+2, 0.005, labels=paste0("n=", length(which(abs.min.gaps.obs<ci95gamma[1])), "/",length(abs.min.gaps.obs)," of ", unique(ind2$ind), "'s chest beats are closer to ", unique(ind1$ind), "'s than \n<5% of the simulated inter-individual time gaps"), col="red", adj=0)
#   
#   #dev.off()
# }
# 
# 
# 
# #  minIIgapsSim <- lapply(starts, function(x) min(abs(x[1] - ind2$full.start)))
# # for (i in 2:nrow(ind1)){
# #   minIIgapsSim <- lapply(starts, function(x) min(abs(x[i] - ind2$full.start)))
# # }
# # 
# # # Check for within-individual escalation:
# # 
# # plot(1:length(indTgaps), indTgaps, las=1, xlab="Chest beat number", ylab="Time gap", bty="l", pch=20)
# # 
# # # Check for btwn-individual escalation:
# # 
# # plot(1:length(minIIgaps), minIIgaps, las=1, xlab="Chest beat number", ylab="Time gap", bty="l", pch=20)
# # 
# # 
# # 
# # plotConseq(rbind(night[[1]]$detections,night[[2]]$detections))
# # 
# # # but make it just a timeline
# # 
# # # compute the metric: nearest-neighbor inter-individual gaps (doesn't matter which individual? or just pairwise?)
# # 
# # # then find the start time and the end time of all chest beating
# # 
# # # then slide each around by 1 second randomly, if it gets to the end, loop it back to the beginning so it's not ordered which is weird for escalation but makes sense for this simulation
# # 
# # # compute the metric
# # 
# # # do that 1000 times
# # 
# # # plot them all, color the real ones - are they in the tail?
