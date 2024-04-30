
setwd("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/")

source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_distance_calcs.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_get_pull_times.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_clock_drift.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_raven_approx_ord.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_cutNbuff_plotConseq.R")
transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

## Consecutive delays analysis

# read in CBs from each night that's done
feb6 <- readRDS("2024.04.07_20230206_individuals_J_U/20240408_161350_CBs.rds")
dec10 <- readRDS("2024.04.15_20221210_individuals_E_J_V_Mdef_Mmb_Dpok_Dsm_Vclap/20240414_132623_CBs.rds") #may not have all of V
jan16 <- readRDS("2024.04.21_20230116_individuals_B1_B2_D1_D2_V/20240421_213648_CBs.rds")
jan17 <- readRDS("2024.04.22_20230117_individuals_D17_D18_K17_K18_V17_V18/20240422_215708_CBs.rds")

# compute consecutive gaps:
# should be based on one PAM...or both their closest so even bias? or which one? ************** FOR LATER

# keep each individual separate, calculate closest in other individual df...

setwd("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/2024.04.07_20230206_individuals_J_U//")
names(jan17)
night1 = list()
night1[["D"]]$detections = rbind(jan17[[1]]$detections,jan17[[2]]$detections)
night1[["D"]]$detections$ind="D"
night1[["K"]]$detections = rbind(jan17[[3]]$detections,jan17[[4]]$detections)
night1[["K"]]$detections$ind="K"
night1[["V"]]$detections = rbind(jan17[[5]]$detections,jan17[[6]]$detections)
night1[["V"]]$detections$ind="V"
night1 <- dec10
# how many CBs per individual?
lapply(night1, function(x) nrow(x$detections))
# put the shorter one as individual 2 and run the function:
simGaps(night1,short=1,long=2,sign="neg")

simGaps <- function(night1,short,long,sign="neg"){
  
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
  pdf(paste0(start.date,"timeline", paste(names(night), collapse=""),".pdf"), width=10, height=3)
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
  dev.off()
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
  title=paste(unique(ind2$ind), "shifted over", unique(ind1$ind)); setx="Minimum time between CBs of different individuals (minutes)"
  hist(abs.gaps.sim, bty="l", las=1, breaks=20, #breaks = quantile(abs.gaps.sim, probs = c(0,0.01,seq(0.1,0.9,0.1),1)),
       xlab=setx, main=title)
  #or plot the gamma distribution:
  pdf(paste0(start.date,"distributionSimGaps", unique(ind1$ind), unique(ind2$ind),".pdf"), width=8, height = 6)
  plot(density(gamma.dens <- rgamma(n=3000,shape=gamma$estimate[1],rate=gamma$estimate[2]), from = 0), yaxs="i", xaxs="i",
       bty="l", xlab=setx, main=title,las=1, xlim=c(0,50), zero.line=FALSE, ylab="Gamma density of simulated values")
  #shade lower 99% CI gamma:
  abline(v=c(ci95gamma[1]-seq(ci95gamma[1],min(abs.gaps.sim),sign)), col=transp("lightblue",0.01),lwd=2)
  #shade lower 99% CI log-normal:
  # abline(v=c(ci95lnorm[1]-seq(ci95lnorm[1],min(abs.gaps.sim),-0.001)), col=transp("lightgreen",0.006),lwd=2)
  #add lines for obs data:
  abs.min.gaps.obs <- abs(as.numeric(minIIgaps)/60)
  points(x=abs.min.gaps.obs, y=rep(0, length(minIIgaps)), col="red", pch=20, xpd=T)
  
  # what fraction are significantly replies:
  length(which(abs.min.gaps.obs<ci95lnorm[1]))/length(abs.min.gaps.obs)
  length(which(abs.min.gaps.obs<ci95gamma[1]))/length(abs.min.gaps.obs)
  
  text(ci95gamma[1]+2, 0.005, labels=paste0("n=", length(which(abs.min.gaps.obs<ci95gamma[1])), "/",length(abs.min.gaps.obs)," of ", unique(ind2$ind), "'s chest beats are closer to ", unique(ind1$ind), "'s than \n<5% of the simulated inter-individual time gaps"), col="red", adj=0)
  
  dev.off()
}



#  minIIgapsSim <- lapply(starts, function(x) min(abs(x[1] - ind2$full.start)))
# for (i in 2:nrow(ind1)){
#   minIIgapsSim <- lapply(starts, function(x) min(abs(x[i] - ind2$full.start)))
# }
# 
# # Check for within-individual escalation:
# 
# plot(1:length(indTgaps), indTgaps, las=1, xlab="Chest beat number", ylab="Time gap", bty="l", pch=20)
# 
# # Check for btwn-individual escalation:
# 
# plot(1:length(minIIgaps), minIIgaps, las=1, xlab="Chest beat number", ylab="Time gap", bty="l", pch=20)
# 
# 
# 
# plotConseq(rbind(night[[1]]$detections,night[[2]]$detections))
# 
# # but make it just a timeline
# 
# # compute the metric: nearest-neighbor inter-individual gaps (doesn't matter which individual? or just pairwise?)
# 
# # then find the start time and the end time of all chest beating
# 
# # then slide each around by 1 second randomly, if it gets to the end, loop it back to the beginning so it's not ordered which is weird for escalation but makes sense for this simulation
# 
# # compute the metric
# 
# # do that 1000 times
# 
# # plot them all, color the real ones - are they in the tail?
