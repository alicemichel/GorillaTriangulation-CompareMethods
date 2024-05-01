source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_get_pull_times.R")
feb6 <- readRDS("2024.04.07_20230206_individuals_J_U/20240408_161350_CBs.rds")
night1 = feb6
# how many CBs per individual?
lapply(night1, function(x) nrow(x$detections))
short=2
long=1


night <- night1[c(names(night1)[long], names(night1)[short], (names(night1)[!(names(night1) %in% c(names(night1)[short], names(night1)[long]))]))]


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

par(mfrow=c(2,1))
#pdf(paste0(start.date,"timeline", paste(names(night), collapse=""),".pdf"), width=10, height=3)
namereps=do.call(c, lapply(full.starts.all, function(x) rep(1, length(x))))
cols <- c("red","blue","turquoise","goldenrod","violet")[as.numeric(as.factor(substr(names(namereps), 1, 1)))]
plot(do.call(c,full.starts.all), namereps, axes=F, pch="|", xlab=NA, ylab=NA, cex=3, main=as.Date(full.starts.all[[1]][1]-86400), xlim=c(cb1-5, cblast+5),
     col=cols)
axis.POSIXct(side = 1, tick=F, pos=0.8)
abline(h=1)
mtext(text = "Chest beats from individuals", cex=1)
legend("top", inset=-0.1, legend = unique(substr(names(namereps), 1, 1)), text.col = unique(cols), bty="n", horiz=TRUE, cex = 1.5)


# slide CBs of individual who lasts less long within this time range:

# for each individual:
cb1 <- min(c(min(ind1$full.start), min(ind2$full.start)))
cblast <- max(c(max(ind1$full.start), max(ind2$full.start)))

minIIgaps <- as.numeric(min(ind2[1,]$full.start - ind1$full.start[(ind2[1,]$full.start - ind1$full.start) > 0]), "secs")
indTgaps <- NA
indTfrom1 <- 0
for (i in 2:nrow(ind2)){
  gapi <- as.numeric(ind2[i,]$full.start - ind1$full.start, "secs")
  minIIgaps <- c(minIIgaps, min(gapi[gapi>0])) #abs(gapi)
  indTgaps <- c(indTgaps, as.numeric(ind2[i,]$full.start - ind2[i-1,]$full.start, "secs"))
  indTfrom1 <- c(indTfrom1, as.numeric(ind2[i,]$full.start - ind2[1,]$full.start, "secs"))
}
minIIgaps
(mean.obs.gaps <- mean(minIIgaps))#[-length(minIIgaps)])) #secs


## is there a good argument to remove the last one in case the other guy has given up??

# ford <- rbind(ind1[,c("ind", "full.start")], ind2[,c("ind", "full.start")])
# ford <- ford[order(ford$full.start),]
# saveRDS(ford, "20230206_UJ_times.rds")

# move by 1 second over the whole thing, starting from 0, ending when the first one is at the last end time and the rest loop around
# for the focal individual:
cb1.1 <- ind2$full.start[1]
seclength <- round(as.numeric(cblast - cb1, "secs"))
#starts <- list()
gaps <- data.frame()
for (i in 0:seclength){ #includes original!
  tmp <- as.numeric(cb1.1, "secs") + 1*i + indTfrom1
  tmp2 <- c()
  for (j in 1:length(tmp)){
    if (tmp[j]>as.numeric(cblast, "secs")){
      tmp[j] <- tmp[j]-as.numeric(cblast, "secs") + as.numeric(cb1, "secs") #amount over, add to cb1 to get new start
    }
    gapj <- tmp[j] - as.numeric(ind1$full.start, "secs")
    tmp2 <- c(tmp2, min(gapj[gapj>0])) #abs(gapj)
  }
  #starts[[i]] <- tmp
  gaps <- rbind(gaps, tmp2) #the first always = the last for ind1 since he's the longest
}
names(gaps) <- paste0("CB",2:(length(gaps)+1))
#abs.gaps.sim <- abs(as.vector(unlist(gaps)))
(mean.sim.gaps <- apply(gaps, 1, FUN= function(x) mean(x)))#[-length(x)])))

as.numeric((ind2$full.start[1]-1 - ind1$full.start[1]), "mins")
## add minutes as 0s to start!!!
plot(c(seq(0,10.6244*60-1, by=1)/60, (10.6244*60-1)/60+ (0:seclength)/60), c(rep(0,10.6244*60), mean.sim.gaps), type="l", bty="l", las=1, xlab="minutes start shifted since obs start of U", ylab="mean sim gaps")


par(mfrow=c(1,1))
hist(mean.sim.gaps, xlim=c(max(0, min(mean.obs.gaps, mean.sim.gaps)), max(c(mean.obs.gaps, mean.sim.gaps))), main=paste(unique(ind2$ind), "after", unique(ind1$ind)), xlab="Simulated inter-individual nearest chest beat time gaps (secs)", color="grey95", border="grey70", xaxs="i", yaxs="i")
# hist(mean.sim.gaps, xlim=c(max(0, min(mean.obs.gaps, mean.sim.gaps)), max(c(mean.obs.gaps, mean.sim.gaps))), main=paste0(unique(ind2$ind), " after ", unique(ind1$ind), "\nexcluding ",unique(ind2$ind),"'s last chest beat"), xlab="Simulated inter-individual nearest chest beat time gaps (secs)", color="grey95", border="grey70", xaxs="i", yaxs="i")

abline(v=mean.obs.gaps, lty=4, col="red", lwd=2)
(pval <- rank(c(mean.obs.gaps, mean.sim.gaps))[1]/length(unique(rank(c(mean.obs.gaps, mean.sim.gaps)))))
text(mean.obs.gaps-5, 250, labels=paste("p = ", round(pval, 3)), col="red")


# tmp <- c(short, long, nrow(ind2), nrow(ind1), mean.obs.gaps, seclength, pval)
# num <- rbind(num, tmp)
# names(num) <- c("FocalInd", "BkgdInd", "num_cbs_focal", "num_cbs_bkgd", "mean_gap", "num_sims", "pval")
