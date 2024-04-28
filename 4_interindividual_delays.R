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
dec10 <- readRDS("2024.04.15_20221210_individuals_E_J_V_Mdef_Mmb_Dpok_Dsm_Vclap/20240414_132623_CBs.rds")
jan16 <- readRDS("2024.04.21_20230116_individuals_B1_B2_D1_D2_V/20240421_213648_CBs.rds")
jan17 <- readRDS("2024.04.22_20230117_individuals_D17_D18_K17_K18_V17_V18/20240422_215708_CBs.rds")

# compute consecutive gaps:
# should be based on one PAM...or both their closest so even bias? or which one? ************** FOR LATER

# keep each individual separate, calculate closest in other individual df...

start.date = "20230206"

ind1 <- feb6[[1]]$detections
ind2 <- feb6[[2]]$detections
ind1$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind1$sound.files, clip.start = ind1$startClip.GPS)
ind2$full.start <- as.POSIXct(start.date, format = "%Y%m%d") + fullTimeFromClipStart(sound.path = ind2$sound.files, clip.start = ind2$startClip.GPS)

# Plot the chest beat timeline colored by individual:

par(mfrow=c(4,1))
#pdf("2024.04.07_20230206_individuals_J_U/timeline.pdf", width=12, height=3)
plot(x=ind1$full.start, y=rep(1, length(ind1$full.start)), axes=F, pch="|", xlab=NA, ylab=NA, cex=3, main=as.Date(ind1$full.start[1]-86400))
axis.POSIXct(side = 1, tick=F)
abline(h=1)
points(x=ind2$full.start, y=rep(1, length(ind2$full.start)), pch="|", col="red", cex=3)
mtext(text = "Time gaps between chest beats from individuals J and U", cex=0.75)
dev.off()
# Get overall first and last CB of all individuals that night:

cb1 <- min(c(min(ind1$full.start), min(ind2$full.start)))
cblast <- max(c(max(ind1$full.start), max(ind2$full.start)))

# slide CBs of individual who lasts less long within this time range:

# for each individual:
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
title="U shifted over J"; setx="Minimum time between CBs of different individuals (minutes)"
hist(abs.gaps.sim, bty="l", las=1, breaks = quantile(abs.gaps.sim, probs = c(0,0.01,seq(0.1,0.9,0.1),1)),
     xlab=setx, main=title)
#or plot the gamma distribution:
pdf("2024.04.07_20230206_individuals_J_U/distributionSimGaps.pdf", width=8, height = 6)
plot(density(gamma.dens <- rgamma(n=3000,shape=gamma$estimate[1],rate=gamma$estimate[2]), from = 0), yaxs="i", xaxs="i",
     bty="l", xlab=setx, main=title,las=1, xlim=c(0,50), zero.line=FALSE, ylab="Gamma density of simulated values")
#shade lower 99% CI gamma:
abline(v=c(ci95gamma[1]-seq(ci95gamma[1],min(abs.gaps.sim),-0.001)), col=transp("lightblue",0.01),lwd=2)
#shade lower 99% CI log-normal:
# abline(v=c(ci95lnorm[1]-seq(ci95lnorm[1],min(abs.gaps.sim),-0.001)), col=transp("lightgreen",0.006),lwd=2)
#add lines for obs data:
abs.min.gaps.obs <- abs(as.numeric(minIIgaps)/60)
points(x=abs.min.gaps.obs, y=rep(0, length(minIIgaps)), col="red", pch=20, xpd=T)

# all are significantly replies:
length(which(abs.min.gaps.obs<ci95lnorm[1]))/length(abs.min.gaps.obs)
length(which(abs.min.gaps.obs<ci95gamma[1]))/length(abs.min.gaps.obs)

text(5, 0.005, labels=paste0("n=", length(which(abs.min.gaps.obs<ci95gamma[1])), "/",length(abs.min.gaps.obs)," of U's chest beats are closer to J's than \n<5% of the simulated inter-individual time gaps"), col="red", adj=0)

dev.off()


 minIIgapsSim <- lapply(starts, function(x) min(abs(x[1] - ind2$full.start)))
for (i in 2:nrow(ind1)){
  minIIgapsSim <- lapply(starts, function(x) min(abs(x[i] - ind2$full.start)))
}

# Check for within-individual escalation:

plot(1:length(indTgaps), indTgaps, las=1, xlab="Chest beat number", ylab="Time gap", bty="l", pch=20)

# Check for btwn-individual escalation:

plot(1:length(minIIgaps), minIIgaps, las=1, xlab="Chest beat number", ylab="Time gap", bty="l", pch=20)



plotConseq(rbind(feb6[[1]]$detections,feb6[[2]]$detections))

# but make it just a timeline

# compute the metric: nearest-neighbor inter-individual gaps (doesn't matter which individual? or just pairwise?)

# then find the start time and the end time of all chest beating

# then slide each around by 1 second randomly, if it gets to the end, loop it back to the beginning so it's not ordered which is weird for escalation but makes sense for this simulation

# compute the metric

# do that 1000 times

# plot them all, color the real ones - are they in the tail?
