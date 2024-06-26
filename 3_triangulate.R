
# Need to do this before, don't need to run necessarily
# source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/2_get_pairwise_lags.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_goriLoc_Damien_edit.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_distance_calcs.R")


# should already be in the correct clock-drift time
# and the right time gaps
# making the lags "real-time"
# so, theoretically, with only those we should be able to triangulate...

setwd("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/2024.04.22_20230117_individuals_D17_D18_K17_K18_V17_V18/")
(lags <- readxl::read_excel(list.files(pattern=".xlsx")))

## Set up for triangulation:

date <- unique(substr(list.files(pattern=".xlsx"), start = 1, stop = 8))[1]
xy <- read.csv("../xy2.csv", row.names=1)[,1:2] #csv of all the pams locations with first column for their name

## Run triangulation and save output as lists and inspect plots:

## Using this to check along with Raven correlations. The new method looks better! And is at least easier...
fieldLocs <- read.csv("../fieldTraingLocs.MethodsComparison.csv")
#points(fieldLocs[,3])
ns = 5
gorilla <- fieldLocs[ns,c("X","Y")] #the gorilla in question, near PAM J

localizedSBs <- list()
optima <- data.frame()
#optima <- read.csv("appendingLocs.csv")
for (sb in unique(lags[!is.na(lags$lag),]$IndID)){
  lags1ind <- lags[lags$IndID==sb,]
  lags1ind$lag <- -(as.numeric(lags1ind$lag))
  temperature = 18
  #pdf(paste0(sb,"loc.pdf"), height=8, width=10)
  localizedSBs[[sb]] <- goriLoc(lags1ind, xy, main=date, temperature = temperature, xextr = 2000)
  #dev.off()
  points(gorilla, pch=25, cex=1.5)
  for (i in 1:length(localizedSBs[[sb]]$hyperbola)){
    points(localizedSBs[[sb]]$hyperbola[[i]], col="blue",cex=0.1)
  }
  
  optima <- rbind(optima, c(date, sb, localizedSBs[[sb]]$optimum))
}
names(optima) <- c("date", "SB","lon","lat")

pdf(paste0(date,"alllocs.pdf"), height=8, width=10)
plot(xy, las=1, xlab=NA, ylab=NA, bty="l", xlim=c(min(xy$lon)-500, max(xy$lon)+500),ylim=c(min(xy$lat)-500, max(xy$lat)+500), pch=20)
points(as.numeric(optima$lon), as.numeric(optima$lat), pch=25)
text(as.numeric(optima$lon), as.numeric(optima$lat), labels = optima$SB, pos = 1)
dev.off()
## Compare to field localizations for accuracy metrics:

## Using this to check along with Raven correlations. The new method looks better! And is at least easier...
fieldLocs <- read.csv("../fieldTraingLocs.MethodsComparison.csv")
#points(fieldLocs[,3])
ns = rep(5,6)
for (sbi in 1:length(unique(lags[!is.na(lags$lag),]$IndID))){
  print(unique(lags[!is.na(lags$lag),]$IndID)[sbi])
  field.prec <- fieldLocs[ns[sbi],]$Precision.m
  gorilla <- fieldLocs[ns[sbi],c("X","Y")]
  dist2gorilla <- distGorMic(gorilla, localizedSBs[[sbi]]$intersection)
  dist2gorillaOpt <- ifelse(!is.na(localizedSBs[[sbi]]$optimum), distGorMic(gorilla, localizedSBs[[sbi]]$optimum), NA)
  plot(dist2gorilla, xaxt="n", xlab=NA, ylim=c(0,50), bty="l", las=1, ylab="Distance to actually found nest site", main=unique(lags[!is.na(lags$lag),]$IndID)[sbi])
  axis(1, at = 1:length(dist2gorilla), labels = localizedSBs[[sbi]]$intersection$pams, las=2)
  cat("Starting with individual...", fieldLocs[ns[sbi],]$Nest.Site.ID, "\n")
  cat("New method best:", "\n")
  print(min(dist2gorilla))
  cat("which was from PAM combo...")
  win <- localizedSBs[[sbi]]$intersection[which(dist2gorilla==min(dist2gorilla)),]$pams
  cat(win)
  points(localizedSBs[[sbi]]$intersection[which(dist2gorilla==min(dist2gorilla)),], col="darkred", pch=16)
  cat("\n3 closest pams would be...")
  top3 <- rownames(xy[order(distGorMic(gorilla, xy)),])[1:3]
  cat(top3, "...same?", all(stringr::str_detect(win,top3)))
  cat("\n\nField best:", "\n")
  print(field.prec)
  cat("New method summary:", "\n")
  print(summary(dist2gorilla))
  cat("New method KDE optimum (weights points within xextr of grid equally):", "\n")
  print(dist2gorillaOpt[[1]])
  cat("\n\n\n")
  
  fieldLocs[ns[sbi],]$new.method.temp.set <- temperature
  fieldLocs[ns[sbi],]$new.method.KDEopt <- dist2gorillaOpt[[1]]
  fieldLocs[ns[sbi],]$new.method.best <- min(dist2gorilla)
  fieldLocs[ns[sbi],]$new.method.best.pams <- win
  fieldLocs[ns[sbi],]$closest.3.pams <- paste0(top3, collapse="")
  fieldLocs[ns[sbi],]$new.method.xcors <- list.files(pattern=".xlsx")[1]
  
}
fieldLocs[ns,]$office.notes <- "another near K and another near V (moving maybe or maybe 2)."
#write.csv(fieldLocs, "../fieldTraingLocs.MethodsComparison.csv")


pams2check <- pam.xy[rownames(pam.xy)%in%c(unique(lags1ind$PAM)),]
pams2check$dist2gorilla <- distGorMic(gorilla, pams2check)
pams2check$time2gorilla <- pams2check$dist2gorilla/340
pams2check$lagExp = pams2check$time2gorilla - min(pams2check$time2gorilla)
pams2check

#fieldLocs <- fieldLocs[,-1]
# pdf("../fieldlocs.pdf", width=10, height=7)
par(mar=c(2,4,2,8))
plot(xy, las=1, xlab=NA, ylab=NA, xlim=range(fieldLocs$X)+c(-100,100), ylim=range(fieldLocs$Y)+c(-100,100), col="grey20", pch=5, cex=1, main="Field-validated localizations to test new method accuracy")  
text(fieldLocs$X, fieldLocs$Y, labels = fieldLocs$Nest.Site.ID, pos=1, cex=0.7)
text(xy, labels = rownames(xy), cex=0.4)
points(fieldLocs$X, fieldLocs$Y, col="red", pch=16)
legend("right", title="Recorded on:", legend=c(paste0(substr(fieldLocs$Nest.Site.ID, start=1, stop=6), "_", as.Date(fieldLocs$date.gorillas.left.site, format = "%m/%d/%y")-1), rep("",4)), bty="n", xpd=T, inset=-0.22, cex=0.9)
# dev.off()
 

