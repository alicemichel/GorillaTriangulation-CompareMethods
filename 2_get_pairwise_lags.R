
#source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/1_combine_multi_cbs_from_Raven_detection_table.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_prep4Raven.R")

# Above should be set from running previous

## Make detection tables into cross-correlated lags

# Depending how many CBs you got and how many individuals there might be,
# tell it which CBs to use for localization that look like the same individuals
# don't use all the ones that you were using for consecutive plots EXCEPT for the ones you know (close, no possible movement)

CBs4raven <- CBs

par(mfrow=c(4,1))
check <- CBs[["D"]]$detections
#check_spectro(check, rownames(check))
Dpok <- c(63,4,32) 
DpokIds <- which(rownames(CBs4raven[["D"]]$detections) %in% Dpok)
CBs4raven[["Dpok"]]$detections <- CBs4raven[["D"]]$detections[DpokIds,]
CBs4raven[["Dpok"]][4:10] <- lapply(CBs4raven[["D"]][4:10],FUN = function(x) x[DpokIds,])
names(CBs4raven[["Dpok"]])[4:10] <- names(CBs4raven[["D"]])[4:10]

Dsm <- c(6,1,8,13,23,28)
DsmIds <- which(rownames(CBs4raven[["D"]]$detections) %in% Dsm)
CBs4raven[["Dsm"]]$detections <- CBs4raven[["D"]]$detections[DsmIds,]
CBs4raven[["Dsm"]][4:10] <- lapply(CBs4raven[["D"]][4:10],FUN = function(x) x[DsmIds,])
names(CBs4raven[["Dsm"]])[4:10] <- names(CBs4raven[["D"]])[4:10]


par(mfrow=c(4,1))
check <- CBs[["M"]]$detections
#check_spectro(check, rownames(check))
Mmb <- c(68,70,75,73,77,79)
MmbIds <- which(rownames(CBs4raven[["M"]]$detections) %in% Mmb)
CBs4raven[["Mmb"]]$detections <- CBs4raven[["M"]]$detections[MmbIds,]
CBs4raven[["Mmb"]][4:10] <- lapply(CBs4raven[["M"]][4:10],FUN = function(x) x[MmbIds,])
names(CBs4raven[["Mmb"]])[4:10] <- names(CBs4raven[["M"]])[4:10]

Mdef <- rownames(check)[!rownames(check) %in% Mmb]
MdefIds <- which(rownames(CBs4raven[["M"]]$detections) %in% Mdef)
CBs4raven[["Mdef"]]$detections <- CBs4raven[["M"]]$detections[MdefIds,]
CBs4raven[["Mdef"]][4:10] <- lapply(CBs4raven[["M"]][4:10],FUN = function(x) x[MdefIds,])
names(CBs4raven[["Mdef"]])[4:10] <- names(CBs4raven[["M"]])[4:10]

dev.off()
par(mfrow=c(4,1))
check <- CBs[["Vclap"]]$detections
check_spectro(check, rownames(check), clipLength = 8)

CBs4raven <- CBs4raven[c("E","J","V","Mdef","Mmb","Dpok","Dsm","Vclap")]
names(CBs4raven) <- inds <- c("E","J","V","Mdef","Mmb","Dpok","Dsm","Vclap")
binds <- data.frame()
for (x in 1:length(CBs4raven)){
  CBs4raven[[x]]$detections$ind <- names(CBs4raven)[x]
  binds <- rbind(binds, CBs4raven[[x]]$detections)
}
plotConseq(binds)


# go back to the previous to remake the plot and CBs if this doesn't work but first try editting CBs


## 1. Export CBs4raven## 1. Export for Raven (or Python):
ravenprep(CBs4raven, len, inds) #, detlim = c(23,NA)
#rm(list = ls())

## 2. Or, do cross-correlation in R:
# lags <- xcorrlags(CBs, len, frq.lim = c(0.2, 0.7), wl = 2048, ovlp = 50)