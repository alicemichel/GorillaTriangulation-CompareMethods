
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/1_combine_multi_cbs_from_Raven_detection_table.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_prep4Raven.R")

# Above should be set from running previous

## Make detection tables into cross-correlated lags

# Depending how many CBs you got and how many individuals there might be,
# tell it which CBs to use for localization that look like the same individuals
# don't use all the ones that you were using for consecutive plots EXCEPT for the ones you know (close, no possible movement)

CBs4raven <- CBs

par(mfrow=c(4,1))
check <- CBs[["D_17"]]$detections
#check_spectro(check, rownames(check))

# go back to the previous to remake the plot and CBs if this doesn't work but first try editting CBs


## 1. Export CBs4raven## 1. Export for Raven (or Python):
ravenprep(CBs4raven, clipLength = 6, inds) #, detlim = c(10,NA,NA,NA,NA)
rm(list = ls())

## 2. Or, do cross-correlation in R:
# lags <- xcorrlags(CBs, len, frq.lim = c(0.2, 0.7), wl = 2048, ovlp = 50)