
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/1_combine_multi_cbs_from_Raven_detection_table.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_prep4Raven.R")

# Above should be set from running previous

## Make detection tables into cross-correlated lags

## 1. Export for Raven (or Python):
ravenprep(CBs, len, detlim = c(23,NA))
rm(list = ls())

## 2. Or, do cross-correlation in R:
# lags <- xcorrlags(CBs, len, frq.lim = c(0.2, 0.7), wl = 2048, ovlp = 50)