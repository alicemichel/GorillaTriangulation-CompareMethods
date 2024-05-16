
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
night1 <- feb6
# how many CBs per individual?
lapply(night1, function(x) nrow(x$detections))
# put the shorter one as individual 2 and run the function:
simGapsListCBs(night1,short=2,long=1,sign="neg")

