
# Need to do this before, don't need to run necessarily
# source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/2_get_pairwise_lags.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_goriLoc_Damien_edit.R")

# should already be in the correct clock-drift time
# and the right time gaps
# making the lags "real-time"
# so, theoretically, with only those we should be able to triangulate...

setwd("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/2024.04.07_20230206_individuals_J_U")
(lags <- readxl::read_excel(list.files(pattern=".xlsx")))

## Set up for triangulation:

date <- unique(substr(list.files(pattern=".xlsx"), start = 1, stop = 8))[1]
xy <- read.csv("../xy2.csv", row.names=1)[,1:2] #csv of all the pams locations with first column for their name

## Run triangulation and save output as lists and inspect plots:

localizedSBs <- list()
for (sb in unique(lags$IndID)){
  lags1ind <- lags[lags$IndID==sb,]
  lags1ind$lag <- -lags1ind$lag
  localizedSBs[[sb]] <- goriLoc(lags1ind, xy, main=date, temperature = 23)
}


## Compare to field localizations for accuracy metrics:

# ## Using this to check along with Raven correlations. The new method looks better! And is at least easier...
# fieldLocs <- read.csv("fieldTraingLocs.csv")
# gorilla <- fieldLocs[19,c("X","Y")] #the gorilla in question, swamp big guy
# pams2check <- pam.xy[rownames(pam.xy)%in%c("J","Q"),]
# pams2check$dist2gorilla <- distGorMic(gorilla, pams2check)
# pams2check$time2gorilla <- pams2check$dist2gorilla/343
# (lagExp = pams2check$time2gorilla[2] - pams2check$time2gorilla[1])
# # Do this for all of them...later...




