
library(monitoR)
ogdir <- "~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/"
setwd(ogdir)

## Using this to check along with Raven correlations. The new method looks better! And is at least easier...
fieldLocs <- read.csv("fieldTraingLocs.csv")
gorilla <- fieldLocs[19,c("X","Y")] #the gorilla in question, swamp big guy
pams2check <- pam.xy[rownames(pam.xy)%in%c("J","Q"),]
pams2check$dist2gorilla <- distGorMic(gorilla, pams2check)
pams2check$time2gorilla <- pams2check$dist2gorilla/343
(lagExp = pams2check$time2gorilla[2] - pams2check$time2gorilla[1])


## Make detection tables into cross-correlated lags

lags <- xcorrlags(CBs, len, frq.lim = c(0.2, 0.7), wl = 2048, ovlp = 50)

## Set up for triangulation:

lags.xy <- 


## Triangulate from lags:

goriLoc(lags.xy)






