
## General setup
ogdir <- "~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/"
setwd(ogdir)
library(warbleR)
library(monitoR)
par(mfrow=c(2,1))

## Strategy:

## 1. Read in detection tables from R
##  a. These could be the tables from the field, created from inputting one PAM and getting detections matching on others.
##  b. Alternatively they could come from Raven. Just put them through your R selection to get them into the same format with the additional R columns.
## 2. Split the table into individual gorillas, by the "Approx Order" column you made, assuming no two gorillas are closest to the same PAM (focus on within grid...).
##  a. Optionally, before splitting them add a buffer to the start of each unique-CB ('buff')
## 3. Loop the above through each gorilla's detections data table.
##  a. Choose a key PAM. Currently it's based on the selected individual IDs, but you will do this using Approx Order to automate it.
##  b. Glue together all the detections from that individual on the key PAM together (set "len" for the unique-CB length).
##  c. Record the time gaps between each unique-CB from that PAM.
##  d. Then loop through the remaining PAMs (actually all even the key PAM, so as to check the gaps didn't get messed):
##    i. Using the first detection on the specific PAM as a base and the time gaps from the key PAM, glue together the same time fragments for each PAM
##   ii. Save as a list to use in the next step (cross-correlation). Print the spectrograms to verify no errors. 

pam.xy <- read.csv("xy2.csv", row.names=1)[,1:2]

dets <- read.csv("20230206_detections_Q.csv")
buff = 1 #Buffer to add to each clip if they're tight
dets$startClip <- dets$startClip - buff

inds <- unique(substr(dets$approxOrd, start = 1, stop = 1))[1:2]

dets <- dets[dets$pam %in% c("J","H","Q","U"),] #for now, use 6 closest in future...within 2nd loop below


pamID <- unique(dets$pam)

len = 5.5

CBs <- as.list(inds)
names(CBs) <- inds
for (c in 1){ #1:length(CBs)
  detections.of.ind <- CBs[[c]]$detections <- dets[substr(dets$approxOrd, start = 1, stop = 1) == inds[c],]
  
  ## pick 1 PAM and fix the time gaps based on ALL its detections. This is the PAM you should have the max number of detections for that individual.
  key.pam <- names(CBs)[c]
  dets.key.pam <- detections.of.ind[detections.of.ind$pam==key.pam,]
  
  first.clip.t.key <- dets.key.pam[1,]$startClip
  first.clip.sound.path.key <- dets.key.pam[1,]$sound.files
  first.clip.sound.path.key.start <- fullTimeFromClipStart(first.clip.sound.path.key,0)
  
  wavs <- list()
  time.gaps <- vector()
  for (i in 1:nrow(dets.key.pam)){
    wavs[[i]] <- read_wave(dets.key.pam[i,]$sound.files, from = dets.key.pam[i,]$startClip, to = (dets.key.pam[i,]$startClip + len))
    clip.sound.path.key.next.start.gap.from.first <- fullTimeFromClipStart(dets.key.pam[i+1,]$sound.files,0) - first.clip.sound.path.key.start
    (time.gaps[i] <- dets.key.pam[i+1,]$startClip - (first.clip.t.key + len) + clip.sound.path.key.next.start.gap.from.first)
  }
  CBs[[c]]$w.keypam <- do.call(tuneR::bind, args = wavs)
  # viewSpec(w, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning", main = key.pam, page.length = duration(w))
  time.gaps <- c(0, time.gaps)[-(length(time.gaps)+1)]
  plot(time.gaps/60, 1:15, bty="l", las=1, type="l", lwd=2, xlab="Elapsed time (min)", ylab="Cumulative # chest beats")
  
  CBs[[c]]$w.pam <- list()
  
  # get 6 closest Pam’s to key Pam by distance
  (pamsToCheck <- pamNeighbors(key.pam, pam.xy, n=10, self=TRUE) )
  (pamsToGet <-  pamsToCheck[pamsToCheck %in% pamID])
  
  for (pam in pamsToGet){ 
    
    # first print the one you're aiming for
    par(mfrow=c(2,1))
    viewSpec(read_wave(inPAMfile, from = clip.start, to = clip.start + len), main = key.pam, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning")
    # prints the clip and saves the clip starts/wave read in into a list
    pam.clip <- pullTime(first.clip.t.key, first.clip.sound.path.key, pam, len) #change len if you have reason to believe they are very far apart
    
    # now for the next ones need to check every start time and what file it falls in... or go up from the starting one? that makes more sense I think since then you won't lose if it loops over into the next day...if you play it right...
    
    
    # must include the first time (negative length)
    time.gaps
    
    
    
    
    
    
    
    
    
    
    # get raw time, then find file within folder for that range, could be two in which case…skip but account for the gap…..paste in empty recording of same length of time from that pam/region. Ensure empty is issue 
    
    
    
    
    
    
    
    wavs <- list()
    for (i in 1:nrow(dets.pam)){
      wavs[[i]] <- read_wave(dets.pam[i,]$sound.files, from = dets.pam[i,]$start.based.on.key, to = (dets.pam[i,]$start.based.on.key + len))
    }
    CBs[[c]]$w.pam[[pam]] <- do.call(bind, args = wavs)
    # viewSpec(w.p, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning", main = pam, page.length = duration(w.p))
    
  }
}

# writeWave(w.p, "testH.wav")



# for cross-correlation, this may be way too much data... python? or
# cut off the amplitude to low (not 0) up to the start of the clipped clip per each pam
