
## General setup
ogdir <- "~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/"
setwd(ogdir)
library(warbleR)
library(monitoR)
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_distance_calcs.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_get_pull_times.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_clock_drift.R")

par(mfrow=c(1,1))

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

dets <- read.csv("20230206_detections_Q.csv")[,1:14]
#dets$approxOrd <- approxOrd(dets)

buff = 1 #Buffer to add to each clip if they're tight
dets$startClip <- dets$startClip - buff

# Convert clip start times to GPS time to correct for clock drift within each hour-long file:
dets$startClip.GPS <- hz2GPStime(clipStart = dets$startClip, soundpath = dets$sound.files)

inds <- unique(substr(dets$approxOrd, start = 1, stop = 1))[1:2]

dets <- dets[dets$pam %in% c("J","H","Q","U"),] #for now, use 6 closest in future...within 2nd loop below
pamID <- unique(dets$pam) # change to reflect sound files present in the folder

start.date <- substr(dets[1,]$sound.files, start = 4, stop = 11)

len = 7 #make >7 if hoots are included, since it'll shift the others way back (hoots up to like 6 secs pre-CB...)
# To get len, take the max selection length. Add the max expected lag based on the max distance. Then add the expected length of a distant CB (2.5sec) to that...

printwindow = 40

CBs <- as.list(inds)
names(CBs) <- inds
par(mfrow=c(length(inds),1))
for (c in 1:length(CBs)){ 
  
  ## isolate the detections of individual "c":
  detections.of.ind <- CBs[[c]]$detections <- dets[substr(dets$approxOrd, start = 1, stop = 1) == inds[c],]
  
  ## pick 1 PAM and fix the time gaps based on ALL its detections. This is the PAM you should have the max number of detections for that individual:
  key.pam <- names(CBs)[c] #here the CBs are handily named by the closest PAM
  dets.key.pam <- detections.of.ind[detections.of.ind$pam==key.pam,]
  
  # Get the start time of the first detection within the file:
  first.clip.t.key <- dets.key.pam[1,]$startClip.GPS
  # Get the file path:
  first.clip.sound.path.key <- dets.key.pam[1,]$sound.files
  # Use the file path to know when the file actually starts since midnight the night before:
  first.clip.sound.path.key.start <- fullTimeFromClipStart(first.clip.sound.path.key,0)
  
  # last file path (for double checking below):
  last.clip.sound.path.key <- dets.key.pam[nrow(dets.key.pam),]$sound.files
  ch <- ifelse(fullTimeFromClipStart(last.clip.sound.path.key,0) - first.clip.sound.path.key.start > 0, "all detections are before midnight \n", "you need to double check if the files got messed up by carrying over to the next day \n")
  cat(ch)
  
  # Read in and combine waves
  wavs <- list()
  time.gaps <- vector()
  for (i in 1:nrow(dets.key.pam)){
    wavs[[i]] <- read_wave(dets.key.pam[i,]$sound.files, from = dets.key.pam[i,]$startClip, to = (dets.key.pam[i,]$startClip + len))
    clip.sound.path.key.next.start.gap.from.first <- fullTimeFromClipStart(dets.key.pam[i+1,]$sound.files,0) - first.clip.sound.path.key.start
    (time.gaps[i] <- dets.key.pam[i+1,]$startClip.GPS - (first.clip.t.key) + clip.sound.path.key.next.start.gap.from.first)
  }
  CBs[[c]]$w.keypam <- do.call(tuneR::bind, args = wavs)
  # Get the real-time start gaps between each detection on the key PAM:
  time.gaps <- c(0, time.gaps)[-(length(time.gaps)+1)]
  
  ## Plot the accumulation of CBs over time - interesting for plotting exchanges iff you select EVERY chest beat:
  plot(time.gaps/60, 1:length(time.gaps), bty="l", pch=21, las=1, type="b", lwd=1.75, xlab="Elapsed time (min)", ylab="Cumulative # chest beats", main=names(CBs)[c])
  
  # get 6 closest Pam’s to key Pam by distance:
  (pamsToCheck <- pamNeighbors(key.pam, pam.xy, n=10, self=TRUE) )
  
  # of those, only those you have sound files for:
  (pamsToGet <-  pamsToCheck[pamsToCheck %in% pamID]) #include key PAM as a check
  
  # loop through each PAM for this gorilla:
  CBs[[c]]$w.pam <- list()
  for (pam in pamsToGet){ 
    
    # Get start times data.frame within the focal PAM based on time gaps in the key PAM:
    CBs[[c]]$w.pam <- mergedClipsFromTimeGaps(clip.start=first.clip.t.key, inPAMfile=first.clip.sound.path.key, gaps=time.gaps, getPAM=pam, clipLength=len, path=".", keyPAMstarts4check = dets.key.pam$startClip, tabulate = TRUE, add=10)
    
    wave="sure"
    if (wave != "NO"){
      
      mergWvFoc <- df2mergedWav(df = CBs[[c]]$w.pam, clipLength=len)
      
      plot.check = "NO"
      if (plot.check!="NO"){

        par(mfrow=c(2,1))
        ## Amalgamated spectro in key PAM:
        viewSpec(CBs[[c]]$w.keypam, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning", main = key.pam, page.length = printwindow)
        ## Amalgamated spectro in focal PAM:
        viewSpec(mergWvFoc, interactive = F, start.time = len, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning", main = pam, page.length = printwindow)
      }
    }
    
    names(CBs[[c]])[names(CBs[[c]])=="w.pam"] <- paste0("w.",pam)
    
    writeWave(mergWvFoc, paste0("merged_", start.date, "_ind_", names(CBs)[c], "_on_", pam, ".wav"))
    
  }
  
}





# for cross-correlation, this may be way too much data... python? or
# cut off the amplitude to low (not 0) up to the start of the clipped clip per each pam...still minutes of data...
## Alternatively could incorporate selection of delay into the function manually instead of wide length, but need to paste in emptiness and gets complicated








