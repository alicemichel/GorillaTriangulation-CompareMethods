
## General setup
ogdir <- "~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/"
setwd(ogdir)
library(warbleR)
library(monitoR)
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_relabel_files.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_distance_calcs.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_get_pull_times.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_clock_drift.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_raven_approx_ord.R")
source("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/CrossCorrMethodsComparison/functions_cut_and_buffer.R")

par(mfrow=c(1,1))

## Strategy:

## 1. Read in detection tables from Raven

## 2. Split the table into individual gorillas, by the "Approx Order" column, assuming no two gorillas are closest to the same top 3 PAMs (focus on within grid...).

## 3. Loop the above through each gorilla's detections data table.
##  a. Choose a key PAM that the individual was closest to.
##  b. Glue together all the detections from that individual on the key PAM together (set "len" for the unique-CB length).
##  c. Record the time gaps between each unique-CB from that PAM -- in "true" (GPS) time.

##  d. Then loop through the remaining PAMs (actually all even the key PAM, so as to check the gaps didn't get messed):
##    i. Using the first detection on the specific PAM as a base and the time gaps from the key PAM, glue together the same time fragments for each PAM. Importantly, the times need to first be converted within each soundfile to the "false" (sample) time, essentially reverse clock-drifting the time PER PAM. Run this by someone else to see if it makes sense...


##   ii. Save as a list to use in the next step (cross-correlation). Print the spectrograms to verify no errors. 


## prep sound files to do detections in Raven on the best PAM(s) per individual (may need to iteratively run approxOrd function to figure out which that is)

#relabelBARfiles(PAM="R")
#grabLaCieBAR(date="S20230116", times=c(16:22), lacienumber = 8, destn=".") #make sure you're in the right directory!

pam.xy <- read.csv("xy2.csv", row.names=1)[,1:2]

#dets.long <- approxOrd(Raven.selections.path = "20230116_V.txt", buffer=2, clipLength = 6)
#saveRDS(dets.long, "dets20240418.rds")
dets.long <- readRDS("dets20240418.rds")

#writexl::write_xlsx(dets.long, "20221210.xlsx")
#dets.edit <- as.data.frame(readxl::read_excel("20221210.xlsx"))
#dets.long <- dets.edit[is.na(dets.edit$check) | dets.edit$check!="nothing",]

unique(dets.long$note)
dets.D <- dets.long[!is.na(dets.long$note) & dets.long$note %in% c("DBNP box weak", "DBNP box"), ] #remove "D, but no other files", for everything but the plot
dets.D$ind="D"
dets.VO <- dets.long[!is.na(dets.long$note) & dets.long$note == "VO" & dets.long$pam=="V", ]
dets.B <- dets.long[!is.na(dets.long$note) & dets.long$note %in% c("BD","BD hoot and ground", "BD box", "BD hoot", "BD, and/or check U", "BD box weak", "BD weak", "BD box, but cut off so maybe not, check N", "BD low"), ]
dets.B$ind="BD"
dets.long <- rbind(dets.B, dets.D, dets.VO)

# take the cut column, subtract from start time, add the buffer
dets.long$min.cut = ifelse(is.na(dets.long$min.cut), 2, dets.long$min.cut)
dets.long$startClip <- cutNbuff(dets.long$start, dets.long$min.cut, buffer=4) #6 for V, 4 for Bs and Ds

# Cleaning particular to each night
dets.long$check.full.st <- (fullTimeFromClipStart(sound.path = dets.long$sound.files, clip.start = dets.long$startClip))
dets.long$check.ind <- substr(dets.long$sound.files, start=1, stop=1)==dets.long$ind

dets <- dets.long

#split by time since there's a huge gap
dets$ind <- ifelse(dets$ind=="D" & dets$check.full.st>70148.70, "D2", ifelse(dets$ind=="D", "D1", dets$ind))
dets$ind <- ifelse(dets$ind=="BD" & dets$check.full.st>70148.70, "B2", ifelse(dets$ind=="BD", "B1", dets$ind))
dets[rownames(dets)==77,]$ind <- "D1"

# Convert clip start times to GPS time to correct for clock drift within each hour-long file:
dets$startClip.GPS <- hz2GPStime(clipStart = dets$startClip, soundpath = dets$sound.files)

# Order the detections
dets <- dets[order(fullTimeFromClipStart(sound.path = dets$sound.files, clip.start = dets$startClip)),]

(inds <- sort(unique(dets$ind)))
#inds <- c("D", "E", "M", "V")
pamID <- unique(substr(list.files(".", pattern = "S20.*.wav"), start = 1, stop = 1)) # reflects sound files present in the folder

start.date <- substr(dets[1,]$sound.files, start = 4, stop = 11)

len = 6 #make >7 if hoots are included, since it'll shift the others way back (hoots up to like 6 secs pre-CB...)
# To get len, take the max selection length. Add the max expected lag based on the max distance. Then add the expected length of a distant CB (2.5sec) to that...

printwindow = 40

# plot the consecutive chest beats of each individual over one another
pdf(paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_consec.pdf"), width=12, height=8)
plotConseq(dets)
dev.off()
plotConseq(dets)

doyouhavemultplpamsperCB="NO"

CBs <- as.list(inds)
names(CBs) <- inds
par(mfrow=c(length(inds),1), mar=c(4,4,2,0))
for (c in 1:length(CBs)){ 
  print(c)
  
  ## isolate the detections of individual "c":
  detections.of.ind <- CBs[[c]]$detections <- dets[dets$ind == inds[c],]
  
  ## remove rows in between files so they don't mess up xcorr due to length diffs, as commanded below
  if (c==1){
    detections.of.ind <- CBs[[c]]$detections <- detections.of.ind[-11,]
  }
  if (c==5){
    detections.of.ind <- CBs[[c]]$detections <- detections.of.ind[-2,]
  }
  
  ## pick 1 PAM and fix the time gaps based on ALL its detections. This is the PAM you should have the max number of detections for that individual:
  key.pam <- substr(names(CBs)[c], start=1, stop=1) #here the CBs are handily named by the closest PAM
  
  dets.key.pam <- detections.of.ind
  ############## IMPORTANT: ABOVE CHANGE BACK IF YOU HAVE MULTIPLE DETECTIONS
  if (doyouhavemultplpamsperCB=="yes"){
    dets.key.pam <- detections.of.ind[detections.of.ind$pam==key.pam,]
  }
  
  # Get the start time of the first detection within the file:
  first.clip.t.key <- dets.key.pam[1,]$startClip.GPS
  # Get the file path:
  first.clip.sound.path.key <- dets.key.pam[1,]$sound.files
  # Use the file path to know when the file actually starts since midnight the night before:
  first.clip.sound.path.key.start <- fullTimeFromClipStart(first.clip.sound.path.key,0)
  
  # last file path (for double checking below):
  last.clip.sound.path.key <- dets.key.pam[nrow(dets.key.pam),]$sound.files
  ch <- ifelse(fullTimeFromClipStart(last.clip.sound.path.key,0) - first.clip.sound.path.key.start >= 0, "all detections are before midnight \n", "you need to double check if the files got messed up by carrying over to the next day \n")
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
  
  # get closest Pam’s to key Pam by distance:
  (pamsToCheck <- pamNeighbors(key.pam, pam.xy, n=6, self=TRUE) )
  
  # of those, only those you have sound files for:
  (pamsToGet <-  pamsToCheck[pamsToCheck %in% pamID]) #include key PAM as a check
  
  # loop through each PAM for this gorilla:
  CBs[[c]]$w.pam <- list()
  for (pam in pamsToGet){ 
    
    #############################################################################################
    #############################################################################################
    #############################################################################################
    #############################################################################################
    
    # Get start times data.frame within the focal PAM based on time gaps in the key PAM:
    CBs[[c]]$w.pam <- mergedClipsFromTimeGaps(clip.start=first.clip.t.key, inPAMfile=first.clip.sound.path.key, gaps=time.gaps, getPAM=pam, clipLength=len, path=".", keyPAMstarts4check = dets.key.pam$startClip, tabulate = TRUE, add=0)
    
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
  }
}

saveRDS(CBs, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_CBs.rds"))



