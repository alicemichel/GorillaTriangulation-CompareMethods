
fullTimeFromClipStart <- function(sound.path, clip.start){
  startFileTxt <- substr(gsub(".*/","", sound.path), start = 13, stop = 25)
  startFileTime <- 3600*as.numeric(substr(startFileTxt, start = 1, stop = 2))+60*as.numeric(substr(startFileTxt, start = 3, stop = 4))+as.numeric(substr(startFileTxt, start = 5, stop = 15))
  return(startFileTime+clip.start)
}

pullTime <- function(clip.start, inPAMfile, getPAM, clipLength=len, path=".", plot=TRUE, clockfix=TRUE){
  fullTimeSecondsToGet <- fullTimeFromClipStart(sound.path = inPAMfile, clip.start = clip.start)
  dateYYYYmmdd <- substr(inPAMfile, start = 4, stop = 11)
  return(pullFromPam(fullTimeSecondsToGet, dateYYYYmmdd, getPAM, clipLength, path, plot, clockfix))
}

whichFile <- function(path, fullTimeSecondsToGet, getPAM, dateYYYYmmdd){
  soundfiles <- list.files(path, pattern = "S20.*.wav")
  kfls <- soundfiles[grep(paste0(getPAM,"_"), soundfiles)]
  daytext = substr(kfls, start = 4, stop = 11) 
  kfls = kfls[daytext == dateYYYYmmdd] #there will be an error if this differs between the PAMs!
  kfls = kfls[nchar(kfls)==81] # exclude cut files from processing.
  # cat("note: files must have the format letter_file!")
  
  t1txt <- substr(kfls, start = 13, stop = 25) 
  t1s <- 3600*as.numeric(substr(t1txt, start = 1, stop = 2))+60*as.numeric(substr(t1txt, start = 3, stop = 4))+as.numeric(substr(t1txt, start = 5, stop = 15))
  
  difs <- fullTimeSecondsToGet - t1s
  m <- which(difs==min(difs[difs>0])) #in the file not the end of the one before
  
  if (T %in% unique(difs>0)){
    kfl <- kfls[m]
    bts <- fullTimeSecondsToGet - t1s[m]
  }
  if (bts>3594.82){
    m2 <- which(abs(difs)==min(abs(difs[difs!=difs[m]])))
    kfl <- kfls[m2]
    bts <- max(0, fullTimeSecondsToGet - t1s[m2])
  }
  
  tmp <- data.frame(kfl, bts)

  return(tmp)
}

pullFromPam <- function(fullTimeSecondsToGet, dateYYYYmmdd, getPAM, clipLength, path=".", plot=TRUE, clockfix=FALSE){
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  tmp <- whichFile(path, fullTimeSecondsToGet, getPAM, dateYYYYmmdd)
  kfl <- tmp[1,1]
  bts <- tmp[1,2]
  
  ## Clock drift correction NOT DONE if FALSE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (clockfix==TRUE){
    bts.clock.redrift <- hz2GPStime(bts, kfl, reverse = TRUE)
  }
  if (clockfix==FALSE){
    bts.clock.redrift <- bts
  }
  
  tmpls <- list(getPAM, kfl, bts.clock.redrift, read_wave(kfl, from = bts.clock.redrift, to = bts.clock.redrift + clipLength))
  if (plot==TRUE){
    viewSpec(tmpls[[4]], main = pam, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning")
  }
  return(tmpls)
  
  options(warn = oldw)
}


mergedClipsFromTimeGaps <- function(clip.start, inPAMfile, gaps, getPAM, clipLength=len, path=".", tabulate=FALSE, keyPAMstarts4check, add=10){
  
  pam.1st.clip <- pullTime(clip.start, inPAMfile, getPAM, clipLength, plot=FALSE, clockfix=FALSE) #FALSE=GPS time
  first.file.start <- fullTimeFromClipStart(pam.1st.clip[[2]],0)
  first.clip.start <- fullTimeFromClipStart(pam.1st.clip[[2]],pam.1st.clip[[3]])
  fl.1st.clip <- which(list.files(path, pattern = "S20.*.wav")[grep(paste0(pam.1st.clip[[1]],"_"), list.files(path, pattern = "S20.*.wav"))] %in% pam.1st.clip[[2]])
  
  all.clip.starts <- first.clip.start + gaps
  starts.in.fls <- all.clip.starts - first.file.start
  
  dateYYYYmmdd <- substr(inPAMfile, start = 4, stop = 11)
  
  ## THIS HAD AN ERROR: ENDED UP AS 1 WHEN SHOULD BE 2 IF THE FIRST CLIP IS IN THE SECOND FILE...
  ## IN TESTING!
  rems <- (starts.in.fls-3594.82)%/%3594.82 + 2 + (fl.1st.clip - 1) #TO ACCOUNT FOR STARTING AFTER THE FIRST ONE...
  #############################################################################################
  
  soundfiles=list.files(path, pattern = "S20.*.wav")
  kfls <- soundfiles[grep(paste0(getPAM,"_"), soundfiles)]
  daytext = substr(kfls, start = 4, stop = 11) 
  kfls = kfls[daytext == dateYYYYmmdd] # exclude other days -- could be a problem for midnight to 1am
  kfls = kfls[nchar(kfls)==81] # exclude cut files from processing.

  fls <- kfls[rems]
  all.clip.starts1 <- all.clip.starts - fullTimeFromClipStart(fls,0)
  
  ## Clock drift correction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Here, go from GPS "true" time back to Hz-based "raw" time within each PAM:
  all.clip.starts.re.drift <- hz2GPStime(clipStart = all.clip.starts1, soundpath = fls, reverse = TRUE)
  
  ## Check that the matching PAM matches:
  if (substr(inPAMfile, start = 1, stop = 1) == pam) {
    if(unique(round(keyPAMstarts4check,8)!=round(all.clip.starts.re.drift,8))){
      stop("There is a problem with the clock drift reversal on the key PAM.") 
    }
  }
  
  # Save data frame of soundfiles and start times
  if (tabulate==TRUE){
    newdets <- data.frame(selection = 1:length(fls), soundfile = fls, startClip = all.clip.starts.re.drift)
    if (add>0){
      extend <- newdets[1,]$startClip - add
      slctn <- 0
      snd <- newdets[1,]$soundfile
      if (extend < 0 ){
        extend <- newdets[nrow(newdets),]$startClip + clipLength + add
        slctn <- nrow(newdets) + 1
        snd <- newdets[nrow(newdets),]$soundfile
      }
    newdets <- rbind(newdets, c(slctn, snd, extend))
    }
    newdets$startClip <- as.numeric(newdets$startClip)
    newdets$selection <- as.numeric(newdets$selection)
    newdets <- newdets[order(newdets$selection),]
    return(newdets)
  }
  
  # Save merged wave sound clips
  else {
    clps <- list()
    for (i in 1:length(gaps)){
      clps[[i]] <- read_wave(fls[i], from = all.clip.starts.re.drift[i], to = (all.clip.starts.re.drift[i] + clipLength))
    }
    # Bind them all up into one:
    return(do.call(tuneR::bind, args = clps))
  }
  
}



df2mergedWav <- function(df, clipLength){
  clps <- list()
  for (i in 1:nrow(df)){
    clps[[i]] <- read_wave(df[i,]$soundfile, from = df[i,]$startClip, to = df[i,]$startClip + clipLength)
  }
  # Bind them all up into one:
  return(do.call(tuneR::bind, args = clps))
}







