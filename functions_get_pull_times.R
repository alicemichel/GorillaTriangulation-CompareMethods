
fullTimeFromClipStart <- function(sound.path, clip.start){
  startFileTxt <- substr(gsub(".*/","", sound.path), start = 13, stop = 25)
  startFileTime <- 3600*as.numeric(substr(startFileTxt, start = 1, stop = 2))+60*as.numeric(substr(startFileTxt, start = 3, stop = 4))+as.numeric(substr(startFileTxt, start = 5, stop = 15))
  return(startFileTime+clip.start)
}

pullTime <- function(clip.start, inPAMfile, getPAM, clipLength=len){
  fullTimeSecondsToGet <- fullTimeFromClipStart(sound.path = inPAMfile, clip.start = clip.start)
  dateYYYYmmdd <- substr(inPAMfile, start = 4, stop = 11)
  return(pullFromPam(fullTimeSecondsToGet, dateYYYYmmdd, getPAM=pam, clipLength=len, path="."))
}


pullFromPam <- function(fullTimeSecondsToGet, dateYYYYmmdd, getPAM=pam, clipLength=len, path="."){
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  soundfiles <- list.files(path, pattern = "S20.*.wav")
  
  kfls <- soundfiles[grep(paste0(getPAM,"_"), soundfiles)]
  daytext = substr(kfls, start = 4, stop = 11) 
  kfls = kfls[daytext == dateYYYYmmdd] # exclude other days -- could be a problem for midnight to 1am
  kfls = kfls[nchar(kfls)==81] # exclude cut files from processing.
  cat("note: files must have the format letter_file!")
  
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
  
  tmpls <- list(getPAM, kfl, bts, read_wave(kfl, from = bts, to = bts + len))
  viewSpec(tmpls[[4]], main = pam, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning")
  return(tmpls)
  
  options(warn = oldw)
}


mergedClipsFromTimeGaps <- function(clip.start, inPAMfile, gaps, getPAM, clipLength=len, path="."){
  
  pam.1st.clip <- pullTime(clip.start, inPAMfile, getPAM, clipLength)
  first.file.start <- fullTimeFromClipStart(pam.1st.clip[[2]],0)
  first.clip.start <- fullTimeFromClipStart(pam.1st.clip[[2]],pam.1st.clip[[3]])
  all.clip.starts <- first.clip.start + gaps
  starts.in.fls <- all.clip.starts - first.file.start
  
  rems <- (starts.in.fls-3594.82)%/%3594.82 + 2
  
  list.files(path, pattern = "S20.*.wav")
  kfls <- soundfiles[grep(paste0(getPAM,"_"), soundfiles)]
  daytext = substr(kfls, start = 4, stop = 11) 
  kfls = kfls[daytext == dateYYYYmmdd] # exclude other days -- could be a problem for midnight to 1am
  kfls = kfls[nchar(kfls)==81] # exclude cut files from processing.
  cat("note: files must have the format letter_file!")
  
  fls <- kfls[rems]
  all.clip.starts1 <- all.clip.starts - fullTimeFromClipStart(fls,0)
  
  clps <- list()
  for (i in 1:length(gaps)){
    clps[[i]] <- read_wave(fls[i], from = all.clip.starts1[i], to = (all.clip.starts1[i] + len))
  }
  return(do.call(tuneR::bind, args = clps))
}




