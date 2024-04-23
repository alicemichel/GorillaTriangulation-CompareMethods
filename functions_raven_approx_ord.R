

approxOrd <- function(Raven.selections.path, clipLength = 7, buffer = 0, cutoff = 0, path=".", pams=rownames(pam.xy), num=nrow(df)) {
  
  dets <- read.table(Raven.selections.path, header = TRUE, sep = "\t")
  dets <- dets[!(rownames(dets) %in% (0:cutoff)),]
  dets$start <- dets$File.Offset..s.
  dets$sound.files <- gsub(".*/","",dets$Begin.Path)
  dets$pam <- substr(dets$sound.files, start=1, stop=1)
  dets$duration <- dets$Sample.Length..samples./44100
  
  df <- dets[,c("pam", "sound.files", "start", "duration")]
  df$approxOrd <- NA
  df$ind <- NA
  df$min.cut <- NA
  df$ordered.cuts <- NA
  df$note <- NA
    
  for (i in 1:num){
    
    #w <- read_wave(df[i,]$sound.files, from = max(0, df[i,]$start - buffer), to = (df[i,]$start + clipLength))
    #viewSpec(w, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning", main = df[i,]$pam)
    
    fullTimeSecondsToGet <- fullTimeFromClipStart(sound.path = df[i,]$sound.files, clip.start = df[i,]$start)
    
    dateYYYYmmdd <- substr(df[i,]$sound.files, start = 4, stop = 11) #there will be an error if this differs between the PAMs!
    
    onOthers <- data.frame(pams)
    onOthers$sound.files <- NA
    onOthers$start <- NA
    errcnt <- 0
    for (getPAM in pams){
      tryCatch({
      tmp <- whichFile(path, fullTimeSecondsToGet, getPAM, dateYYYYmmdd)
      onOthers[onOthers$pams==getPAM,]$sound.files <- tmp[[1]]
      onOthers[onOthers$pams==getPAM,]$start <- tmp[[2]]
      }, error=function(e){0}) #cat("ERROR :",conditionMessage(e), "\n")
    }
    rownames(onOthers) <- onOthers$pams
    # reorder the df so that the detection PAM comes first, and remove missing data (PAM not in folder, this should not be any later on, but it might be if the time isn't found, actually you'd probably get an error then in the function whichFile)
    onOthers1 <- as.data.frame(rbind(onOthers[rownames(onOthers)==df[i,]$pam,-1], onOthers[!is.na(onOthers$start) & rownames(onOthers)!=df[i,]$pam,-1]))
    cat(paste0("removed"), sum(is.na(onOthers$start)), "PAMs without soundfiles to search \n  printing soundclips to select start times next for selection", i, "out of", nrow(df), "..")
    # now need to print these and select in them to see who is first! and loop through all selections...
    
    onOthers1$cut <- NA
    par(mfrow=c(2,1))
    for (j in rownames(onOthers1)){
      wav <- readWave(onOthers1[j,]$sound.files, header=TRUE)
      Hz.file.length <- wav$samples / wav$sample.rate
      tmpwv <- read_wave(onOthers1[j,]$sound.files, from = max(0, onOthers1[j,]$start - buffer), to = min(onOthers1[j,]$start + clipLength,Hz.file.length))
      viewSpec(tmpwv, interactive = F, units="seconds", frq.lim = c(0.1, 0.9), wl=2048, ovlp = 90, wn="hanning", main = rownames(onOthers1[j,]))
      
      usr <- readline("start of chest beat:")
      if (usr=="skip"){
        break
      }
      onOthers1[j,]$cut <- as.numeric(usr)
    }
    ordered <- onOthers1[order(onOthers1$cut),]
    
    df[i,]$approxOrd <- paste0(rownames(ordered), collapse = "")
    df[i,]$ind <- substr(df[i,]$approxOrd, start=1, stop=1)
    df[i,]$min.cut <- ordered[1,]$cut
    df[i,]$ordered.cuts <- paste(ordered$cut, collapse = "_")
    df[i,]$note <- readline("ENTER NOTES HERE:")
    
    #write to csv as you go:
    write.csv(df, paste0(path, "/ordering.csv"))
    
  }
  return(df)
}


check_spectro <- function(df, rowid, buffer=0, clipLength=4){
  for (i in rowid){
    w <- read_wave(df[rownames(df)==i,]$sound.files, from = max(0,df[rownames(df)==i,]$start - buffer), to = df[rownames(df)==i,]$start + clipLength)
    viewSpec(w, interactive = F, units="seconds", frq.lim = c(0.1, 0.9), wl=2048, ovlp = 90, wn="hanning", main = paste(i, "\n", df[rownames(df)==i,]$ordered.cuts, "\n", df[rownames(df)==i,]$approxOrd))
  }
}


