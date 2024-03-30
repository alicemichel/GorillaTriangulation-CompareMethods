

hz2GPStime <- function(clipStart, soundpath, reverse=FALSE) {
  
  if (length(clipStart) == length(soundpath)){
    n = length(clipStart)
  }
  
  else{
    stop("The number of clips and files do not match.")
  }
  
  GPS.start.clips <- c()
  
  for (i in 1:n){
    
    pathStart <- fullTimeFromClipStart(soundpath[i], 0)
    
    endFileTxt <- substr(gsub(".*/","", soundpath[i]), start = 42, stop = 54)
    pathEnd <- 3600*as.numeric(substr(endFileTxt, start = 1, stop = 2))+60*as.numeric(substr(endFileTxt, start = 3, stop = 4))+as.numeric(substr(endFileTxt, start = 5, stop = 15))
    
    GPS.file.length <- pathEnd - pathStart
    
    wav <- readWave(soundpath[i], header=TRUE)
    Hz.file.length <- wav$samples / wav$sample.rate
    
    GPS.Hz.ratio <- GPS.file.length / Hz.file.length
    
    GPS.start.clip <- GPS.Hz.ratio * clipStart[i]
    # Check the worst by setting clipStart to Hz.file.length
    # --> worst it could be off is by 0.087seconds if at the very end of the clip
    
    if (reverse==TRUE){
      GPS.start.clip <- Hz.file.length / GPS.file.length * clipStart[i]
    }
    
    GPS.start.clips <- c(GPS.start.clips, GPS.start.clip)
  }
  return(GPS.start.clips)
}
