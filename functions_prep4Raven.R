
ravenprep <- function(CBdetsList, clipLength, savewvs=TRUE, detbegin=rep(0, length(CBdetsList)), detlim=rep(NA, length(CBdetsList))){ #change var names so they don't match what's in code
  subDir = paste0(gsub("-",".",Sys.Date()), "_", substr(CBdetsList[[1]]$detections$sound.files[1], start=4, stop=11), "_individuals_", paste0(names(CBdetsList), collapse="_"))
  if (!file.exists(subDir)){
    dir.create(file.path(subDir))
  }
  wavlist <- list()
  for (k in inds){ 
    ki=which(inds %in% k)
    detlim[ki] <- ifelse(is.na(detlim[ki]), sum(!is.na(CBdetsList[[k]][[4]]$soundfile)), detlim[ki])
    CBdetsList[[k]][-(1:3)] = lapply(CBdetsList[[k]][-(1:3)], FUN = function(x) x[1:detlim[ki],])
    CBdetsList[[k]][-(1:3)] = lapply(CBdetsList[[k]][-(1:3)], FUN = function(x) x[detbegin[ki]:nrow(x),])
    wavlist[[k]] <- lapply(CBdetsList[[k]][-(1:3)], FUN = function(x) {df2mergedWav(x, clipLength)})
    if (length(unique(lapply(wavlist[[k]], FUN = function(x) duration(x))))!=1){
      stop(paste("ERROR in", k, " -- The files are not the same length! \n Consider reducing the length of each clip or implementing a limit on the first or last detections (former is less recommended due to initial buffer)."))
    }
    for (w in 1:length(wavlist[[k]])){
      writeWave(wavlist[[k]][[w]], paste0(subDir,"/", names(CBdetsList)[ki], "_", substr(names(wavlist[[k]][w]), start=3, stop=3), ".wav"))
    }
  }
  if (savewvs==TRUE){
    timetxt=paste0("fileCreated_",substr(Sys.time(), start=12, stop=19))
    timetxt=gsub(":","",timetxt)
    flnm <- paste0(subDir, "/", substr(CBdetsList[[1]]$detections$sound.files[1], start=4, stop=11), "_RavenLags_", timetxt, ".xlsx")
    pamset=lapply(wavlist, function(x) substr(names(x), start=3, stop=3))
    writexl::write_xlsx(data.frame(IndID=rep(names(CBdetsList), lengths(pamset)), PAM=Reduce(c, pamset), lag=NA), flnm)
  }
  return(list(indWavs=wavlist, indDets=CBdetsList))
}



## Problem: overlapping into next file... what if the sound is too close to the end of the file that the length doesn't cover it!
## Solutions: check it! And shift length if that is still OK.
## Or, if it's too tight, manually cut off the last few CBs that are at the end. Or remove the one ID. In the J example of 0206 it should be #23

