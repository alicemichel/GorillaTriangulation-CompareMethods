
ravenprep <- function(CBdetsList, clipLength, savewvs=TRUE, detbegin=rep(0, length(CBs)), detlim=rep(NA, length(CBs))){ #change var names so they don't match what's in code
  subDir = paste0(Sys.Date(), "_", substr(CBdetsList[[1]]$detections$sound.files[1], start=4, stop=11), "_individuals_", paste0(names(CBdetsList), collapse="_"))
  if (!file.exists(subDir)){
    dir.create(file.path(subDir))
  }
  wavlist <- list()
  for (k in inds){ 
    detlim <- ifelse(is.na(detlim), nrow(CBdetsList[[k]][[4]]), detlim)
    CBdetsList[[k]][-(1:3)] = lapply(CBdetsList[[k]][-(1:3)], FUN = function(x) x[1:detlim[k],])
    CBdetsList[[k]][-(1:3)] = lapply(CBdetsList[[k]][-(1:3)], FUN = function(x) x[detbegin[k]:nrow(x),])
    wavlist[[k]] <- lapply(CBdetsList[[k]][-(1:3)], FUN = function(x) {df2mergedWav(x, clipLength)})
    if (length(unique(lapply(wavlist[[k]], FUN = function(x) duration(x))))!=1){
      stop("The files are not the same length! \n Consider reducing the length of each clip or implementing a limit on the first or last detections (former is less recommended due to initial buffer).")
    }
    lapply(wavlist[[k]], FUN = function(x) writeWave(x, paste0(subDir,"/", names(CBdetsList)[k], "_", substr(names(x), start=3, stop=3), ".wav")))
  }
  
  if (savewvs==TRUE){
    flnm <- paste0(subDir, "/_Raven_lags_", substr(Sys.time(), start=12, stop=19), ".csv")
    write.csv(data.frame(IndID=names(CBdetsList), PAM=substr(names(wavlist[[k]]), start=3, stop=3), lag=NA), flnm)
  }
  return(wavlist)
}

ravenprep(CBs, len, detlim = 23)


## Problem: overlapping into next file... what if the sound is too close to the end of the file that the length doesn't cover it!
## Solutions: check it! And shift length if that is still OK.
## Or, if it's too tight, manually cut off the last few CBs that are at the end. Or remove the one ID. In the J example of 0206 it should be #23



xcorrlags <- function(CBs, len, frq.lim = c(0.2, 0.7), wl = 2048, ovlp = 99){
  wavlist <- list()
  lags <- data.frame()
  for (k in inds){ 
    wavlist[[k]] <- lapply((CBs[[k]])[-(1:3)], FUN = function(x) {df2mergedWav(x, len)})
    
    ## Grab the key PAM templates:
    key.tmp.path <- file.path(tempdir(), "tmp.key.wav")
    writeWave(CBs[[k]]$w.keypam, key.tmp.path)
    template <- makeCorTemplate(key.tmp.path, frq.lim = frq.lim, wl = wl, ovlp = ovlp, name = "Key")
    templateCutoff(template)[1] <- 0.3
    
    ## Loop thru the PAMs and do correlations:
    for (p in 1:length((wavlist[[k]]))){
      
      extended = duration(wavlist[[k]][[p]]) - duration(CBs[[k]]$w.keypam)
      
      fp <- file.path(".", "tmp.wav")
      writeWave(wavlist[[k]][[3]], fp)
      (cscores <- corMatch(fp, template)) #ok if template length == file length?
      
      (cdetects <- findPeaks(cscores))
      
      getDetections(cdetects)[,3]==cdetects@peaks[[1]][,2]
      mid.exp <- extended + (length(CBs[[k]]$w.keypam)/44100)/2
      mid.cor <- cdetects@peaks[[1]][,2]
      lag <- mid.cor - mid.exp # + = before template
      pamname=substr(names(wavlist[[k]])[p], start=3, stop=3)
      cat(paste("Min lag from", k, "to", pamname, "is:", round(min(abs(lag)),3), "seconds"))
      # this should be the lag that's good to go...
      lags <- rbind(lags, c(k, pamname, lag[abs(lag)==min(abs(lag))]))
    }
  }
  names(lags) <- c("key.pam","pam", "lag")
  return(lags)
}







# Get samplingrate - assume delayed signal has same samplingrate
sr = j@samp.rate

# Extract signal - assume mono with signal in left
orig = j2@left
delayed = q2@left

# Zero pad, at least half. nextn() selects "factor rich" length
nlength = nextn(max(length(orig), length(delayed))*2)
orig = c(orig, rep(0, nlength-length(orig)))
delayed = c(delayed, rep(0, nlength-length(delayed)))

# Calculate cross correlation
ccor = fft(fft(orig)*Conj(fft(delayed)), inverse=TRUE) #Cong = complex conjugate for complex values #You find the complex conjugate simply by changing the sign of the imaginary part of the complex number. To find the complex conjugate of 4+7i we change the sign of the imaginary part. Thus the complex conjugate of 4+7i is 4 - 7i
(lmax = abs(which.max(Re(ccor)) - length(orig))) #Re = real part
(delay = lmax/sr)


## need to constrain frequency range


# Calculate cross correlation, try for lags up to 1000ms
ccor = ccf(orig,delayed, plot=FALSE, lag.max=sr)
cor = ccor$acf[,,1]
lag = ccor$lag[,,1]
lagmax = lag[which.max(cor)]
(delay = abs(lagmax/sr))
