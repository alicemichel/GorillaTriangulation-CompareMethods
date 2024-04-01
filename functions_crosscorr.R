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
      
      fp <- file.path(tempdir(), "tmp.wav")
      writeWave(wavlist[[k]][[p]], fp)
      (cscores <- corMatch(fp, template)) #ok if template length == file length?
      
      (cdetects <- findPeaks(cscores))
      
      getDetections(cdetects)[,3]==cdetects@peaks[[1]][,2]
      mid.exp <- extended + (length(CBs[[k]]$w.keypam)/44100)/2
      mid.cor <- cdetects@peaks[[1]][,2]
      lag <- mid.cor - mid.exp # + = before template
      pamname=substr(names(wavlist[[k]])[p], start=3, stop=3)
      cat(paste("Min lag from", k, "to", pamname, "is:", round(min(abs(lag)),3), "seconds"))
      # this should be the lag that's good to go...
      lags <- rbind(lags, c(k, pamname, min(lag)))
    }
  }
  return(lags)
}