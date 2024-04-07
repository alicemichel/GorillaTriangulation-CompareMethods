





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





working <- function(x){
  
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
  
}


## In this example, we get the lag J-Q. This was not a combo in the field that led to good triangulation, but let's see if we can get a lag from the ordered cross-correlation method. Compare it to the expected lag based on the true field-validated location.

q <- read_wave("20230206_192154__Q_clip_start_T1492.049489.wav") 
len = 4
exten <- 2 #add 2 secs to beginning of one/both
stq = (length(q)/44100-4)
(gap = stq - len)
q1 <- cutWave(q, from=0, to=len)
q2 <- cutWave(q, from=stq, to=stq+len)
q <- bind(q1, q2)
viewSpec(q, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning")
st = 1492.049489
q.full <- read_wave("Q_S20230206T185702.506263+0100_E20230206T195657.382118+0100_+1.26056+17.20221.wav", from = st, to = st + len)
identical(q.full, q1)
q.ext <- bind(read_wave("Q_S20230206T185702.506263+0100_E20230206T195657.382118+0100_+1.26056+17.20221.wav", from = st - exten, to = st), q, read_wave("Q_S20230206T185702.506263+0100_E20230206T195657.382118+0100_+1.26056+17.20221.wav", from = st + len, to = st + len + exten))
# viewSpec(q.ext, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning")

j <- read_wave("20230206_192154__J_clip_start_T1479.476155.wav")
stj = (length(j)/44100-4)
j1 <- cutWave(j, from=0, to=len)
j2 <- cutWave(j, from=stq, to=stq+len)
j <- bind(j1, j2)
viewSpec(j, interactive = F, frq.lim = c(0.1, 0.9), wl=2048, ovlp=99, wn="hanning")

q.fp <- file.path(tempdir(), "2023-02-06_192154_UTC_Q.wav")
writeWave(q.ext, q.fp)

j.fp <- file.path(tempdir(), "2023-02-06_192154_UTC_J.wav")
writeWave(j, j.fp)


# iteratively make each of them the template for the other (should be the same result though, so randomly pick alternatively)
q.temp <- makeCorTemplate(q.fp, frq.lim = c(0.1, 0.9), wl = 2048, ovlp = 99, name = "q")
j.temp <- makeCorTemplate(j.fp, frq.lim = c(0.1, 0.9), wl = 2048, ovlp = 99, name = "j")
# 1/10 of all points

ctemps <- combineCorTemplates(q.temp, j.temp)
templateCutoff(ctemps)[1:2] <- 0.3
ctemps

(cscores <- corMatch(q.fp, ctemps))

(cdetects <- findPeaks(cscores[2])) #2nd to get where q (being longer) fits to j (shorter template)

getDetections(cdetects)[,3]==cdetects@peaks[[1]][,2]
mid.exp <- length(q.ext)/44100/2
mid.cor <- cdetects@peaks[[1]][,2]
(lag <- mid.cor - mid.exp) # + = before template
# this should be the lag that's good to go...

