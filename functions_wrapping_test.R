
metricCalc <- function(dat, shift){##calculates all the U->J and J->U intervals, then averages them
  dat <- dat[order(dat$full.start),]
  x <- dat$ind==shift
  JUs <- x[-1]-x[-length(x)]==-1 #adjacent intervals where the first is J and the next is U
  UJs <- x[-1]-x[-length(x)]==1 #adjacent intervals where the first is U and the next is J (0: #adjacent intervals where the first is J and the next is J or U and U)
  UJJUs <- c(dat$full.start[which(JUs)+1]-dat$full.start[JUs], #the +1 brings it back to starting at 1 #EXCLUDES THE LAST ONE BC IT IS U THEN U!
             dat$full.start[which(UJs)+1]-dat$full.start[UJs]) #skips the same JJ UU switch at the end before the last unanswered J continues
  #only concern - double counting?
  return(mean(UJJUs))
}

randomization <- function(dat, fun, n, shift){
  obsMetric <- fun(dat, shift)
  simulMetric <- numeric(0)
  for(i in 1:n){
    cesure <- as.POSIXct(runif(1, min(dat$full.start), max(dat$full.start))) #pick a random time in the interval
    dat2 <- dat
    dat2$full.start[dat2$ind==shift] <- dat2$full.start[dat2$ind==shift]+(cesure-min(dat$full.start)) #add each start time for the shifting ind to bring it to starting at the random time above
    dat2$full.start[dat2$ind==shift & dat2$full.start>max(dat$full.start)] <- dat2$full.start[dat2$ind==shift & dat2$full.start>max(dat$full.start)]-diff(range(dat$full.start)) #if it goes over the last time for either, subtract the full range of time of CBs - ugh that works yes
    dat2 <- dat2[order(dat2$full.start),] #reorder them so the ones that wrapped around go to the start so the function works
    simulMetric <- c(simulMetric, as.numeric(fun(dat2, shift), units="secs"))
  }	
  return(list(as.numeric(obsMetric), as.numeric(simulMetric)))
}

wrap <- function(dat, shift=ind, fun=metricCalc, n=1000){
  res <- randomization(dat, fun=metricCalc, n=1000, shift=ind)
  pVal <- mean(res[[2]]<res[[1]])
  par(mfrow=c(1,1))
  hist(res[[2]], main=paste(shift, "and", unique(dat$ind[dat$ind!=shift])), 
       xlab="Mean adjacent chest beat time gap (secs)", col="grey95", border="grey70")
  abline(v=res[[1]], lty=4, col="red", lwd=2)
  text(res[[1]]+(max(res[[2]])-min(res[[2]]))/10, 140, labels=paste("p = ", round(pVal, 3)), col="red")
return(pVal)
  }
