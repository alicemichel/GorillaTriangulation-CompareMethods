

#dat <- readRDS("~/Library/CloudStorage/Box-Box/Remaining/Students/Alice Michel/Wrapping test/20230206_UJ_times.rds")
dat <- readRDS("~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/20230206_UJ_times.rds")

# try it with this weird situation:
#dat$ind = c(rep("J",32), rep("U",5))

# which would look like this in a timeline:
plot(dat$full.start[dat$ind=="J"], rep(1, length(dat$full.start[dat$ind=="J"])), col="red", pch=16, xlim=c(min(dat$full.start), max(dat$full.start)))
points(dat$full.start[dat$ind=="U"], rep(1, length(dat$full.start[dat$ind=="U"])), pch=16)

metricCalc <- function(dat){##calculates all the U->J and J->U intervals, then averages them
	dat <- dat[order(dat$full.start),]
	x <- dat$ind=="J"
	JUs <- x[-1]-x[-length(x)]==-1 #adjacent intervals where the first is J and the next is U
	UJs <- x[-1]-x[-length(x)]==1 #adjacent intervals where the first is U and the next is J (0: #adjacent intervals where the first is J and the next is J or U and U)
	UJJUs <- c(dat$full.start[which(JUs)+1]-dat$full.start[JUs], #the +1 brings it back to starting at 1 #EXCLUDES THE LAST ONE BC IT IS U THEN U!
	           dat$full.start[which(UJs)+1]-dat$full.start[UJs]) #skips the same JJ UU switch at the end before the last unanswered J continues
  #only concern - double counting?
	return(mean(UJJUs))
}

randomization <- function(dat, fun = metricCalc, n=1000, shift="J"){
	obsMetric <- fun(dat)
	simulMetric <- numeric(0)
	for(i in 1:n){
	  cesure <- as.POSIXct(runif(1, min(dat$full.start), max(dat$full.start))) #pick a random time in the interval
	  #test <- as.POSIXct(sample(c(min(dat$full.start):max(dat$full.start)), 1)) same how I would do it - curious if this is doing the runif under the hood?
	  dat2 <- dat
	  dat2$full.start[dat2$ind==shift] <- dat2$full.start[dat2$ind==shift]+(cesure-min(dat$full.start)) #add each start time for the shifting ind to bring it to starting at the random time above
	  dat2$full.start[dat2$ind==shift & dat2$full.start>max(dat$full.start)] <- dat2$full.start[dat2$ind==shift & dat2$full.start>max(dat$full.start)]-diff(range(dat$full.start)) #if it goes over the last time for either, subtract the full range of time of CBs - ugh that works yes
	  dat2 <- dat2[order(dat2$full.start),] #reorder them so the ones that wrapped around go to the start so the function works
	  simulMetric <- c(simulMetric, as.numeric(fun(dat2), units="secs"))
	}	
	return(list(as.numeric(obsMetric), as.numeric(simulMetric)))
}

res <- randomization(dat, fun=metricCalc, n=1000, shift="J")
(pVal <- mean(res[[2]]<res[[1]]))##0.029 #==what fraction of the sims are less than the obs?

shift="J"
par(mfrow=c(1,1))
hist(res[[2]], main=paste(shift, "and", unique(dat$ind[dat$ind!=shift])), 
     xlab="Mean adjacent chest beat time gap (secs)", col="grey95", border="grey70")
abline(v=res[[1]], lty=4, col="red", lwd=2)
text(res[[1]]+(max(res[[2]])-min(res[[2]]))/10, 140, labels=paste("p = ", round(pVal, 3)), col="red")





