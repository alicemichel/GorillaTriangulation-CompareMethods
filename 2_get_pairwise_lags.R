
library(monitoR)
ogdir <- "~/Library/CloudStorage/Box-Box/AliceMichel/Research/Lac Tele/FieldSeason2/00 Analysis/Office Triangulation/20230206_new_idea/"
setwd(ogdir)

# To summarize, for template matching using spectrogram cross-correlation one would typically use the
# following functions in the order given:
#
# 1. makeCorTemplate to make the template(s)
# 2. combineCorTemplates to combine templates together in a single template list
# 3. corMatch to calculate scores
# 4. findPeaks to find peaks and identify detections
# 5. plot to see the scores and detections
# 6. getDetections to get the (numeric) detection results
# 7. templateCutoff to change template cutoffs in the detection list object (iteratively with plot and getDetections)


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

