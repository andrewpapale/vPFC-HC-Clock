# 2025-08-01 AndyP
# Get N Peaks and Peak Prominence from SCEPTIC

library(tidyverse)
library(pracma)
library(breath)

mmc <- read_csv('/Volumes/Users/Andrew/mmclock_fmri_decay_factorize_selective_psequate_mfx_sceptic_trial_outputs_by_timestep.csv')
exp <- read_csv('/Volumes/Users/Andrew/explore_decay_factorize_selective_psequate_fixedparams_fmri_mfx_sceptic_trial_outputs_by_timestep.csv')

v <- mmc

nR = nrow(v)

y0 <- NULL
trial <- NULL
peaks <- NULL
id <- NULL
for (iR in 1:nR){
  y <- cbind(v$V_1[iR],v$V_2[iR],v$V_3[iR],v$V_4[iR],v$V_5[iR],v$V_6[iR],v$V_7[iR],v$V_8[iR],v$V_9[iR],v$V_10[iR]
             ,v$V_11[iR],v$V_12[iR],v$V_13[iR],v$V_14[iR],v$V_15[iR],v$V_16[iR],v$V_17[iR],v$V_18[iR],v$V_19[iR],v$V_20[iR]
             ,v$V_21[iR],v$V_22[iR],v$V_23[iR],v$V_24[iR])
  #y0 <- rbind(y0,y)
  peaks1 <- pracma::findpeaks(as.numeric(y),minpeakheight = 0.3, minpeakdistance = 4)
  nP <- nrow(peaks1)
  if (!is.null(peaks1)){
    P <- NULL
    for (iP in 1:nP){
      vL <- y[peaks1[iP,3]]
      vR <- y[peaks1[iP,4]]
      P[iP] <- min(peaks[iP,1]-vL,peaks[iP,1]-vR)
    }
    peaks0 <- data.frame(ix = peaks1[,2],peak=peaks1[,1],prominence = P)
    peaks0 <- peaks0 %>% filter(prominence > 0.33)
    
    trial0 <- as.matrix(rep(v$asc_trial[iR],nrow(peaks0)))
    id0 <- as.matrix(rep(v$id[iR],nrow(peaks0)))
    
    trial <- rbind(trial,trial0)  
    peaks <- rbind(peaks,peaks0)
    id <- rbind(id,id0)
  }
}


# iR = sample.int(nR,1,replace=FALSE)
# y <- cbind(v$V_1[iR],v$V_2[iR],v$V_3[iR],v$V_4[iR],v$V_5[iR],v$V_6[iR],v$V_7[iR],v$V_8[iR],v$V_9[iR],v$V_10[iR]
#            ,v$V_11[iR],v$V_12[iR],v$V_13[iR],v$V_14[iR],v$V_15[iR],v$V_16[iR],v$V_17[iR],v$V_18[iR],v$V_19[iR],v$V_20[iR]
#            ,v$V_21[iR],v$V_22[iR],v$V_23[iR],v$V_24[iR])
# y1 <- data.frame(y=t(y)) 
# peaks <- pracma::findpeaks(as.numeric(y),minpeakheight = 0.3, minpeakdistance = 4)
# nP <- nrow(peaks)
# P <- NULL
# for (iP in 1:nP){
#   vL <- y[peaks[iP,3]]
#   vR <- y[peaks[iP,4]]
#   P[iP] <- min(peaks[iP,1]-vL,peaks[iP,1]-vR)
# }
# peaks <- data.frame(ix = peaks[,2],peak=peaks[,1],prominence = P)
# peaks <- peaks %>% filter(prominence > 0.33)
# if (!isempty(peaks)){
#   ggplot(y1,aes(x=sort(as.numeric(rownames(y1))),y=y)) + geom_point() + 
#     geom_line() + geom_point(data=peaks,aes(x=ix,y=peak,size=prominence*3))
# } else {
#   ggplot(y1,aes(x=sort(as.numeric(rownames(y1))),y=y)) + geom_point() + 
#     geom_line()
# }

