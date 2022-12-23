# 2022-12-23 AndyP
# Explore_MEDuSA combine by id/run
# m_dir full path of directory for MEDuSA interpolated files
# /Users/andypapale/vmPFC/Explore_MEDuSA/interpolated/rt_aligned
# w_dir write directory to output concatenated data
# /Volumes/Users/Andrew/MEDuSA_data_Explore
# name_to_save - name of file to save
# 'fb-vmPFC'

load_MEDuSA_data_Explore <- function(m_dir,w_dir,name_to_save){
  
  
  library(tidyverse)
  
  mf <- list.files(path=m_dir,pattern='^sub',full.names=TRUE)
  nF <- length(mf)
  md <- NULL
  for (iF in 1:nF){
    m0 <- mf[iF]
    m1 <- sub("*.*run","",m0)
    m2 <- as.integer(sub("_interpolated.csv.gz","",m1))
    m4 <- substr(m0,77,77)
    message(paste('file ', iF, 'out of', nF))
    currF <- read_csv(mf[iF])
    currF <- currF %>% mutate(run=m2,id=m4)
    md <- rbind(md,currF)
  }
  
  save(md,file=paste0(w_dir,'/',name_to_save,'.Rdata'))
}