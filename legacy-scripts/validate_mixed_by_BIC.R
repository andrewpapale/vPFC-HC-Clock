# 2022-03-29 AndyP
# Test BIC of Uncompressed - Compressed Models

library(stringr)
library(pracma)
library(wesanderson)

plot_BIC_diff <- function(U_file='',C_file='',struct='',toalign='',toprocess=''){
  load(U_file)
  Ufit <- ddf$fit_df
  rm(ddf)
  load(C_file)
  Cfit <- ddf$fit_df
  rm(ddf)
  epoch_label = paste("Time relative to",toalign, "[s]")
  
  if (strcmp(struct,'HC')){
    pal1 <- palette(c('#F21A00','#3B9AB2'))
    Ufit <- Ufit %>% select(evt_time,HC_region,BIC) %>% rename(BICu=BIC)
    Cfit <- Cfit %>% select(evt_time,HC_region,BIC) %>% rename(BICc=BIC)
    D <- inner_join(Cfit,Ufit,by=c('evt_time','HC_region'))
    s <- str_split(C_file,'-')
    if (length(s[[1]])==9){
      LorR <-substr(s[[1]][9],1,1)
      totest <- paste0(s[[1]][6],'-',s[[1]][7])
    } else if (length(s[[1]])==8){
      LorR <- substr(s[[1]][8],1,1)
      totest <- s[[1]][6]
    } else {
    }
    pdf(paste0(s[[1]][5],'-',totest,'-',toprocess,'-',LorR,'-BICuminusc','.pdf'),width=9,height=3.5)
    gg <- ggplot(D, aes(x=evt_time,y=BICu-BICc,group=HC_region,color=HC_region)) + 
      geom_line(size=1) + scale_color_manual(values = pal1,labels=c('AH','PH')) + 
      xlab(epoch_label) + ylab('BIC: Uncompressed - Compressed')
    print(gg)
    dev.off()
  } else if (strcmp(struct,'vmPFC')) {
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    Ufit <- Ufit %>% select(evt_time,network,BIC) %>% rename(BICu=BIC)
    Cfit <- Cfit %>% select(evt_time,network,BIC) %>% rename(BICc=BIC)
    D <- inner_join(Cfit,Ufit,by=c('evt_time','network'))
    s <- str_split(C_file,'-')
    if (length(s[[1]])==8){
      totest <- paste0(s[[1]][6],'-',s[[1]][7])
    } else if (length(s[[1]])==7){
      totest <- s[[1]][6]
    } else {
    }
    pdf(paste0(s[[1]][5],'-',totest,'-',toprocess,'-BICuminusc','.pdf'),width=9,height=3.5)
    gg <- ggplot(D, aes(x=evt_time,y=BICu-BICc,group=network,color=network)) + 
      geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
      xlab(epoch_label) + ylab('BIC: Uncompressed - Compressed')
    print(gg)
    dev.off()
  } else if (strcmp(struct,'vmPFC-HC')){
    pal = wes_palette("FantasticFox1", 3, type = "discrete")
    Ufit <- Ufit %>% select(evt_time,network,HC_region,BIC) %>% rename(BICu=BIC)
    Cfit <- Cfit %>% select(evt_time,network,HC_region,BIC) %>% rename(BICc=BIC)
    D <- inner_join(Cfit,Ufit,by=c('evt_time','network','HC_region'))
    s <- str_split(C_file,'-')
    if (length(s[[1]])==11){
      LorR <-substr(s[[1]][11],1,1)
      totest <- paste0(s[[1]][6],'-',s[[1]][7])
    } else if (length(s[[1]])==10){
      totest <- s[[1]][6]
      LorR <-substr(s[[1]][10],1,1)
    } else{
      browser()
    }
    pdf(paste0(s[[1]][5],'-',totest,'-',toprocess,'-',LorR,'-BICuminusc','.pdf'),width=9,height=3.5)
    gg <- ggplot(D, aes(x=evt_time,y=BICu-BICc,group=network,color=network)) + facet_wrap(~HC_region)+
    geom_line(size=1) + scale_color_manual(values = pal,labels=c('DMN','CTR','LIM')) + 
      xlab(epoch_label) + ylab('BIC: Uncompressed - Compressed')
    print(gg)
    dev.off()
  }
}


setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback/')
toalign <- 'feedback'
toprocess <- 'network'
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='vmPFC',toalign=toalign,toprocess=toprocess)
setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
toalign <- 'clock'
toprocess <- 'network'
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='vmPFC',toalign=toalign,toprocess=toprocess)
setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC/')
toalign <- 'feedback'
toprocess <- 'region'
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'-L','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-L','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='HC',toalign=toalign,toprocess=toprocess)
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='HC',toalign=toalign,toprocess=toprocess)
setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC/')
toalign <- 'clock'
toprocess <- 'region'
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'-L','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-L','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='HC',toalign=toalign,toprocess=toprocess)
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='HC',toalign=toalign,toprocess=toprocess)
setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction/')
toalign <- 'feedback'
toprocess <- 'network-by-HC'
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'-L','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-L','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='vmPFC-HC',toalign=toalign,toprocess=toprocess)
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='vmPFC-HC',toalign=toalign,toprocess=toprocess)
setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction/')
toalign <- 'clock'
toprocess <- 'network-by-HC'
#U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-network-by-HC','-L','.Rdata')
#U_file <- Sys.glob(paste0('*',U_file))
#C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-network-by-HC','-L','.Rdata')
#C_file <- Sys.glob(paste0('*',C_file))
#plot_BIC_diff(U_file=U_file,C_file=C_file,struct='vmPFC-HC',toalign=toalign)
U_file <- paste0('-','uncompressed','-',toalign,'-','base-wPE','-network-by-HC','-R','.Rdata')
U_file <- Sys.glob(paste0('*',U_file))
C_file <- paste0('-','compressed','-',toalign,'-','base-wPE','-',toprocess,'-R','.Rdata')
C_file <- Sys.glob(paste0('*',C_file))
plot_BIC_diff(U_file=U_file,C_file=C_file,struct='vmPFC-HC',toalign=toalign,toprocess=toprocess)


