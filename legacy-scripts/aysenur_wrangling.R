# 2023-08-15 AndyP edited AO script

##Hippocampal decons:  #does not have "trial", uses "run_trial" so we add "trial". It also does not have decon_median values - can be calculated if need be.
load("/Volumes/Users/Andrew/MMClock_10_2023/HC_fb_Aug2023.Rdata")
hc <- hc  %>% dplyr::select(-atlas_value0, -HC_region, -bin_num) #evt_time [-6:9]
#Get bins first since atlas_values are different contralaterally so we cannot get combined decons based on them
hc <- hc %>% mutate(bin = paste0("Bin", cut(atlas_value, breaks = 12, labels = 1:12)))
hippo1 <- hc %>% group_by(id, run, trial, evt_time, bin, side) %>% summarise(decon_mean1 = mean(decon_mean, na.rm=TRUE)) %>% ungroup() 
hippo2 <- hippo1 %>% group_by(id,run,trial,evt_time,bin) %>% summarize(decon_mean2 = mean(decon_mean1,na.rm=TRUE))%>% ungroup()
hippo2 <- hippo2 %>% mutate(HC_region = case_when(bin=='Bin1'~'PH',
                                                  bin=='Bin2'~'PH',
                                                  bin=='Bin3'~'PH',
                                                  bin=='Bin4'~'PH',
                                                  bin=='Bin5'~'AH',
                                                  bin=='Bin6'~'AH',
                                                  bin=='Bin7'~'AH',
                                                  bin=='Bin8'~'AH',
                                                  bin=='Bin9'~'AH',
                                                  bin=='Bin10'~'AH',
                                                  bin=='Bin11'~'AH',
                                                  bin=='Bin12'~'AH'))
hippo3 <- hippo2 %>% group_by(id,run,trial,evt_time,HC_region) %>% summarize(decon_mean3 = mean(decon_mean2,na.rm=TRUE)) %>% ungroup()
hippo4 <- hippo3 %>% group_by(evt_time,HC_region) %>% summarize(decon_mean4 = mean(decon_mean3,na.rm=TRUE),decon_sd = sd(decon_mean3,na.rm=TRUE),N=length(decon_mean3)) %>% ungroup()
ggplot(hippo4,aes(x=evt_time,y=decon_mean4)) + geom_line() + geom_point() + facet_wrap(~HC_region) + geom_errorbar(aes(ymin=decon_mean4-(decon_sd/sqrt(N)),ymax=decon_mean4+(decon_sd/sqrt(N))))

#hippo1 <- hippo1 %>% pivot_wider(names_from = side, values_from = decon_mean)
#hippo2 <- hippo1 %>% rowwise() %>% mutate(decon_mean = mean(c(l, r), trim=0, na.rm = TRUE)) %>% dplyr::select(-r, -l)
#hippo2 <- hippo2 %>% pivot_wider(names_from = "bin", values_from = "decon_mean")
#I had used last 3 bins to define subregions. However, I am starting to use 1/3 and 2/3 based on a recent discussion with Michael on 7/2023
#hippo3 <- hippo2 %>% rowwise %>% mutate(PH = mean(c(Bin1, Bin2, Bin3, Bin4), trim=0, na.rm = TRUE), AH = mean(c(Bin5, Bin6, Bin7, Bin8, Bin9, Bin10, Bin11, Bin12), trim=0, na.rm = TRUE)) %>% ungroup()
#hippo4 <- hippo3 %>% pivot_longer(names_to = "atlas_value", values_to = "decon_mean", (5:18))

#if this is a new session, load behavioral data, get lags, merge with signals (code above to load the df). If not, you already have the df