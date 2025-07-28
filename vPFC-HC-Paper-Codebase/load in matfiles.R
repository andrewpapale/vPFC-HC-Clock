library(tidyverse)
library(R.matlab)
setwd("~/Documents/Projects/Zebra/recompute bsocial model r2/")
files <- list.files("simplified subj mat files/", full.names = TRUE)
test_mat <- readMat(files[[1]])
comp_robcor <- function(y, yhat) {
  rob_r = robust::covRob(data.frame(y = as.vector(y), yhat = as.vector(yhat)), corr = TRUE, estim = "donostah")
  rob_r_val = rob_r$cov[[2,1]]
}

load_mat <- function(fname) {
  m <- readMat(fname) 
  model_r2 <- m$r2[[1]]
  y <- t(as.matrix(m$y)) %>% as.data.frame() 
  acc <- m$acc[[1]]
  bacc <- m$bacc[[1]]
  names(y) <- paste0("timebin_", 1:ncol(y))
  y$asc_trial <- 1:nrow(y)
  subj_id <- gsub("_trialwise.mat", "",fname)
  subj_id <- gsub("simplified subj mat files//", "",subj_id, fixed = TRUE)
  y$id <- subj_id
  y_long <- gather(y, key = "key", value = "value", -id, -asc_trial) %>% mutate(timebin = gsub("timebin_", "", key) , timebin = as.numeric(timebin)) 
  y_trialwise <- group_by(y_long, id, asc_trial) %>% mutate(max_y = max(value)) %>% dplyr::filter(max_y == value)
  y_trialwise <- dplyr::select(y_trialwise, id, asc_trial, timebin) %>% arrange(id, asc_trial)
  y_trialwise <- rename(y_trialwise, y = timebin)
  
  yhat <- t(as.matrix(m$yhat)) %>% as.data.frame()
  names(yhat) <- paste0("timebin_", 1:ncol(yhat))
  yhat$asc_trial <- 1:nrow(y)
  yhat$id <- subj_id
  yhat_long <- gather(yhat, key = "key", value = "value", -id, -asc_trial) %>% mutate(timebin = gsub("timebin_", "", key) , timebin = as.numeric(timebin)) 
  yhat_trialwise <- group_by(yhat_long, id, asc_trial) %>% mutate(max_y = max(value)) %>% dplyr::filter(max_y == value)
  dups_info <- yhat_trialwise %>% group_by(id, asc_trial) %>% summarise(count = n()) %>% ungroup()
  dup_trials <- dplyr::filter(dups_info, count > 1)
  if(nrow(dup_trials) > 0) {
    yhat_trialwise <- dplyr::filter(yhat_trialwise, !asc_trial %in% dup_trials$asc_trial)
  }
  yhat_trialwise <- dplyr::select(yhat_trialwise, id, asc_trial, timebin) %>% arrange(id, asc_trial)
  yhat_trialwise <- rename(yhat_trialwise, yhat = timebin)
  
  info <- full_join(y_trialwise, yhat_trialwise) 
  
  cor <- cor(info$y, info$yhat, use = "complete.obs")
  for_robcor <- info
  for_robcor <- na.omit(for_robcor)
  rob_r <- comp_robcor(for_robcor$y, for_robcor$yhat)
  info$r <- cor
  r2 <- cor^2
  rob_r2 <- rob_r^2
  info$rob_r <- rob_r
  info$r2 <- r2
  info$rob_r2 <- rob_r2
  info$vba_r2 <- model_r2
  info$acc <- acc
  info$bacc <- bacc
  return(info)
}
allsubjs<- plyr::ldply(.data = files, .fun = load_mat)
allsubjs_statinfo <- dplyr::select(allsubjs, id, r, rob_r, r2, rob_r2, vba_r2, acc, bacc) %>% distinct()

pdata <- read_csv("../data/bsocial-clock/raw and prepped/pdata_df.csv")
trust_pdata <- read_csv("../raw_data/agree_factors_6-21-22.csv")

stats_w_pdata <- left_join(allsubjs_statinfo %>% mutate(id = gsub("_1", "", id) %>% as.numeric), dplyr::select(pdata,id, pid5_callousness, pid5_manipulativeness))
stats_w_pdata2 <- left_join(allsubjs_statinfo %>% mutate(id = gsub("_1", "", id) %>% as.numeric), dplyr::select(trust_pdata,id, exploit2f, callous2f))


a <- ggplot(allsubjs_statinfo, aes(x = r)) + geom_histogram() + ggtitle("Correlation between observed and predicted time bin") + theme_bw()
b <- ggplot(allsubjs_statinfo, aes(x = r2)) + geom_histogram() + ggtitle("Amount of variance explained") + theme_bw()
c <- ggplot(allsubjs_statinfo, aes(x = rob_r)) + geom_histogram() + ggtitle("Robust correlation between observed and\npredicted time bin") + theme_bw()
d <- ggplot(allsubjs_statinfo, aes(x = rob_r2)) + geom_histogram() + ggtitle("Amount of variance explained - robust r") + theme_bw()
abcd <- cowplot::plot_grid(a,b,c,d, ncol = 2)
e <- ggplot(allsubjs_statinfo, aes(x = r, vba_r2)) + geom_point() + geom_smooth(method = "lm") + ggtitle("Association between r and vba r2")+ theme_bw()
f <- ggplot(allsubjs_statinfo, aes(x = acc)) + geom_density()+ ggtitle("Classification accuracy")+ theme_bw()
# honestly looks very good
pdf("bsocial re computed r vs. r2.pdf", width = 10, height = 10)
plot(abcd)
plot(e)
dev.off()
cor(allsubjs_statinfo$r, allsubjs_statinfo$vba_r2) # r = .77

# bsocial_trust_r2 <- read_csv("~/Downloads/trust_cleaned.csv") %>% dplyr::select(id, R2) %>% distinct()
# ggplot(bsocial_trust_r2, aes(x = R2)) + geom_histogram()
# 
# allsubjs_statinfo <- mutate(allsubjs_statinfo, id = gsub("_1", "", id), id = as.numeric(id))
# 
# bsocial_trust_r2 <- rename(bsocial_trust_r2, trust_vba_r2 = R2)
# allsubjs_statinfo <- left_join(allsubjs_statinfo, bsocial_trust_r2)
# 
# ggplot(allsubjs_statinfo, aes(x = vba_r2, trust_vba_r2)) + geom_point() + geom_smooth(method = "lm")
# cor(allsubjs_statinfo$vba_r2, allsubjs_statinfo$trust_vba_r2, use = "complete.obs") # r = .24
# cor(allsubjs_statinfo$r2, allsubjs_statinfo$trust_vba_r2, use = "complete.obs") # r = .21
# ggplot(allsubjs_statinfo, aes(x = r2, trust_vba_r2)) + geom_point() + geom_smooth(method = "lm")
# 
# 
# df <- read_csv("~/Downloads/trust_cleaned.csv") %>% dplyr::select(id, R2) %>% rename(trust_vba_r2 = R2) %>% mutate(task = "trust")
# 
# orderinfo <- readLines("~/Downloads/all_bsocial_mrdirs.txt")
# orderinfo_df <- data.frame(dir = as.vector(orderinfo))
# orderinfo_df <- mutate(orderinfo_df, id = str_sub(dir, start = 3, 8))
# 
# orderinfo_df <- separate(orderinfo_df, dir, into = c("root", "top", "mid", "lower"), sep = "\\/")
# orderinfo_df <- dplyr::select(orderinfo_df, top, mid, lower, id) %>% mutate(lower2 = lower) %>% separate(lower2, into = c("runnum"), sep = "-")
# orderinfo_df <- mutate(orderinfo_df, runnum = as.numeric(runnum), task = case_when(str_detect(lower, "trust") ~ "trust",
#                                                       str_detect(lower, "clock") ~ "clock", 
#                                                       str_detect(lower, "spott") ~ "spott"))
# 
# 
# orderinfo_df_sub <- dplyr::filter(orderinfo_df, !is.na(task))
# orderinfo_df_sub <- dplyr::filter(orderinfo_df_sub, !grepl("gre", lower),
#                                   !grepl("SBRef", lower))
# order_df_summarized <- group_by(orderinfo_df_sub, id, top, mid, task) %>% summarise(lowest_run_win_task = min(runnum)) %>% ungroup()
# order_df_summarized <- group_by(order_df_summarized, id, top, mid) %>% mutate(n_tasks_during_scan = n()) %>% ungroup()
# order_df_summarized <- group_by(order_df_summarized, id, top, mid) %>% mutate(order = rank(lowest_run_win_task)) %>% ungroup()
# clock_order_info <- mutate(order_df_summarized %>% dplyr::filter(task == "clock"), is_clock_first = if_else(order == 1, 1, 0),
#                               is_clock_second = if_else(order == 2, 1, 0),
#                               is_clock_third = if_else(order == 3, 1, 0))
# 
# trust_order_info <- mutate(order_df_summarized %>% dplyr::filter(task == "trust"), is_trust_first = if_else(order == 1, 1, 0),
#                            is_trust_second = if_else(order == 2, 1, 0),
#                            is_trust_third = if_else(order == 3, 1, 0))
# 
# allsubjs_statinfo <- allsubjs_statinfo %>% mutate(id = gsub("_1", "", id), id = as.numeric(id))
# clock_order_info <- clock_order_info  %>% mutate(id = as.numeric(id))
# allsubjs_statinfo <- left_join(allsubjs_statinfo, clock_order_info)
# # not perfect but can at least see what happens
# 
# ggplot(allsubjs_statinfo,aes(x= r2, color = as.factor(is_clock_first), fill = as.factor(is_clock_first))) + geom_density(alpha = .1)
# 
# trust_order_info <- trust_order_info  %>% mutate(id = as.numeric(id))
# df <- left_join(df, trust_order_info)
# 
# ggplot(df,aes(x= trust_vba_r2, color = as.factor(is_trust_first), fill = as.factor(is_trust_first))) + geom_density(alpha = .1)
# 
# order_df_summarized <- dplyr::select(order_df_summarized, id, mid, task, n_tasks_during_scan, order) %>% rename(subfolder = mid)
# write_csv(order_df_summarized, "bsocial order of tasks - pulled from bierka.csv")
