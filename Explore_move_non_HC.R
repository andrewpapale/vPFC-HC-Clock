# 2025-05-20 AndyP
# Explore move non-HC aside

demo <- readRDS('/Volumes/Users/Andrew/MEDuSA_data_Explore/explore_n146.rds')
HC <- demo %>% filter(registration_group == 'HC')


my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}


ids <- HC$registration_redcapid


for (iD in 1:length(ids)){
  print(paste0("/Users/dnplserv/gPPI/Explore",'/sub-',ids[iD]))
  
  dir0 <- list.files(paste0("/Users/dnplserv/gPPI/Explore",'/sub-',ids[iD]),recursive=TRUE,full.names = TRUE)
  
  for (iD0 in 1:length(dir0)){
    temp <- str_split(dir0[iD0],'/')
    temp0 <- paste0(temp[[1]][6],'/',temp[[1]][7],'/',temp[[1]][8])
    print(temp0)
    my.file.rename(from = dir0[iD0],
                   to = paste0("/Users/dnplserv/gPPI/Explore_HC_only/",temp0))

  }
}
