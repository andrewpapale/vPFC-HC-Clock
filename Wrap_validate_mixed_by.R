# 2022-03-17 AndyP
# Wrapper function for validate_mixed_by vmPFC, HC, and vmPFC/HC functions
# 

doHC <- FALSE
dovmPFC <- FALSE
doHC_vmPFC <- TRUE

doclock <- TRUE
dofeedback <- TRUE
dortvmax <- FALSE

load_existing <- FALSE
reanalyze = TRUE

ncores = 26
# same repo throughout
repo_directory <- "~/clock_analysis"
HC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'
vmPFC_cache_dir = '~/vmPFC/MEDUSA Schaefer Analysis'

# these are the same for all 3 functions
totest <-  c('base','base-wPE')
behavmodel <- c('compressed','uncompressed')
hc_LorR <- c('LR') # vmPFC does not have this as an input

source('~/vmPFC/validate_mixed_by_HC.R')
source('~/vmPFC/validate_mixed_by_vmPFC.R')
source('~/vmPFC/validate_mixed_by_vmPFC_HC.R')


if (doHC==TRUE && dofeedback==TRUE){
  toprocess <- c('axis','region') # for HC
  toalign <- 'feedback'
  iC <-1
  for (i in toprocess){
    for (j in totest){
      for (k in behavmodel){
        for (z in hc_LorR){
          if (iC==1){
            reload=TRUE
            Q <- NULL
          } else {
            reload=TRUE
          }
          tryCatch({
            if (reanalyze==FALSE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)>1){
                message('multiple ddf_files found')
                stopifnot(length(ddf_file)==1)
              } else if (length(ddf_file)==1){
                skip=TRUE
              } else if (length(ddf_file)==0){
                skip=FALSE
              }
            } else {
              skip = FALSE
            }
            if (load_existing==TRUE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'-',z,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)!=1){
                message('multiple ddf_files found')
              }
            } else {
              ddf_file <- ''
            }
            if (!skip){
              Q <- validate_mixed_by_HC(ncores=ncores,toalign=toalign,toprocess=i,totest=j,behavmodel=k,repo_directory=repo_directory,
                                        HC_cache_dir=HC_cache_dir,hc_LorR=z,reload=reload,replot=TRUE,load_existing=load_existing,ddf_file=ddf_file,Q=NULL)
            }
          },
          error=function(cond){
            message(paste0("error on:",toprocess,' ', totest, ' ', behavmodel, ' ', hc_LorR))
            message(cond)
            return(NA)
          })
          iC<-iC+1;
        }
      }
    }
  }
}
if (doHC==TRUE && doclock==TRUE){
  toprocess <- c('axis','region') # for HC
  toalign <- 'clock'
  iC <- 1
  for (i in toprocess){
    for (j in totest){
      for (k in behavmodel){
        for (z in hc_LorR){
          if (iC==1){
            reload=TRUE
            Q<-NULL
          } else {
            reload=TRUE
          }
          tryCatch({
            if (reanalyze==FALSE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)>1){
                message('multiple ddf_files found')
                stopifnot(length(ddf_file)==1)
              } else if (length(ddf_file)==1){
                skip=TRUE
              } else if (length(ddf_file)==0){
                skip=FALSE
              }
            } else {
              skip = FALSE
            }
            if (load_existing==TRUE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'-',z,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)!=1){
                message('multiple ddf_files found')
              }
            } else {
              ddf_file <- ''
            }
            if (!skip){
              Q <- validate_mixed_by_HC(ncores=ncores,toalign=toalign,toprocess=i,totest=j,behavmodel=k,repo_directory=repo_directory,
                                        HC_cache_dir=HC_cache_dir,hc_LorR=z,reload=reload,replot=TRUE,load_existing=load_existing,ddf_file=ddf_file,Q=NULL)
            }
          },
          error=function(cond){
            message(paste("error on:",toprocess,' ', totest, ' ', behavmodel, ' ', hc_LorR))
            message(cond)
            return(NA)
          })
          iC <- iC + 1
        }
      }
    }
  }
}
if (dovmPFC==TRUE && dofeedback==TRUE){
  toprocess <- c('network','symmetry_group') # for vmPFC
  toalign <- 'feedback'
  iC <- 1
  for (i in toprocess){
    for (j in totest){
      for (k in behavmodel){
        if (iC==1){
          reload=TRUE
          Q <- NULL
        } else {
          reload=TRUE
        }
        tryCatch({
          if (reanalyze==FALSE){
            setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback/')
            ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
            ddf_file <- Sys.glob(paste0('*',ddf_file))
            print(ddf_file)
            if (length(ddf_file)>1){
              message('multiple ddf_files found')
              stopifnot(length(ddf_file)==1)
            } else if (length(ddf_file)==1){
              skip=TRUE
            } else if (length(ddf_file)==0){
              skip=FALSE
            }
          } else {
            skip = FALSE
          }
          if (load_existing==TRUE){
            setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback/')
            ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'-','.Rdata')
            ddf_file <- Sys.glob(paste0('*',ddf_file))
            print(ddf_file)
            if (length(ddf_file)!=1){
              message('multiple ddf_files found')
            }
          } else {
            ddf_file <- ''
          }
          if (!skip){
            Q <- validate_mixed_by_vmPFC(ncores=ncores,toalign=toalign,toprocess=i,totest=j,behavmodel=k,repo_directory=repo_directory,
                                         cache_dir=vmPFC_cache_dir,reload=FALSE,replot=TRUE,load_existing=load_existing,ddf_file=ddf_file,d=NULL)
          }
        },
        error=function(cond){
          message(paste("error on:",toprocess,' ', totest, ' ', behavmodel, ' ', hc_LorR))
          message(cond)
          return(NA)
        })
        iC <- iC + 1
      }
    }
  }
}

if (dovmPFC==TRUE && doclock==TRUE){
  toprocess <- c('network','symmetry_group') # for vmPFC
  toalign <- 'clock'
  iC <- 1
  for (i in toprocess){
    for (j in totest){
      for (k in behavmodel){
        if (iC==1){
          reload=TRUE
          Q <- NULL
        } else {
          reload=TRUE
        }
        tryCatch({
          if (reanalyze==FALSE){
            setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
            ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
            ddf_file <- Sys.glob(paste0('*',ddf_file))
            print(ddf_file)
            if (length(ddf_file)>1){
              message('multiple ddf_files found')
              stopifnot(length(ddf_file)==1)
            } else if (length(ddf_file)==1){
              skip=TRUE
            } else if (length(ddf_file)==0){
              skip=FALSE
            }
          } else {
            skip = FALSE
          }
          if (load_existing==TRUE){
            setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock/')
            ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
            ddf_file <- Sys.glob(paste0('*',ddf_file))
            print(ddf_file)
            if (length(ddf_file)!=1){
              message('multiple ddf_files found')
              stopifnot(length(ddf_file)==1)
            }
          } else {
            ddf_file <- ''
          }
          if (!skip){
            Q <- validate_mixed_by_vmPFC(ncores=ncores,toalign=toalign,toprocess=i,totest=j,behavmodel=k,repo_directory=repo_directory,
                                         cache_dir=vmPFC_cache_dir,reload=reload,replot=TRUE,load_existing=load_existing,ddf_file=ddf_file, d=NULL)
          }
        },
        error=function(cond){
          message(paste("error on:",toprocess,' ', totest, ' ', behavmodel, ' ', hc_LorR))
          message(cond)
          return(NA)
        })
        iC <- iC + 1
      }
    }
  }
}
if (doHC_vmPFC==TRUE && dofeedback==TRUE){
  toprocess <- c('network-by-HC','symmetry-by-HC') # for vmPFC/HC
  toalign <- 'feedback'
  iC <- 1
  for (i in toprocess){
    for (j in totest){
      for (k in behavmodel){
        for (z in hc_LorR){
          if (iC==1){
            reload=TRUE
            Q <- NULL
          } else {
            reload=TRUE
          }
          tryCatch({
            if (reanalyze==FALSE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)>1){
                message('multiple ddf_files found')
                stopifnot(length(ddf_file)==1)
              } else if (length(ddf_file)==1){
                skip=TRUE
              } else if (length(ddf_file)==0){
                skip=FALSE
              }
            } else {
              skip = FALSE
            }
            if (load_existing==TRUE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_feedback_HC_interaction/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'-',z,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)!=1){
                message('multiple ddf_files found')
              }
            } else {
              ddf_file <- ''
            }
            if (!skip){
              Q <- validate_mixed_by_vmPFC_HC(ncores=ncores,toalign=toalign,toprocess=i,totest=j,behavmodel=k,repo_directory=repo_directory,
                                              vmPFC_cache_dir=vmPFC_cache_dir,HC_cache_dir=HC_cache_dir,hc_LorR=z,
                                              reload=reload,replot=TRUE,load_existing=load_existing,ddf_file=ddf_file,Q=NULL)
            }
          },
          error=function(cond){
            message(paste("error on:",toprocess,' ', totest, ' ', behavmodel, ' ', hc_LorR))
            message(cond)
            return(NA)
          })
          iC <- iC + 1
        }
      }
    }
  }
}
if (doHC_vmPFC==TRUE && doclock==TRUE){
  toprocess <- c('network-by-HC','symmetry-by-HC') # for vmPFC/HC
  toalign <- 'clock'
  iC <- 1
  for (i in toprocess){
    for (j in totest){
      for (k in behavmodel){
        for (z in hc_LorR){
          if (iC==1){
            reload=TRUE
            Q <- NULL
          } else {
            reload=TRUE
          }
          tryCatch({
            if (reanalyze==FALSE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)>1){
                message('multiple ddf_files found')
                stopifnot(length(ddf_file)==1)
              } else if (length(ddf_file)==1){
                skip=TRUE
              } else if (length(ddf_file)==0){
                skip=FALSE
              }
            } else {
              skip = FALSE
            }
            if (load_existing==TRUE){
              setwd('~/vmPFC/MEDUSA Schaefer Analysis/validate_mixed_by_clock_HC_interaction/')
              ddf_file <- paste0('-',k,'-',toalign,'-',j,'-',i,'-',z,'.Rdata')
              ddf_file <- Sys.glob(paste0('*',ddf_file))
              print(ddf_file)
              if (length(ddf_file)!=1){
                message('multiple ddf_files found')
              }
            } else {
              ddf_file <- ''
            }
            if (!skip){
              Q <- validate_mixed_by_vmPFC_HC(ncores=ncores,toalign=toalign,toprocess=i,totest=j,behavmodel=k,repo_directory=repo_directory,
                                              vmPFC_cache_dir=vmPFC_cache_dir,HC_cache_dir=HC_cache_dir,
                                              hc_LorR=z,reload=reload,replot=TRUE,load_existing=load_existing,ddf_file=ddf_file,Q=NULL)
            }
          },
          error=function(cond){
            message(paste("error on:",toprocess,' ', totest, ' ', behavmodel, ' ', hc_LorR))
            message(cond)
            return(NA)
          })
          iC <- iC + 1
        }
      }
    }
  }
}