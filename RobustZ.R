# RobustZ
# 2022-06-29 AndyP

RobustZ <- function(x){

mu <- 0
lastmu <- Inf
keep <- rep(1,length.out=length(x))

while (abs(mu-lastmu)>0){
  
  m <- median(x[keep==1],na.rm=TRUE)
  q <- quantile(abs(x-m),na.rm=TRUE)
  lo <- q[[1]]
  hi <- q[[3]]
  
  keep <- abs(x-m)>lo & abs(x-m)<hi
  
  mu <- mean(x[keep==1],na.rm=TRUE)
  sigma <- sd(x[keep==1],na.rm=TRUE)
  
  lastmu <- mu
  
}

z <- (x-mu)/sigma
return(z)

}