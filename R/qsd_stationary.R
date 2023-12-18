qsd_stationary <- function(yinput,scale.factor,parm){
  
  #source("~//Dropbox/Phil Yates/Phil New Functions/sandwich.mat.R")
  #source("July 2021 Current Code/sandwich.mat.R")
  #source("~//Dropbox/Phil Yates/Phil New Functions/inform.mat.R")
  
  mu0 <- as.numeric(parm[1])
  sigma0 <- as.numeric(parm[2])
  delta0 <- as.numeric(parm[3])
  
  # Estimating the quantile at t.star so (t.star-t)=0
  
  grad.sigma0 <- 1/delta0*((-log(0.99))^(-delta0)-1)
  grad.delta0 <- -sigma0/delta0^2*((-log(0.99))^(-delta0)-1)-sigma0/delta0*(log(-log(0.99))*(-log(0.99))^(-delta0))
  
  grad <- c(1,grad.sigma0,grad.delta0)
  tgrad <- t(grad)
  
  sd <- sqrt(tgrad%*%sandwich_stationary.mat(yinput,scale.factor,parm)%*%grad)
  #sd <- sqrt(tgrad%*%inform.mat(t,yinput,h,scale.factor,parm)%*%grad)
  
  as.numeric(sd*scale.factor)
}
  