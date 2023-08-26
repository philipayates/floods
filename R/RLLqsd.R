RLLqsd <- function(t,yinput,wts,h,scale.factor,parm){
  
  source("~//Dropbox/Phil Yates/Phil New Functions/RLL.sandwich.mat.R")
  
  beta0 <- as.numeric(parm[1])
  beta1 <- as.numeric(parm[2])
  gamma0 <- as.numeric(parm[3])
  delta0 <- as.numeric(parm[4])
  
  # Estimating the quantile at t.star so (t.star-t)=0
  
  grad.gamma0 <- exp(gamma0)/delta0*((-log(0.99))^(-delta0)-1)
  grad.delta0 <- -exp(gamma0)/delta0^2*((-log(0.99))^(-delta0)-1)-exp(gamma0)/delta0*(log(-log(0.99))*(-log(0.99))^(-delta0))
  
  grad <- c(1,0,grad.gamma0,grad.delta0)
  tgrad <- t(grad)
  
  sd <- sqrt(tgrad%*%RLL.sandwich.mat(t,yinput,wts,h,scale.factor,parm)%*%grad)
  
  as.numeric(sd*scale.factor)
}
  