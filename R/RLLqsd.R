RLLqsd <- function(t,yinput,wts,h,scale.factor,parm){
  
  source("RLL.sandwich.mat.R")
  
  beta0 <- as.numeric(parm[1])
  beta1 <- as.numeric(parm[2])
  gamma0 <- as.numeric(parm[3])
  c.phi <- 0.8
  b.phi <- -c.phi^(-1)*log(1-(1/2)^c.phi)*(1-(1/2)^c.phi)*2^(c.phi-1)
  a.phi <- -b.phi*log(-log(1-(1/2)^c.phi))
  delta0 <- (1-exp(-exp((as.numeric(parm[4])-a.phi)/b.phi)))^(1/c.phi)-(1/2)
  
  # Estimating the quantile at t.star so (t.star-t)=0
  
  grad.gamma0 <- exp(gamma0)/delta0*((-log(0.99))^(-delta0)-1)
  grad.delta0 <- -exp(gamma0)/delta0^2*((-log(0.99))^(-delta0)-1)-exp(gamma0)/delta0*(log(-log(0.99))*(-log(0.99))^(-delta0))
  
  grad <- c(1,0,grad.gamma0,grad.delta0)
  tgrad <- t(grad)
  
  sd <- sqrt(tgrad%*%RLL.sandwich.mat(t,yinput,wts,h,scale.factor,parm)%*%grad)
  
  as.numeric(sd*scale.factor)
}
  
