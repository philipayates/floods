l.cv <- function(t,y,h,scale.factor){
  
  source("~//Dropbox/Phil Yates/July 2021 Current Code/Update/RLLgev.R")
  
  n <- length(y)
  CV.test <- matrix(NA,nrow=n,ncol=5)
  colnames(CV.test)=c("t","Beta_0","Beta_1","Gamma_0","Phi_0")
  
  for(i in 1:n){
    CV.test[i,] <- RLLgev(t[i],y[-i],h,scale.factor)$parm.est[1:5]
  }
  
  c.phi <- 0.8
  b.phi <- -c.phi^(-1)*log(1-(1/2)^c.phi)*(1-(1/2)^c.phi)*2^(c.phi-1)
  a.phi <- -b.phi*log(-log(1-(1/2)^c.phi))
  
  location <- scale.factor*CV.test[,2]
  scale <- scale.factor*exp(CV.test[,4])
  shape <- (1-exp(-exp((CV.test[,5]-a.phi)/b.phi)))^(1/c.phi)-(1/2)
  
  
  like.cv <- NULL
  for(i in 1:n){
    like.cv[i] <- dgev(y[i],loc=location[i],scale=scale[i],shape=shape[i])
  }

  l.cv <- sum(log(like.cv[like.cv>0]))
  
return(l.cv)
}