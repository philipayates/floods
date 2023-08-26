SimC=read.csv("/Users/philipayates/Dropbox/Phil Yates/sim-modelC.csv")
#Study first X column
library(dplyr)
library(ggplot2)
library(FAdist)
#ggplot(SimC, aes(x=t, y=X)) + geom_point() + ylab("Simulated Observations from GEV")

local_like=function(parm,t,h,y){
  n=length(y)
  #Bounds of smoothing window
  lowerl=floor(n*(t-h))+1
  upperl=floor(n*(t+h))
  #Account for out-of-range subscripts
  istar=max(1,lowerl):min(n,upperl)
  #Weights for smoothing window
  wstar=(1-((istar/n-t)/h)^2)
  #Likelihood for smoothing window
  location <- parm[1]
  scale <- exp(parm[2])
  shape <- parm[3]
  z <- (y[istar]-location)/scale
  x <- 1+shape*z
  # Negative of local likelihood for optim function
  lflow <- 0
  for(i in 1:length(x)){
    if(x[i]>0){
      lflow <- lflow+(wstar[i]*(log(scale)+(1+1/shape)*log(x[i])+x[i]^(-1/shape)))
    }
  }
  return(lflow)
}

local_like_R=function(parm,t,h,y,rweight=1){
  n=length(y)
  #Bounds of smoothing window
  lowerl=floor(n*(t-h))+1
  upperl=floor(n*(t+h))
  #Account for out-of-range subscripts
  istar=max(1,lowerl):min(n,upperl)
  #Weights for smoothing window
  wstar=rweight*(1-((istar/n-t)/h)^2)
  
  #Likelihood for smoothing window
  #eta is location (betas), tau is scale (gammas), kappa is shape (deltas)
  location <- parm[1]
  scale <- exp(parm[2])
  shape <- parm[3]
  z <- (y[istar]-location)/scale
  x <- 1+shape*z
  # Negative of local likelihood for optim function
  lflow <- 0
  for(i in 1:length(x)){
    if(x[i]>0){
      lflow <- lflow+(wstar[i]*(log(scale)+(1+1/shape)*log(x[i])+x[i]^(-1/shape)))
    }
  }
  return(lflow)
}

t <- SimC[,1]
eta.B <- 0.25*t+0.1*dnorm(t,0.5,0.1)+0.1*dnorm(t,1,0.1)
tau.B <- exp(-1.5+0.2*t)
kappa.B <- rep(-0.1,length(t))

trunsLL=function(tseq,yinput){
  outparm=matrix(ncol=6,nrow=length(tseq))
  conv <- NULL
  colnames(outparm)=c("eta","tau","kappa","Q99","Median","Mean")
  for(i in 1:length(tseq)) {
    conv[i]=optim(par=c(eta.B[i],log(tau.B[i]),kappa.B[i]),control=list(maxit=2000),local_like,t=tseq[i],h=0.3,y=yinput)$convergence
    outparm[i,1:3]=optim(par=c(eta.B[i],log(tau.B[i]),kappa.B[i]),control=list(maxit=2000),local_like,t=tseq[i],h=0.3,y=yinput)$par
    loc <- outparm[i,1]
    scale <- exp(outparm[i,2])
    shape <- outparm[i,3]
    outparm[i,4]=qgev(0.99,location=loc,scale=scale,shape=shape)
    outparm[i,5]=loc+scale*(log(2)^(-shape)-1)/shape
    outparm[i,6]=loc+scale*(gamma(1-shape)-1)/shape
  }
  return(cbind(tseq,conv,outparm))
}

trunsRLL=function(tseq,yinput){
  outparm=matrix(ncol=6,nrow=length(tseq))
  #conv <- NULL
  outwt=matrix(ncol=length(yinput),nrow=length(tseq))
  n=length(yinput)
  colnames(outparm)=c("eta","tau","kappa","Q99","Median","Mean")
  for(i in 1:length(tseq)) {
    diff=1
    wstar=1
    t0=tseq[i]
    #print(t0)
    h0=0.3
    lowerl=floor(n*(t0-h0))+1
    upperl=floor(n*(t0+h0))
    istar=max(1,lowerl):min(n,upperl)
    parm_init=c(eta.B[i],log(tau.B[i]),kappa.B[i])
    ictr=1
    while(max(abs(diff))>1.e-4){
      parmhat=optim(par=parm_init,control=list(maxit=10000),local_like_R,t=t0,h=h0,y=yinput,rweight=wstar)
      parm_new=parmhat$par
      muhat=parm_new[1]+ parm_new[2]*(gamma(1-parm_new[3])-1)/parm_new[3]
      resid=yinput[istar]-muhat
      sixs=8*median(abs(resid))
      wstar=(1-(pmin(1,abs(resid/sixs))^2)^2)  
      diff=parm_init-parm_new
      parm_init=parm_new
      ictr=ictr+1
    }
    
    outparm[i,1:3]=parm_new
    outparm[i,4]=qgev(0.99,location=outparm[i,1],scale=exp(outparm[i,2]),shape=outparm[i,3])
    loc <- outparm[i,1]
    scale <- exp(outparm[i,2])
    shape <- outparm[i,3]
    outparm[i,5]=loc+scale*(log(2)^(-shape)-1)/shape
    outparm[i,6]=loc+scale*(gamma(1-shape)-1)/shape
    outwt[i,istar]=wstar
    outlist=list(est=cbind(tseq,outparm),wts=outwt)
  }
  return(outlist)
}

outCX <- NULL
for(i in 1:1000){
  sim <- i
  truns <- trunsLL(t,SimC[,i+1])
  out.truns <- cbind(rep(sim,100),truns)
  outCX <- rbind(outCX,out.truns)
}

colnames(outCX)[1] <- "sim"
write.table(outCX,"/Users/philipayates/Dropbox/Phil Yates/Sim Results/outCX-h30.csv",sep=",",row.names=FALSE)

outCXR <- NULL
for(i in 1:1000){
  sim <- i
  truns.RLL <- trunsRLL(t,SimC[,i+1])
  out.truns.RLL <- cbind(rep(sim,100),truns.RLL$est)
  outCXR <- rbind(outCXR,out.truns.RLL)
}

colnames(outCXR)[1] <- "sim"
write.table(outCXR,"/Users/philipayates/Dropbox/Phil Yates/Sim Results/outCXR-h30.csv",sep=",",row.names=FALSE)

#write.table(outAXR[[1]],"outAXR.csv",sep=",",row.names=FALSE)
#write.table(outAXR[[2]],"outAXRwt.csv",sep=",",row.names=FALSE)

