LLgev <- function(tseq,yinput,h,scale.factor){

  library(evd)
  # scale factor divides y by some factor to help the functions run
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
  location <- parm[1]+parm[2]*(istar-n*t)
  scale <- exp(parm[3])
  c.phi <- 0.8
  b.phi <- -c.phi^(-1)*log(1-(1/2)^c.phi)*(1-(1/2)^c.phi)*2^(c.phi-1)
  a.phi <- -b.phi*log(-log(1-(1/2)^c.phi))
  shape <- (1-exp(-exp((parm[4]-a.phi)/b.phi)))^(1/c.phi)-(1/2)
  dgevvec=dgev(y[istar],loc=location,scale=scale,shape=shape)
  dgevvec[dgevvec==0]=NA
  dflow=log(dgevvec)
  # Negative of local likelihood for optim function
  lflow=as.numeric(-dflow[!is.na(dflow)]%*%wstar[!is.na(dflow)])
  return(lflow)
}

yinput <- yinput/scale.factor
outparm=matrix(ncol=7,nrow=length(tseq))
colnames(outparm)=c("Beta_0","Beta_1","Gamma_0","Phi_0","Q99","Median","Mean")
for(i in 1:length(tseq)) {
  n <- length(yinput)
  t <- tseq[i]
  lowerl <- floor(n*(t-h))+1
  upperl <- floor(n*(t+h))
  istar <- max(1,lowerl):min(n,upperl)
  mle.fit <- as.numeric(fgev(yinput[istar],nsloc=rep((yinput[max(istar)]-yinput[min(istar)])/(max(istar)-min(istar)),length(yinput[istar])),std.err=FALSE)$estimate)
  if (mle.fit[4] > 0) {
    mle.fit[4] <- min(mle.fit[4],0.45)
  } else {
    mle.fit[4] <- max(mle.fit[4],-0.45)
  }
  c.phi <- 0.8
  b.phi <- -c.phi^(-1)*log(1-(1/2)^c.phi)*(1-(1/2)^c.phi)*2^(c.phi-1)
  a.phi <- -b.phi*log(-log(1-(1/2)^c.phi))
  ip4 <- a.phi+b.phi*log(-log(1-(mle.fit[4]+1/2)^c.phi))
  initial.parm <- c(mle.fit[1],mle.fit[2],log(mle.fit[3]),ip4)
  outparm[i,1:4]=optim(par=initial.parm,control=list(maxit=20000),local_like,t=tseq[i],h=h,y=yinput)$par
  location <- scale.factor*outparm[i,1]
  scale <- scale.factor*exp(outparm[i,3])
  shape <- (1-exp(-exp((outparm[i,4]-a.phi)/b.phi)))^(1/c.phi)-(1/2)
  outparm[i,5]=qgev(0.99,loc=location,scale=scale,shape=shape)
  outparm[i,6]=location+scale*(log(2)^(-shape)-1)/shape
  outparm[i,7]=location+scale*(gamma(1-shape)-1)/shape
}
return(cbind(t=tseq,outparm))
}