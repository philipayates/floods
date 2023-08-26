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
  #location <- parm[1]+parm[2]*(istar/n-t)
  scale <- exp(parm[3])
  #scale <- parm[3]
  shape <- parm[4]
  dgevvec=dgev(y[istar],loc=location,scale=scale,shape=shape)
  dgevvec[dgevvec==0]=NA
  dflow=log(dgevvec)
  # Negative of local likelihood for optim function
  lflow=as.numeric(-dflow[!is.na(dflow)]%*%wstar[!is.na(dflow)])
  return(lflow)
}

yinput <- yinput/scale.factor
outparm=matrix(ncol=7,nrow=length(tseq))
colnames(outparm)=c("Beta_0","Beta_1","Gamma_0","Delta_0","Q99","Median","Mean")
for(i in 1:length(tseq)) {
  n <- length(yinput)
  t <- tseq[i]
  lowerl <- floor(n*(t-h))+1
  upperl <- floor(n*(t+h))
  istar <- max(1,lowerl):min(n,upperl)
  mle.fit <- as.numeric(fgev(yinput[istar],nsloc=rep((yinput[max(istar)]-yinput[min(istar)])/(max(istar)-min(istar)),length(yinput[istar])),std.err=FALSE)$estimate)
  #initial.parm <- mle.fit
  initial.parm <- c(mle.fit[1],mle.fit[2],log(mle.fit[3]),mle.fit[4])
  outparm[i,1:4]=optim(par=initial.parm,control=list(maxit=20000),local_like,t=tseq[i],h=h,y=yinput)$par
  location <- outparm[i,1]
  scale <- exp(outparm[i,3])
  #scale <- outparm[i,3]
  shape <- outparm[i,4]
  outparm[i,5]=scale.factor*qgev(0.99,loc=location,scale=scale,shape=shape)
  outparm[i,6]=scale.factor*(location+scale*(log(2)^(-shape)-1)/shape)
  outparm[i,7]=scale.factor*(location+scale*(gamma(1-shape)-1)/shape)
}
return(cbind(t=tseq,outparm))
}