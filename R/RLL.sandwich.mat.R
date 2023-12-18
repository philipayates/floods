RLL.sandwich.mat <- function(t,yinput,wts,h,scale.factor,parm){
  y <- yinput/scale.factor
  beta0 <- as.numeric(parm[1])
  beta1 <- as.numeric(parm[2])
  gamma0 <- as.numeric(parm[3])
  c.phi <- 0.8
  b.phi <- -c.phi^(-1)*log(1-(1/2)^c.phi)*(1-(1/2)^c.phi)*2^(c.phi-1)
  a.phi <- -b.phi*log(-log(1-(1/2)^c.phi))
  delta0 <- (1-exp(-exp((as.numeric(parm[4])-a.phi)/b.phi)))^(1/c.phi)-(1/2)
  
  n <- length(y)
  lowerl=floor(n*(t-h))+1
  upperl=floor(n*(t+h))
  
  istar <- max(1,lowerl):min(n,upperl)
  
  y.star <- y[istar]
  wstar <- as.numeric(wts[istar])
  
  # s1 is score of beta0
  # s2 is score of beta1
  # s3 is score of gamma0
  # s4 is score of delta0
  k1 <- (1+delta0*(y.star-beta0-beta1*(istar-n*t))/exp(gamma0))^(-1/delta0)
  score.num <- y.star-beta0-beta1*(istar-n*t)
  score.den <- exp(gamma0)+delta0*(y.star-beta0-beta1*(istar-n*t))
  d_k1_d_delta0 <- -k1/delta0*(log(k1)+score.num/score.den)
  s1 <- (delta0+1-k1)/score.den
  s2 <- (delta0+1-k1)*(istar-n*t)/score.den
  s3 <- (delta0+1-k1)*score.num/score.den-1
#  s4 <- log(k1)*(1-(delta0+1-k1)/delta0)-(delta0+1-k1)/delta0*(score.num/score.den)
  s4 <- log(k1)-((delta0+1-k1)/delta0)*(log(k1)+score.num/score.den)
  
  s11 <- sum(wstar^2*s1*s1,na.rm=TRUE)
  s12 <- sum(wstar^2*s1*s2,na.rm=TRUE)
  s13 <- sum(wstar^2*s1*s3,na.rm=TRUE)
  s14 <- sum(wstar^2*s1*s4,na.rm=TRUE)
  s22 <- sum(wstar^2*s2*s2,na.rm=TRUE)
  s23 <- sum(wstar^2*s2*s3,na.rm=TRUE)
  s24 <- sum(wstar^2*s2*s4,na.rm=TRUE)
  s33 <- sum(wstar^2*s3*s3,na.rm=TRUE)
  s34 <- sum(wstar^2*s3*s4,na.rm=TRUE)
  s44 <- sum(wstar^2*s4*s4,na.rm=TRUE)
  
  score.mat <- matrix(c(s11,s12,s13,s14,
                s12,s22,s23,s24,
                s13,s23,s33,s34,
                s14,s24,s34,s44),ncol=4,byrow=T)
  
  # Elements of J matrix
  J11 <- -sum(wstar*(delta0+1)*(delta0-k1)/score.den^2,na.rm=TRUE)
  J12 <- -sum(wstar*(delta0+1)*(delta0-k1)*(istar-n*t)/score.den^2,na.rm=TRUE)
  J13 <- -sum(wstar*(-k1*score.num-exp(gamma0)*(delta0+1-k1))/score.den^2,na.rm=TRUE)
  J14 <- -sum(wstar*(score.den*(1-d_k1_d_delta0)-(delta0+1-k1)*score.num)/score.den^2,na.rm=TRUE)
  J22 <- -sum(wstar*(delta0+1)*(delta0-k1)*(istar-n*t)^2/score.den^2,na.rm=TRUE)
  J23 <- -sum(wstar*(-k1*(istar-n*t)*score.num-exp(gamma0)*(istar-n*t)*(delta0+1-k1))/score.den^2,na.rm=TRUE)
  J24 <- -sum(wstar*(score.den*(istar-n*t)*(1-d_k1_d_delta0)-(istar-n*t)*(delta0+1-k1)*score.num)/score.den^2,na.rm=TRUE)
  J33 <- -sum(wstar*(-k1*score.num^2-(delta0+1-k1)*score.num*exp(gamma0))/score.den^2,na.rm=TRUE)
  J34 <- -sum(wstar*(score.den*score.num*(1-d_k1_d_delta0)-(delta0+1-k1)*score.num^2)/score.den^2,na.rm=TRUE)
#  J44 <- -sum(wstar*(d_k1_d_delta0/k1*((k1-1)/delta0)+log(k1)*(delta0*d_k1_d_delta0+1-k1)/delta0^2+(delta0*d_k1_d_delta0+1-k1)/delta0^2*score.num/score.den+((k1-1)/delta0-1)*(-(score.num/score.den)^2)))
  J44 <- -sum(wstar*(d_k1_d_delta0/k1*(1+k1/delta0)*(1-k1)+log(k1)*(delta0*d_k1_d_delta0+1-k1)/delta0^2+(delta0*d_k1_d_delta0+1-k1)/delta0^2*score.num/score.den+((k1-1)/delta0-1)*(-(score.num/score.den)^2)),na.rm=TRUE)
  
  J <- matrix(c(J11,J12,J13,J14,
                J12,J22,J23,J24,
                J13,J23,J33,J34,
                J14,J24,J34,J44),ncol=4,byrow=T)
  
  V <- solve(J)%*%score.mat%*%solve(J)
  V
}
