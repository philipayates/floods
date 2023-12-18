sandwich_stationary.mat <- function(yinput,scale.factor,parm){
  y.star <- yinput/scale.factor
  mu0 <- as.numeric(parm[1])
  sigma0 <- as.numeric(parm[2])
  delta0 <- as.numeric(parm[3])
  n <- length(y.star)
  
  # s1 is score of mu0
  # s2 is score of sigma0
  # s3 is score of delta0
  k1 <- (1+delta0*(y.star-mu0)/sigma0)^(-1/delta0)
  score.num <- y.star-mu0
  score.den <- sigma0+delta0*(y.star-mu0)
  d_k1_d_delta0 <- -k1/delta0*(log(k1)+score.num/score.den)
  s1 <- (delta0+1-k1)/score.den
  s2 <- -1/sigma0 + ((delta0+1-k1)*score.num)/(sigma0*score.den)
  s3 <- log(k1)-((delta0+1-k1)/delta0)*(log(k1)+score.num/score.den)
  
  s11 <- sum(s1*s1)
  s12 <- sum(s1*s2)
  s13 <- sum(s1*s3)
  s22 <- sum(s2*s2)
  s23 <- sum(s2*s3)
  s33 <- sum(s3*s3)
  
  score.mat <- matrix(c(s11,s12,s13,
                s12,s22,s23,
                s13,s23,s33),ncol=3,byrow=T)
  
  # Elements of J matrix
  J11 <- -sum((delta0+1)*(delta0-k1)/score.den^2)
  J12 <- -sum((-score.num*k1/sigma0 - (delta0+1-k1))/score.den^2)
  J13 <- -sum((score.den*(1-d_k1_d_delta0)-(delta0+1-k1)*score.num)/score.den^2)
  J22 <- -sum(sigma0^(-2) - k1*score.num^2/(sigma0*score.den)^2 - score.num*(delta0+1-k1)*(sigma0+score.den)/(sigma0*score.den)^2)
  J23 <- -sum(score.num*(1-d_k1_d_delta0)/(sigma0*score.den) - sigma0*score.num^2*(delta0+1-k1)/(sigma0*score.den)^2)
  J33 <- -sum((d_k1_d_delta0/k1*(1+k1/delta0)*(1-k1)+log(k1)*(delta0*d_k1_d_delta0+1-k1)/delta0^2+(delta0*d_k1_d_delta0+1-k1)/delta0^2*score.num/score.den+((k1-1)/delta0-1)*(-(score.num/score.den)^2)))
  J <- matrix(c(J11,J12,J13,
                J12,J22,J23,
                J13,J23,J33),ncol=3,byrow=T)
  V <- solve(J)%*%score.mat%*%solve(J)
  #browser()
  V
}
