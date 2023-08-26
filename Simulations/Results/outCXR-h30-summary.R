outCXR.h30 <- read.csv("~/Dropbox/Phil Yates/Sim Results/outCXR-h30.csv")

library(tidyverse)
outCXR.h30.summary <- outCXR.h30 %>%
group_by(tseq) %>%
summarize(var.Median=var(Median),var.Mean=var(Mean),var.Q99=var(Q99),mean.Q99=mean(Q99),mean.Median=mean(Median),mean.mean=mean(Mean))

library(FAdist)
true.Q99 <- NULL
true.Median <- NULL
true.Mean <- NULL

# Calculate eta.C, tau.C, kappa.C
tseq <- outCXR.h30.summary$tseq
eta.C <- max(0.25*tseq+0.1*dnorm(tseq,0.5,0.1)+0.1*dnorm(tseq,1,0.1))
tau.C <- max(exp(-1.5+0.2*tseq))
kappa.C <- -0.1

# Calculate the 0.99-quantile for Model C (the "contaminated" version of Model B)
cmmgev <- function(x,p,pi,eta,tau,kappa){
  cmmgev <- pi*pgev(x,location=eta[1],scale=tau[1],shape=kappa[1])+
    (1-pi)*pgev(x,location=eta[2],scale=tau[2],shape=kappa[2])-p
  cmmgev
}

qmixgev <- function(p,pi,eta,tau,kappa){
  q1 <- qgev(p,location=eta[1],scale=tau[1],shape=kappa[1])
  q2 <- qgev(p,location=eta[2],scale=tau[2],shape=kappa[2])
  qa <- min(q1,q2)
  qb <- max(q1,q2)
  qmm <- uniroot(cmmgev,interval=c(qa,qb),p=p,pi=pi,eta=eta,tau=tau,kappa=kappa)$root
  qmm
}


for(i in 1:100){
  t <- outCXR.h30.summary$tseq[i]
  eta.B <- 0.25*t+0.1*dnorm(t,0.5,0.1)+0.1*dnorm(t,1,0.1)
  tau.B <- exp(-1.5+0.2*t)
  kappa.B <- -0.1
  if(i<=99){
  true.Q99[i] <- qmixgev(0.99,0.9,eta=c(eta.B,eta.C),tau=c(tau.B,tau.C),
                         kappa=c(kappa.B,kappa.C))
  true.Median[i] <- qmixgev(0.5,0.9,eta=c(eta.B,eta.C),tau=c(tau.B,tau.C),
                            kappa=c(kappa.B,kappa.C))
  true.Mean[i] <- 0.9*(eta.B+tau.B*(gamma(1-kappa.B)-1)/kappa.B)+0.1*(eta.C+tau.C*(gamma(1-kappa.C)-1)/kappa.C)
  }
  if(i>99){
    true.Q99[i] <- qgev(0.99,location=eta.B,scale=tau.B,shape=kappa.B)
    true.Median[i] <- eta.B+tau.B*(log(2)^(-kappa.B)-1)/kappa.B
    true.Mean[i] <- eta.B+tau.B*(gamma(1-kappa.B)-1)/kappa.B
  }
}

outCXR.h30.summary <- outCXR.h30.summary %>%
mutate(bias.Q99=mean.Q99-true.Q99,MSE.Q99=var.Q99+bias.Q99^2,bias.Median=mean.Median-true.Median,MSE.Median=var.Median+bias.Median^2,bias.Mean=mean.mean-true.Mean,MSE.Mean=var.Mean+bias.Mean^2)
write_csv(outCXR.h30.summary,"/Users/philipayates/Dropbox/Phil Yates/Sim Results/Summaries/outCXR-h30.csv")


