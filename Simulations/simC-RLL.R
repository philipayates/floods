library(evd)

est.Q100 <- NULL
est.mean <- NULL
est.median <- NULL

# Simulation C, t=0.5, h=0.3
for(i in 2:1001){
  RLL.sim <- RLLgev(tseq=0.5,yinput=sim.modelC[,i],0.3,1)$parm.est
  est.Q100 <- c(est.Q100,RLL.sim[1,6])
  est.median <- c(est.median,RLL.sim[1,7])
  est.mean <- c(est.mean,RLL.sim[1,8])
}

t <- 0.5
eta.B <- 0.25*t+0.1*dnorm(t,0.5,0.1)+0.1*dnorm(t,1,0.1)
tau.B <- exp(-1.5+0.2*t)
kappa.B <- rep(-0.1,length(t))

ts <- seq(0,1,0.001)
eta.max <- max(0.25*ts+0.1*dnorm(ts,0.5,0.1)+0.1*dnorm(ts,1,0.1))
tau.max <- max(exp(-1.5+0.2*ts))
kappa.max <- -0.1

cmmgev <- function(x,p,pi,eta,tau,kappa){
  cmmgev <- pi*pgev(x,loc=eta[1],scale=tau[1],shape=kappa[1])+
    (1-pi)*pgev(x,loc=eta[2],scale=tau[2],shape=kappa[2])-p
  cmmgev
}

qmixgev <- function(p,pi,eta,tau,kappa){
  q1 <- qgev(p,loc=eta[1],scale=tau[1],shape=kappa[1])
  q2 <- qgev(p,loc=eta[2],scale=tau[2],shape=kappa[2])
  qa <- min(q1,q2)
  qb <- max(q1,q2)
  qmm <- uniroot(cmmgev,interval=c(qa,qb),p=p,pi=pi,eta=eta,tau=tau,kappa=kappa)$root
  qmm
}

true.Q100 <- qmixgev(0.99,0.9,eta=c(eta.B,eta.max),tau=c(tau.B,tau.max),
                     kappa=c(kappa.B,kappa.max))
true.median <- qmixgev(0.5,0.9,eta=c(eta.B,eta.max),tau=c(tau.B,tau.max),
                       kappa=c(kappa.B,kappa.max))
true.mean <- 0.9*(eta.B+tau.B*(gamma(1-kappa.B)-1)/kappa.B)+
  0.1*(eta.max+tau.max*(gamma(1-kappa.max)-1)/kappa.max)

# Mean
MSE.mean <- var(est.mean)+(mean(est.mean)-true.mean)^2
bias.mean <- mean(est.mean)-true.mean

# Median 
MSE.median <- var(est.median)+(mean(est.median)-true.median)^2
bias.median <- mean(est.median)-true.median

# Q100
MSE.Q100 <- var(est.Q100)+(mean(est.Q100)-true.Q100)^2
bias.Q100 <- mean(est.Q100)-true.Q100