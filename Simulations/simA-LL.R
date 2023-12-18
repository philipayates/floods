est.Q100 <- NULL
est.mean <- NULL
est.median <- NULL

# Simulation A, t=0.95, h=0.3
for(i in 2:1001){
  LL.sim <- LLgev(tseq=0.95,yinput=sim.modelA[,i],0.3,1)
  est.Q100 <- c(est.Q100,LL.sim[1,6])
  est.median <- c(est.median,LL.sim[1,7])
  est.mean <- c(est.mean,LL.sim[1,8])
}

t <- 0.95
eta.A <- 0.25*t
tau.A <- exp(-1.5+0.2*t)
kappa.A <- rep(-0.1,length(t))

true.Q100 <- qgev(0.99,loc=eta.A,scale=tau.A,shape=kappa.A)
true.median <- eta.A+tau.A*(log(2)^(-kappa.A)-1)/kappa.A
true.mean <- eta.A+tau.A*(gamma(1-kappa.A)-1)/kappa.A

# Mean
MSE.mean <- var(est.mean)+(mean(est.mean)-true.mean)^2
bias.mean <- mean(est.mean)-true.mean

# Median 
MSE.median <- var(est.median)+(mean(est.median)-true.median)^2
bias.median <- mean(est.median)-true.median

# Q100
MSE.Q100 <- var(est.Q100)+(mean(est.Q100)-true.Q100)^2
bias.Q100 <- mean(est.Q100)-true.Q100