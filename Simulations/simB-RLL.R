est.Q100 <- NULL
est.mean <- NULL
est.median <- NULL

# Simulation B, t=0.95, h=0.3
for(i in 2:1001){
  RLL.sim <- RLLgev(tseq=0.95,yinput=sim.modelB[,i],0.3,1)$parm.est
  est.Q100 <- c(est.Q100,RLL.sim[1,6])
  est.median <- c(est.median,RLL.sim[1,7])
  est.mean <- c(est.mean,RLL.sim[1,8])
}

t <- 0.95
eta.B <- 0.25*t+0.1*dnorm(t,0.5,0.1)+0.1*dnorm(t,1,0.1)
tau.B <- exp(-1.5+0.2*t)
kappa.B <- rep(-0.1,length(t))

true.Q100 <- qgev(0.99,loc=eta.B,scale=tau.B,shape=kappa.B)
true.median <- eta.B+tau.B*(log(2)^(-kappa.B)-1)/kappa.B
true.mean <- eta.B+tau.B*(gamma(1-kappa.B)-1)/kappa.B

# Mean
MSE.mean <- var(est.mean)+(mean(est.mean)-true.mean)^2
bias.mean <- mean(est.mean)-true.mean

# Median 
MSE.median <- var(est.median)+(mean(est.median)-true.median)^2
bias.median <- mean(est.median)-true.median

# Q100
MSE.Q100 <- var(est.Q100)+(mean(est.Q100)-true.Q100)^2
bias.Q100 <- mean(est.Q100)-true.Q100