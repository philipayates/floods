outAX.h15 <- read.csv("~/Dropbox/Phil Yates/Sim Results/outAXR-h30.csv")

library(tidyverse)
outAX.h15.summary <- outAX.h15 %>%
group_by(tseq) %>%
summarize(var.Median=var(Median),var.Mean=var(Mean),var.Q99=var(Q99),mean.Q99=mean(Q99),mean.Median=mean(Median),mean.mean=mean(Mean))

library(FAdist)
true.Q99 <- NULL
true.Median <- NULL
true.Mean <- NULL

for(i in 1:100){
  t <- outAX.h15.summary$tseq[i]
  eta.A <- 0.25*t
  tau.A <- exp(-1.5+0.2*t)
  kappa.A <- -0.1
  true.Q99[i] <- qgev(0.99,location=eta.A,scale=tau.A,shape=kappa.A)
  true.Median[i] <- eta.A+tau.A*(log(2)^(-kappa.A)-1)/kappa.A
  true.Mean[i] <- eta.A+tau.A*(gamma(1-kappa.A)-1)/kappa.A
}

outAX.h15.summary <- outAX.h15.summary %>%
mutate(bias.Q99=mean.Q99-true.Q99,MSE.Q99=var.Q99+bias.Q99^2,bias.Median=mean.Median-true.Median,MSE.Median=var.Median+bias.Median^2,bias.Mean=mean.mean-true.Mean,MSE.Mean=var.Mean+bias.Mean^2)
write_csv(outAX.h15.summary,"/Users/philipayates/Dropbox/Phil Yates/Sim Results/Summaries/outAXR-h30.csv")


