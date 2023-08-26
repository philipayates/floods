outBXR.h15 <- read.csv("~/Dropbox/Phil Yates/Sim Results/outBXR-h15.csv")

library(tidyverse)
outBXR.h15.summary <- outBXR.h15 %>%
group_by(tseq) %>%
summarize(var.Median=var(Median),var.Mean=var(Mean),var.Q99=var(Q99),mean.Q99=mean(Q99),mean.Median=mean(Median),mean.mean=mean(Mean))

library(FAdist)
true.Q99 <- NULL
true.Median <- NULL
true.Mean <- NULL

for(i in 1:100){
  t <- outBXR.h15.summary$tseq[i]
  eta.B <- 0.25*t+0.1*dnorm(t,0.5,0.1)+0.1*dnorm(t,1,0.1)
  tau.B <- exp(-1.5+0.2*t)
  kappa.B <- -0.1
  true.Q99[i] <- qgev(0.99,location=eta.B,scale=tau.B,shape=kappa.B)
  true.Median[i] <- eta.B+tau.B*(log(2)^(-kappa.B)-1)/kappa.B
  true.Mean[i] <- eta.B+tau.B*(gamma(1-kappa.B)-1)/kappa.B
}

outBXR.h15.summary <- outBXR.h15.summary %>%
mutate(bias.Q99=mean.Q99-true.Q99,MSE.Q99=var.Q99+bias.Q99^2,bias.Median=mean.Median-true.Median,MSE.Median=var.Median+bias.Median^2,bias.Mean=mean.mean-true.Mean,MSE.Mean=var.Mean+bias.Mean^2)
write_csv(outBXR.h15.summary,"/Users/philipayates/Dropbox/Phil Yates/Sim Results/Summaries/outBXR-h15.csv")


