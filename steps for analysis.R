# This analysis uses h=0.45 for the Congaree River

library(tidyverse)

# Read in the Congaree data
y <- peak$Peak_Flow
t <- (peak$Year-(min(peak$Year)-1))/131

# Run the robust likelihood estimation
congaree.RLL <- RLLgev(t,y,0.45,50000)

# 1%-chance flood standard error estimation
Q100.sd <- NULL
for(i in 1:131){
  Q100.sd[i] <- RLLqsd(t[i],y,congaree.RLL$wts[i,],0.45,50000,
                       congaree.RLL$parm.est[i,2:5])
}

# Save Results
congaree <- as.data.frame(congaree.RLL$parm.est) %>%
  mutate(Q100.sd=Q100.sd,LL.Q100=Q100-qnorm(0.975)*Q100.sd,
         UL.Q100=Q100+qnorm(0.975)*Q100.sd,
         Year=peak$Year,Peak=peak$Peak_Flow) %>%
  mutate(LL.Q100=ifelse(LL.Q100<=0,0,LL.Q100))

# Save the parameter estimates (for data visualization purposes)
congaree.est <- as.data.frame(congareeRLL$parm.est) %>%
  mutate(Year=1892:2022)

# Reshape robust local likelihood estimation weights for boxplots
congaree.wts <- cbind(congaree.RLL$wts,row.num=1:131)
congaree.wts.long <- as.data.frame(congaree.wts) %>%
  pivot_longer(-row.num,values_to="Weight") %>%
  mutate(Year=rep(1892:2022,131)) %>%
  select(-name,-row.num) %>%
  filter(!is.na(Weight))

# Stationary Analysis
library(evd)
mle.parms <- fgev(congareeRLL$Peak/50000)
mle.parms.est <- mle.parms$estimate
Q100.est <- qgev(0.99,loc=mle.parms.est[1]*50000,scale=mle.parms.est[2]*50000,
                 shape=mle.parms.est[3])
Q100.est.SE <- qsd_stationary(congareeRLL$Peak,50000,mle.parms.est)
Q100.LL <- Q100.est-qnorm(0.975)*Q100.est.SE
Q100.UL <- Q100.est+qnorm(0.975)*Q100.est.SE

# Repeat these steps for the Illinois River using h=0.2:
t <- illinois$Index/max(illinois$Index)
y <- illinois$Peak
illinois.RLL <- RLLgev(t,y,0.2,35000)

# Repeat these steps for the Winooski River using h=0.425:
t <- winooski$Index/max(winooski$Index)
y <- winooski$Peak
winooski.RLL <- RLLgev(t,y,0.425,15000)


