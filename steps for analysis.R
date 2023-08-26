# This analysis uses h=0.3.

library(tidyverse)

# Read in the Congaree data
y <- peak$Peak_Flow

# Run the robustlikelihood estimation
t <- (1:131)/131
congareeRLL <- RLLgev(t,y,0.3,50000)

# Save the parameter estimates (for data visualization purposes)
congaree.est <- as.data.frame(congareeRLL$parm.est) %>%
  mutate(Year=1892:2022)

# 1%-chance flood standard error estimation
Q99.SE <- NULL
for(i in 1:131){
  Q99.SE[i] <- RLLqsd(congareeRLL$parm.est[i,1],y,congareeRLL$wts[,i],0.3,50000,
                      congareeRLL$parm.est[i,2:5])
}

congaree.est <- congaree.est %>%
  mutate(Q99.SE=Q99.SE)

