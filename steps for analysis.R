# This analysis uses h=0.3.

library(tidyverse)

# Read in the Congaree data
y <- peak$Peak_Flow

# Run the robustlikelihood estimation
t <- (1:131)/131
congareeRLL <- RLLgev(t,y,0.35,50000)

# Save the parameter estimates (for data visualization purposes)
congaree.est <- as.data.frame(congareeRLL$parm.est) %>%
  mutate(Year=1892:2022)

# 1%-chance flood standard error estimation
Q99.SE <- NULL
for(i in 1:131){
  Q99.SE[i] <- RLLqsd(congareeRLL$parm.est[i,1],y,congareeRLL$wts[,i],0.35,50000,
                      congareeRLL$parm.est[i,2:5])
}

congaree.est <- congaree.est %>%
  mutate(Q99.SE=Q99.SE)

# Boxplots of weights by flood year
congareeRLL.wts <- cbind(congareeRLL$wts,row.num=1:131)
congareeRLL.wts.long <- as.data.frame(congareeRLL.wts) %>%
  pivot_longer(-row.num) %>%
  mutate(Year=as.numeric(substr(name,2,4))+1891) %>%
  select(-name,-row.num) %>%
  filter(!is.na(value))
ggplot(congareeRLL.wts.long,aes(x=factor(Year),y=value))+geom_boxplot()+
  xlab("Flood Year")+ylab("Robust Local Likelihood Weights")+theme_classic()+
  scale_x_discrete(breaks=c("1900","1920","1940","1960","1980","2000","2020"))

# Check condition number of sandwich estimator
cond.num <- NULL
for(i in 1:131){
  A <- RLL.sandwich.mat(t[i],y,congareeRLL$wts[,i],0.35,50000,
                        congareeRLL$parm.est[i,2:5])
  cond.num[i] <- kappa(A,exact=TRUE)
}

# Add condition number of sandwich estimator to Congaree data
congaree.est <- congaree.est %>%
  mutate(cond.num=cond.num)

# The n for each flood year's robust likelihood estimation
n.robust <- NULL
for(i in 1:131){
  n.robust[i] <- sum(!is.na(congareeRLL.wts[,i]))
}

congaree.est <- congaree.est %>%
  mutate(n.robust=n.robust)
