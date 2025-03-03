rm(list=ls())
require(ggplot2)
require(dyn)
require(zoo)

####################################
# Random walk with a drift         #
# y[t] = mu + y[t-1] + e[t]        #
#                                  # 
####################################

rho    <- 1
sigma2 <- 1
mu     <- c(0.0,0.0)

# simulation
set.seed(234)
nobs    = 100
RWsim    = matrix(0,nobs,2)
for(k in 1:2){
  for(t in 2:nobs){
    RWsim[t,k]=mu[k] + rho*RWsim[t-1,k]+rnorm(1,mean = 0, sd = sqrt(sigma2))
  }
}


RWsimTab <- data.frame(obs = 1:nobs, RWsim)

# plot time series
ggplot(RWsimTab, aes(x=obs)) + 
  geom_line(aes(y=X1),size=1, col=1)+
  geom_line(aes(y=X2),size=1, col=2)+
  theme(text = element_text(size=12))+
  labs(title="Random walk simulation", y="", x="observation", caption="")+
  geom_hline(yintercept=0, size=1) 

# plot a scatter plot
ggplot(RWsimTab, aes(x=X1, y=X2)) + geom_point() + geom_smooth(method=lm, se=FALSE)

# regression
spurious = lm(RWsimTab$X1 ~ 0 +  RWsimTab$X2)
summary(spurious)

# residuals?
RWsimTab$resid=resid(spurious)
ggplot(RWsimTab, aes(x=obs)) + 
  geom_line(aes(y=resid),size=1)+
  theme(text = element_text(size=12))+
  labs(title="Residuals from OLS", y="", x="observation", caption="")+
  geom_hline(yintercept=0, size=1) 


# test statistic distribution
set.seed(234)
nobs    = 100
nexper  = 5000
RWsim   = matrix(0,nobs,2)
tvec = vector(mode = "numeric",length = nexper)

for (n in 1:nexper){
  for(k in 1:2){
    for(t in 2:nobs){
      RWsim[t,k]=mu[k]+rho*RWsim[t-1,k]+rnorm(1,mean = 0, sd = sqrt(sigma2))
    }
  }
  RWsimTab <- data.frame(obs = 1:nobs, RWsim)
  
  spurious = lm(RWsimTab$X1 ~ 0 +  RWsimTab$X2)
  tvec[n] <- coef(summary(spurious))[, "t value"]
  
  
}

tvector <- data.frame(obs = 1:nexper,tstat = tvec) # returns the density data 

# visualize the results
p <- ggplot(tvector, aes(x=tstat)) + 
  geom_density(alpha=0.25, size = 1.5)

p = p +   stat_function(fun = function(x) dnorm(x, mean = 0, sd = 1),
                        color = "orange", linetype = "dotted", size = 2)

p = p+ geom_vline(aes(xintercept=mean(tstat)),
              color="blue", linetype="dashed", size=2)
p = p+ geom_vline(aes(xintercept=-1.96),
              color="red", linetype="dashed", size=0.5)
p = p+ geom_vline(aes(xintercept=+1.96),
              color="red", linetype="dashed", size=0.5)



p + coord_cartesian(xlim = c(-10,10))

