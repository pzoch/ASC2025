rm(list=ls())
require(ggplot2)
require(zoo)
require(urca)
require(fpp2)

# generate AR(2)
v_rho  <- c(1,0)
sigma2 <- 1

# simulation
set.seed(220)
nobs      = 100
x    = matrix(0,nobs,1)
x[1] = rnorm(1,mean = 0, sd = sqrt(sigma2))
x[2] = v_rho[1]*x[1] +  rnorm(1,mean = 0, sd = sqrt(sigma2))

  for(t in 3:nobs){
    x[t]=v_rho[1]*x[t-1]+v_rho[2]*x[t-2] + rnorm(1,mean = 0, sd = sqrt(sigma2))
  }


ts.plot(x)

# Test DF
########################################
dx <- diff(x) 

ts.plot(dx)
df_test1 <- ur.df(x, type="none", lags=0)
ts.plot(df_test1@res)

ggAcf(df_test1@res)
ggPacf(df_test1@res)


df_test2 <- ur.df(dx, type="none", lags=0)
df_test3 <- ur.df(x, type="drift", lags=0)
df_test4 <- ur.df(x, type="trend", lags=0)

summary(df_test1)
summary(df_test2)
summary(df_test3)
summary(df_test4)


# Test KPSS
########################################

kpss_test1 <- ur.kpss(x, type = "mu")
kpss_test2 <- ur.kpss(x, type = "tau")
summary(kpss_test1)
summary(kpss_test2)


# Power and size of UR tests
########################################


# process y(t) = rho * y(t-1) + (1-rho) * alpha + e(t) with e(t) ~ N(0,sigma2)
set.seed(777)

T    <- 4 * 25   # quarters: Polish data 4 * 25 (years), US data 4 * 70 (years)
N    <- 100   # experiments
rho    <- 0.5  # rho 
sigma2 <- 1.0  # variance
alpha  <- 1.0  # constant

count <- matrix(NA,N,2)
for (n in 1:N) {
  y<- e <- rnorm(T, mean = 0, sd = sqrt(sigma2))
  # start from unconditional mean
  y[1] <- alpha + 1/sqrt((1-rho^2))*e[1]
  for(i in 2:T) y[i] <- (1-rho) * alpha + rho*y[i-1] + e[i]
  
  z<-ur.df(y, type="drift", lags = 0)
  count[n,1] <- z@teststat[1] < z@cval[1,2] 
  z<-ur.kpss(y, type = "mu")
  count[n,2] <- z@teststat < z@cval[2]
}
Table        <- colMeans(count)*100
names(Table) <- c("ADF","KPSS")
Table

