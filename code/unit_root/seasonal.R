rm(list=ls())
set.seed(13) 

require(data.table)
require(ggplot2)
require(dyn)
require(zoo)
require(urca)
require(uroot)
require(fpp2)
require(lmtest)
require(eurostat)
require(lubridate)


no_periods = 20 * 4
scale = 0.1

# simulate deterministic seasonality 
white_noise <- arima.sim(model = list(order = c(0, 0, 0)), n = no_periods, sd = scale)
x_det <- ts(rep(c(-1,0,1,2),no_periods / 4)+white_noise,frequency = 4)


# simulate ar(1) seasonality 
model <- Arima(ts(rnorm(no_periods, sd = scale),freq=4), order=c(0,0,0), seasonal=c(1,0,0),fixed=c(Phi_1=0.5,mu=0))
x_ar <- simulate(model, no_periods=no_periods)

# simulate unit root seasonality 
model <- Arima(ts(rnorm(no_periods, sd = scale),freq=4), order=c(0,0,0), seasonal=c(0,1,0))
x_ur <- simulate(model, no_periods=no_periods)



x_det <- ts(x_det,start = c(1980, 1), frequency =4)
x_ar  <- ts(x_ar,start = c(1980, 1), frequency =4) 
x_ur  <- ts(x_ur,start = c(1980, 1), frequency =4)

horizon = 12
# generate dummies
d1 = ts(rep(c(1,0,0,0),no_periods / 4 + horizon /4),start = c(1980, 1), frequency =4)
d2 = ts(rep(c(0,1,0,0),no_periods / 4 + horizon /4),start = c(1980, 1), frequency =4)
d3 = ts(rep(c(0,0,1,0),no_periods / 4 + horizon /4),start = c(1980, 1), frequency =4)

dummies     = cbind(d1,d2,d3)
dummies_reg = window(dummies,start = c(1980, 1),end = c(1999, 4))
dummies_fcast = window(dummies,start = c(2000, 1))


#model with deterministic seasonality
model_det  <- Arima(x_det, order = c(0,0,0), xreg = dummies_reg,method = "ML",)
autoplot(forecast(model_det,h=length(dummies_fcast),xreg  = dummies_fcast))


#model with unit root seasonality
model_ur  <- Arima(x_ur, order=c(0,0,0), seasonal=c(0,1,0))
autoplot(forecast(model_ur),h=horizon)



#distinguishing seasonality?
autoplot(x_det) + autolayer(x_ar) + autolayer(x_ur)

autoplot(x_ar) + autolayer(x_ur)

par(mfrow=c(1,3))



spec.pgram(coredata(x_det),log='no')
spec.pgram(coredata(x_ar),log='no')
spec.pgram(coredata(x_ur),log='no')

par(mfrow=c(2,3))
acf(x_det,lag.max=24) 
acf(x_ar,lag.max=24) 
acf(x_ur,lag.max=24) 
pacf(x_det,lag.max=24)  
pacf(x_ar,lag.max=24)  
pacf(x_ur,lag.max=24)  



# Box-Jenkins AirPassengers

autoplot(AirPassengers)

par(mfrow=c(1,3))
plot(window(AirPassengers,start = c(1951, 1),end = c(1951, 12)),type='l')
plot(window(AirPassengers,start = c(1953, 1),end = c(1953, 12)),type='l')
plot(window(AirPassengers,start = c(1955, 1),end = c(1955, 12)),type='l')

lAP <- log(AirPassengers)
autoplot(lAP)
par(mfrow=c(1,3))
plot(window(lAP,start = c(1951, 1),end = c(1951, 12)),type='l')
plot(window(lAP,start = c(1953, 1),end = c(1953, 12)),type='l')
plot(window(lAP,start = c(1955, 1),end = c(1955, 12)),type='l')

par(mfrow=c(1,3))
acf(lAP,lag.max=24) 
pacf(lAP,lag.max=24) 
spec.pgram(coredata(lAP),log='no',span=5)

df_test <-  ur.df(lAP, type="none", lags=0)

summary(df_test)

df_test <-  ur.df(lAP, type="trend", lags=12)
summary(df_test)

dlAP12 <- diff(lAP,lag=12)
autoplot(dlAP12)

par(mfrow=c(1,3))
acf(dlAP12,lag.max=24) 
pacf(dlAP12,lag.max=24) 
spec.pgram(coredata(dlAP12),log='no',span=5)

df_test <-  ur.df(dlAP12, type="none", lags=12)
summary(df_test)


HEGY_test <- hegy.test(lAP, deterministic = c(1,1,1),lag.method = c("AIC"), maxlag = 24)
CH_test <- ch.test(lAP,  type ="trigonometric", sid ="all",lag1=TRUE)

dlAP <- diff(lAP)
autoplot(dlAP)
par(mfrow=c(1,3))
acf(dlAP,lag.max=24) 
pacf(dlAP,lag.max=24) 
spec.pgram(coredata(dlAP),log='no',span=5)

HEGY_test <- hegy.test(dlAP, deterministic = c(1,1,1),lag.method = c("AIC"), maxlag = 24)
CH_test <- ch.test(dlAP,  type ="trigonometric", sid ="all",lag1=TRUE)


df_test3 <- ur.df(dlAP, type="drift", lags=0)
df_test4 <- ur.df(dlAP, type="trend", lags=0)


dSlAP <- diff(lAP,12)
autoplot(dSlAP)
par(mfrow=c(1,3))
acf(dSlAP,lag.max=24) 
pacf(dSlAP,lag.max=24) 
spec.pgram(coredata(dSlAP),log='no',span=5)

HEGY_test <- hegy.test(dSlAP, deterministic = c(1,1,1),lag.method = c("AIC"), maxlag = 24)
CH_test <- ch.test(dSlAP,  type ="trigonometric", sid ="all",lag1=TRUE)

df_test5 <- ur.df(dSlAP, type="trend", lags=12)

ddSlAP <- diff(dlAP,12)
autoplot(ddSlAP)
par(mfrow=c(1,3))
acf(ddSlAP,lag.max=24) 
pacf(ddSlAP,lag.max=24) 
spec.pgram(coredata(ddSlAP),log='no',span=5)

HEGY_test <- hegy.test(ddSlAP, deterministic = c(1,1,1),lag.method = c("AIC"), maxlag = 24)
CH_test <- ch.test(ddSlAP,  type ="trigonometric", sid ="all",lag1=TRUE)

model_s  <- arima(lAP, order = c(1,1,1), seasonal=c(1,1,1), method = "ML")
res_model_s<- residuals(model_s)
Box.test(res_model_s, lag = 1,  type = "Ljung")
Box.test(res_model_s, lag = 12, type = "Ljung")
checkresiduals(model_s)


model_auto  <- auto.arima(lAP, method = "ML",max.d=1,max.p=24,max.q=24,ic="aic")
checkresiduals(model_auto)
autoplot(forecast(model_auto))


         