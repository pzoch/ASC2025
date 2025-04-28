rm(list=ls())
set.seed(13) 

require(data.table)
require(ggplot2)
require(dyn)
require(zoo)
require(urca)
require(fpp2)
require(lmtest)
require(eurostat)
require(lubridate)

# load data

datasetY <- get_eurostat("namq_10_gdp")

logY <-subset(datasetY, geo == c("PL") & na_item =="B1G" & unit=="CP_MNAC" & s_adj == "NSA") 
logY <- na.omit(ts(rev(logY$values),start = c(1995, 1), frequency =4))
dlogY <- (diff(logY, lag = 1) * 100) # growth rate


dlogY   <- window(dlogY,start = c(2004, 4), end = c(2019, 4))
data_testing.dlogY  <- window(dlogY, start = c(2019, 1))
data_training.dlogY    <- window(dlogY, end = c(2018, 4))

autoplot(dlogY)



#SARMA
# ACF and PACF
par(mfrow=c(1,2))
acf(data_training.dlogY,lag.max=12) 
pacf(data_training.dlogY,lag.max=12)  


data_training.dlogY4 <- diff(data_training.dlogY,lag=4)
par(mfrow=c(1,2))
acf(data_training.dlogY4,lag.max=12) 
pacf(data_training.dlogY4,lag.max=12)  


model_s  <- arima(data_training.dlogY, order = c(1,0,1), seasonal=c(1,0,1), method = "ML")
resids_s <- model_s$residuals
checkresiduals(model_s)
# how do we interpret these results?
par(mfrow=c(1,2))
acf(resids_s,lag.max=24) 
pacf(resids_s,lag.max=24)  
coeftest(model_s)

model_s_aic = auto.arima(data_training.dlogY,max.d=0,max.D=0,max.p=4,max.q=4,ic="aic")
model_s_bic = auto.arima(data_training.dlogY,max.d=0,max.D=0,max.p=4,max.q=4,ic="bic")
checkresiduals(model_s_aic)
checkresiduals(model_s_bic)
autoplot(forecast(model_s,h=4)) + autolayer(dlogY)
autoplot(forecast(model_s_aic,h=4)) + autolayer(dlogY)
autoplot(forecast(model_s_bic,h=4)) + autolayer(dlogY)

model_s_aic = auto.arima(data_training.dlogY,max.d=0,max.D=1,max.p=4,max.q=4,ic="aic")
model_s_bic = auto.arima(data_training.dlogY,max.d=0,max.D=1,max.p=4,max.q=4,ic="bic")
checkresiduals(model_s_aic)
checkresiduals(model_s_bic)
autoplot(forecast(model_s,h=4)) + autolayer(dlogY)
autoplot(forecast(model_s_aic,h=4)) + autolayer(dlogY)
autoplot(forecast(model_s_bic,h=4)) + autolayer(dlogY)


model  <- arima(data_training.dlogY, order = c(1,0,1), method = "ML")
autoplot(forecast(model,h=4)) + autolayer(dlogY)
accuracy(model_s)
