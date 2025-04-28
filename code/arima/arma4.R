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

datasetPi <- get_eurostat("prc_hicp_manr")

tempPi <- subset(datasetPi, geo == c("PL") & coicop=="CP00" & unit=="RCH_A") 
tempPi <- na.omit(ts(rev(tempPi$values),start = c(1997, 1), frequency =12))


Pi   <- window(tempPi,start = c(2002, 1), end = c(2019, 12))


data_training.Pi <- window(Pi, end = c(2018, 12))
data_testing.Pi <- window(Pi, start = c(2019, 1))
autoplot(Pi)



# Box-Jenkins procedure

  # ACF and PACF
  par(mfrow=c(1,2))
  acf(data_training.Pi,lag.max=24) 
  pacf(data_training.Pi,lag.max=24)
  df_test <- ur.df(data_training.Pi, type="none", lags=4)
  
  
  
  # is there a problem?
  model  <- arima(data_training.Pi, order = c(1,1,2), method = "ML")
  
  data_training.dPi = diff(data_training.Pi)
  par(mfrow=c(1,2))
  acf(data_training.dPi,lag.max=24) 
  pacf(data_training.dPi,lag.max=24)
  
  # more specifications
  model_plusar  <- arima(data_training.Pi, order = c(2,1,2), method = "ML")
  model_plusma  <- arima(data_training.Pi, order = c(1,1,3), method = "ML")
  coeftest(model)
  coeftest(model_plusar)
  coeftest(model_plusma)
  
  resids <- model$residuals
  Box.test(resids,lag=24, type="Ljung-Box",fitdf = 2 + 1 )
  checkresiduals(model)
  
  
  # how do we interpret these results?
  par(mfrow=c(1,2))
  acf(resids,lag.max=24) 
  pacf(resids,lag.max=24)  
  coeftest(model)
  
# try with AIC etc
  # information criteria
  model_aic = auto.arima(data_training.Pi,max.d=1,max.p=24,max.q=24,ic="aic",seasonal = FALSE,stepwise=FALSE, approximation=FALSE)
  model_bic = auto.arima(data_training.Pi,max.d=1,max.p=24,max.q=24,ic="bic",seasonal = FALSE,stepwise=FALSE, approximation=FALSE)

  checkresiduals(model_aic)
  checkresiduals(model_bic)
  
  
  
# general to specific 
  model_general   <- arima(data_training.Pi, order = c(1,1,2), method = "ML")
  model_specific  <- arima(data_training.Pi, order = c(1,1,0), method = "ML")
  lgen = logLik(model_general)
  lspec =logLik(model_specific)
  dfs <- length(coef(model_general)) - length(coef(model_specific))
  teststat<--2*(as.numeric(lspec-lgen))
  pchisq(teststat,df=dfs,lower.tail=FALSE)
  
  
  
  
  
  ### conclusion? ###


ar.coef <- model[["coef"]][1:2]
ma.coef = NULL


#MA representation
par(mfrow=c(1,1))
IRF <- ARMAtoMA(ar = ar.coef,ma = ma.coef,24)
IRF <- ts(c(1,IRF),start=0)
plot.ts(IRF)


# forecasts
autoplot(forecast(model,h=12)) + autolayer(Pi)
autoplot(forecast(model_aic,h=12)) + autolayer(Pi)
autoplot(forecast(model_bic,h=12)) + autolayer(Pi)


###SARMA
# ACF and PACF

par(mfrow=c(1,2))
acf(data_training.Pi,lag.max=24) 
pacf(data_training.Pi,lag.max=24)  


data_training.Pi12 <- diff(data_training.Pi,lag=12)
par(mfrow=c(1,2))
acf(data_training.Pi12,lag.max=24) 
pacf(data_training.Pi12,lag.max=24)  


model_s  <- arima(data_training.Pi, order = c(2,0,0), seasonal=c(1,1,0), method = "ML")
resids_s <- model_s$residuals

# how do we interpret these results?
par(mfrow=c(1,2))
acf(resids_s,lag.max=24) 
pacf(resids_s,lag.max=24)  
coeftest(model_s)

model_s_aic = auto.arima(data_training.Pi,max.d=0,max.p=24,max.q=24,ic="aic")
model_s_bic = auto.arima(data_training.Pi,max.d=0,max.p=24,max.q=24,ic="bic")
checkresiduals(model_s_aic)
checkresiduals(model_s_bic)
autoplot(forecast(model_s,h=12)) + autolayer(Pi)
autoplot(forecast(model_s_aic,h=12)) + autolayer(Pi)
autoplot(forecast(model_s_bic,h=12)) + autolayer(Pi)



