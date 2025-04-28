
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

T = 240
ma_vec = c(0.0,0.0,0.0)
ar_vec = c(0.6,0.2,0.1)
y_arma <-ts(arima.sim(model = list(ma = ma_vec, ar = ar_vec),n=T,sd = 1),frequency=12,start=c(2000,1))
autoplot(y_arma)
par(mfrow=c(1,2))
acf(y_arma,lag.max=24) 
pacf(y_arma,lag.max=24)     


data_training.y_arma <- window(y_arma, end = c(2018, 12))
data_testing.y_arma <- window(y_arma, start = c(2019, 1))



# Box-Jenkins procedure

  # ACF and PACF
  par(mfrow=c(1,2))
  acf(data_training.y_arma,lag.max=24) 
  pacf(data_training.y_arma,lag.max=24)
  df_test <- ur.df(data_training.y_arma, type="none", lags=4)

  # is there a problem?
  model  <- arima(data_training.y_arma, order = c(2,0,0), method = "ML")
  
  # more specifications
  model_plusar  <- arima(data_training.Pi, order = c(3,0,0), method = "ML")
  model_plusma  <- arima(data_training.Pi, order = c(2,0,1), method = "ML")
  coeftest(model)
  coeftest(model_plusar)
  coeftest(model_plusma)
  
  resids <- model$residuals
  Box.test(resids,lag=24, type="Ljung-Box",fitdf = 2 + 0 )
  checkresiduals(model)
  
  
  # how do we interpret these results?
  par(mfrow=c(1,2))
  acf(resids,lag.max=24) 
  pacf(resids,lag.max=24)  
  coeftest(model)
  
# try with AIC etc
  # information criteria
  model_aic = auto.arima(data_training.y_arma,max.d=0,max.p=24,max.q=24,ic="aic",seasonal = FALSE,stepwise=FALSE, approximation=FALSE)
  model_bic = auto.arima(data_training.y_arma,max.d=0,max.p=24,max.q=24,ic="bic",seasonal = FALSE,stepwise=FALSE, approximation=FALSE)

  checkresiduals(model_aic)
  checkresiduals(model_bic)
  
  
  
# general to specific 
  model_general   <- arima(data_training.y_arma, order = c(3,0,1), method = "ML")
  model_specific  <- arima(data_training.y_arma, order = c(2,0,0), method = "ML")
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
autoplot(forecast(model,h=12)) + autolayer(y_arma)
autoplot(forecast(model_aic,h=12)) + autolayer(y_arma)
autoplot(forecast(model_bic,h=12)) + autolayer(y_armai)



