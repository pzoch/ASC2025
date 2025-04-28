rm(list = ls())
require(zoo)
require(ggfortify)
require(lubridate)
require(urca)
require(knitr)
require(stringr)
require(data.table)
require(dynlm)
require(lmtest)
# load data
input_dir  =  "C:/users/Piotr/Dropbox/WNE/ASC_2025/data/"

data_pl     = read.csv(file = paste0(input_dir, "GDP_POLAND.csv"))
colnames(data_pl) = c("date","Y")

data_pl$date = strptime(data_pl$date, "%m/%d/%Y" )
data_pl$Y    = zoo(data_pl$Y,order.by = data_pl$date)
data_pl$date = as.POSIXct(data_pl$date)

Y         =  data_pl$Y
Y_short   =  window(data_pl$Y, end = parse_date_time("01/01/2015", orders = "%m-%d-%Y"))

# transform to logs
logY          = log(Y)

logY <- ts(logY,start = c(1995, 1), frequency =4)
dlogY <- (diff(logY, lag = 1) * 100) # growth rate

dlogY   <- window(dlogY,start = c(2004, 4), end = c(2019, 4))
data_test.dlogY  <- window(dlogY, start = c(2019, 1))
data_training.dlogY    <- window(dlogY, end = c(2018, 4))

autoplot(data_training.dlogY)


# Box-Jenkins procedure

  # ACF and PACF
  par(mfrow=c(1,2))
  acf(data_training.dlogY,lag.max=8) 
  pacf(data_training.dlogY,lag.max=8)  

  # is there a problem?

  model  <- arima(data_training.dlogY, order = c(4,0,4), method = "ML")
  
  # more specifications
  model_plusar  <- arima(data_training.dlogY, order = c(4,0,0), method = "ML")
  model_plusma  <- arima(data_training.dlogY, order = c(3,0,1), method = "ML")
  
  coeftest(model)
  coeftest(model_plusar)
  coeftest(model_plusma)

  resids <- model$residuals
  Box.test(resids,lag=8, type="Ljung-Box",fitdf = 3 + 0 )
  checkresiduals(model)
  
  
  # how do we interpret these results?
  par(mfrow=c(1,2))
  acf(resids,lag.max=8) 
  pacf(resids,lag.max=8)  
  coeftest(model)
  
  # try with AIC etc
  # information criteria
  model_aic = auto.arima(data_training.dlogY,max.d=0,max.p=8,max.q=8,ic="aic",seasonal = FALSE)
  model_bic = auto.arima(data_training.dlogY,max.d=0,max.p=8,max.q=8,ic="bic",seasonal = FALSE)

  checkresiduals(model_aic)
  checkresiduals(model_bic)
  
  
  
  # general to specific modeling
  model_general  <- arima(data_training.dlogY, order = c(4,0,4), method = "ML")
  model_specific  <- arima(data_training.dlogY, order = c(2,0,3), method = "ML")
  lgen = logLik(model_general)
  lspec =logLik(model_specific)
  dfs <- length(coef(model_general)) - length(coef(model_specific))
  teststat<--2*(as.numeric(lspec-lgen))
  pchisq(teststat,df=dfs,lower.tail=FALSE)
  
  checkresiduals(model_specific)
  ### conclusion? ###

  
ar.coef <- model_specific[["coef"]][1:2]
ma.coef <- model_specific[["coef"]][3:5]

#MA representation
par(mfrow=c(1,1))
IRF <- ARMAtoMA(ar = ar.coef,ma = ma.coef,12)
IRF <- ts(c(1,IRF),start=0)
plot.ts(IRF)


# forecasts
autoplot(forecast(model_specific,h=4))

