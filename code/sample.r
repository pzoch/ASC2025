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
require(forecast)
# load data
input_dir  =  "D:/Dropbox (Personal)/WNE/ASC_2025/data/"

data     = read.csv(file = paste0(input_dir, "USEPUINDXD.csv"))
colnames(data) = c("date","EPU")

data$date = as.Date(data$date, format = "%Y-%m-%d")
data$EPU    = zoo(data$EPU,order.by = data$date)

data$date = as.POSIXct(data$date)

EPU         =  data$EPU
EPU <- ts(EPU,start = c(1985, 1), frequency =12)
data_test.EPU  <- window(EPU, start = c(2025, 1))

data_training.EPU    <- window(EPU, end = c(2024, 12))

autoplot(data_training.EPU)


# plot acf and pacf
acf(coredata(data_training.EPU), lag.max = 12, plot = TRUE)
pacf(coredata(data_training.EPU), lag.max = 12, plot = TRUE)

# Test ADF
########################################


adf_test <- ur.df(data_training.EPU, type="drift", lags=12)
summary(adf_test)

y <- data_training.EPU
dy <- diff(y)
y_lag1 <- lag(y, -1)[-1]

# Create lagged differences correctly - shift backward in time
dy_lags <- embed(dy, 13)[, -1] 
colnames(dy_lags) <- paste0("dy_lag", 1:12)

dy_clean <- dy[-(1:12)]
y_lag1_clean <- y_lag1[-(1:12)]

# Align all series properly
dy_clean <- dy[-c(1:12)]  # remove first 12 observations
y_lag1_clean <- y_lag1[-c(1:12)]  # remove first 12 observations

# Create the ADF regression: Δy_t = α + βy_{t-1} + Σγ_iΔy_{t-i} + ε_t
adf_regression <- lm(dy_clean ~ y_lag1_clean + dy_lags)
bgtest(adf_regression, order = 24)



adf_test <- ur.df(data_training.EPU, type="drift", lags=12)
summary(adf_test)
# Manually recreate the ADF regression
y <- data_training.EPU
dy <- diff(y)
y_lag1 <- lag(y, -1)[-1]  # lagged level (remove first NA)

dy_lags <- sapply(1:12, function(i) lag(dy, -i))
dy_lags <- dy_lags[-c(1:12), ]  # remove NAs

dy_clean <- dy[-c(1:12)]
y_lag1_clean <- y_lag1[-c(1:12)]

# Create the ADF regression: Δy_t = α + βy_{t-1} + Σγ_iΔy_{t-i} + ε_t
adf_regression <- lm(dy_clean ~ y_lag1_clean + dy_lags)

# Now apply Breusch-Godfrey test
bgtest(adf_regression, order = 24)


summary(ur.kpss(y, type = "mu",lags = "short"))


# conclusion? 


f_hw1    <- hw((y), h=12,"alpha"=0.45,"beta"=0.1)
f_hw2    <- hw((y), h=12,"alpha"=0.25,"beta"=0.1)
f_hw3    <- hw((y), h=12,"alpha"=0.15,"beta"=0.1)
f_hw_opt   <- hw((y), h=12)

autoplot(f_hw1) +
  autolayer(fitted(f_hw1), series="Fitted - Holt-Winters") +
  autolayer(data_test.EPU, series="Data - test") +
  ylab("Inflation") + xlab("Year")


# ACF and PACF
par(mfrow=c(1,2))
acf(y,lag.max=36) 
pacf(y,lag.max=36)  

par(mfrow=c(1,2))
acf(dy,lag.max=36) 
pacf(dy,lag.max=36)  
# is there a problem?


model_bic = auto.arima(y,max.d=1,max.p=24,max.q=24,ic="bic")
model_aic = auto.arima(y,max.d=1,max.p=24,max.q=24,ic="aic")

checkresiduals(model_bic)
checkresiduals(model_aic)

model_general   <- model_aic
model_specific  <- arima(y, order = c(1,1,1), method = "ML")
lgen = logLik(model_general)  
lspec =logLik(model_specific)

dfs <- length(coef(model_general)) - length(coef(model_specific))
teststat<--2*(as.numeric(lspec-lgen))
pchisq(teststat,df=dfs,lower.tail=FALSE)  

forecast_arima <- forecast(model_aic, h = length(data_test.EPU))
forecast_hw  <- forecast(f_hw_opt, h = length(data_test.EPU))
accuracy_arima <- accuracy(forecast_aic, data_test.EPU)
accuracy_hw<- accuracy(forecast_hw, data_test.EPU)

print(accuracy_arima)
print(accuracy_hw)

