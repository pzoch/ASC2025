rm(list=ls())

require(zoo)
require(ggfortify)
require(forecast)
require(smooth)
require(pracma)
require(TSstudio)
require(tidyverse)
require(dplyr)

# load data
input_dir  = "C:/users/Piotr/Dropbox/WNE/ASC_2025/data/"
data_pl_sa     = read.csv(file = paste0(input_dir, "GDP_POLAND.csv"))
data_pl_nsa     = read.csv(file = paste0(input_dir, "GDP_POLAND_NSA.csv"))
colnames(data_pl_sa) = c("date","Y")
colnames(data_pl_nsa) = c("date","Y")

data_pl_sa$date = strptime(data_pl_sa$date, "%m/%d/%Y" )
data_pl_sa$Y    = zoo(data_pl_sa$Y,order.by = data_pl_sa$date)
data_pl_sa$date = as.POSIXct(data_pl_sa$date)

data_pl_nsa$date = strptime(data_pl_nsa$date, "%m/%d/%Y" )
data_pl_nsa$Y    = zoo(data_pl_nsa$Y,order.by = data_pl_nsa$date)
data_pl_nsa$date = as.POSIXct(data_pl_nsa$date)

Y_SA         =  data_pl_sa$Y
Y_NS         =  data_pl_nsa$Y
# some plots
autoplot(Y_SA)
autoplot(Y_NS)
# transform to logs

logY_SA  <- ts(log(Y_SA),start = c(1995, 1), frequency =4)
dlogY_SA <- (diff(logY_SA, lag = 1) * 100) # growth 
logY_NS  <- ts(log(Y_NS),start = c(1995, 1), frequency =4)
dlogY_NS <- (diff(logY_NS, lag = 1) * 100) # growth rate

p1 <- autoplot(logY_SA)
p2 <- autoplot(logY_NS)

gridExtra::grid.arrange(p1, p2, nrow = 1)

# SES - do it on differences first
data_test.dSA    <- window(dlogY_SA, start = c(2015, 1))
data_train.dSA   <- window(dlogY_SA, end = c(2014, 4))

data_test.dNS <- window(dlogY_NS, start = c(2015, 1))
data_train.dNS   <- window(dlogY_NS, end = c(2014, 4))

GDP_forecast.dSA <- ses((data_train.dSA), h=length(data_test.dSA),alpha=0.3)
GDP_forecast.dNS <- ses((data_train.dNS), h=length(data_test.dNS),alpha=0.3)

autoplot(GDP_forecast.dSA) +
  ylab("dlog GDP") + xlab("Year")

autoplot(GDP_forecast.dNS) +
  autolayer(fitted(GDP_forecast.dNS), series="Fitted - NS") +
  autolayer(data_test.dNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")


# SES - do it on levels now
data_test.lSA    <- window(logY_SA, start = c(2015, 1))
data_train.lSA   <- window(logY_SA, end = c(2014, 4))

data_test.lNS <- window(logY_NS, start = c(2015, 1))
data_train.lNS   <- window(logY_NS, end = c(2014, 4))

GDP_forecast.lSA <- ses((data_train.lSA), h=length(data_test.lSA),alpha=0.3)
GDP_forecast.lNS <- ses((data_train.lNS), h=length(data_test.lNS),alpha=0.3)

autoplot(GDP_forecast.lSA) +
  autolayer(fitted(GDP_forecast.lSA), series="Fitted - SA") +
  autolayer(data_test.lSA, series="Data - test") +
  ylab("log GDP") + xlab("Year")

autoplot(GDP_forecast.lNS) +
  autolayer(fitted(GDP_forecast.lNS), series="Fitted - NS") +
  autolayer(data_test.lNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

# Holt - do it on levels
data_test.hlSA <- window(logY_SA, start = c(2015, 1))
data_train.hlSA   <- window(logY_SA, end = c(2014, 4))

data_test.hlNS <- window(logY_NS, start = c(2015, 1))
data_train.hlNS   <- window(logY_NS, end = c(2014, 4))

GDP_forecast.hlSA <- holt((data_train.hlSA), h=length(data_test.hlSA),alpha=0.3,beta=0.3)
GDP_forecast.hlNS <- holt((data_train.hlNS), h=length(data_test.hlNS),alpha=0.3,beta=0.3)

autoplot(GDP_forecast.hlSA) +
  autolayer(fitted(GDP_forecast.hlSA), series="Fitted - SA") +
  autolayer(data_test.hlSA, series="Data - test") +
  ylab("log GDP") + xlab("Year")

autoplot(GDP_forecast.hlNS) +
  autolayer(fitted(GDP_forecast.hlNS), series="Fitted - NS") +
  autolayer(data_test.hlNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

## play a bit with beta to see what happens!

# Holt-Winters - do it on levels
data_test.hlSA <- window(logY_SA, start = c(2015, 1))
data_train.hlSA   <- window(logY_SA, end = c(2014, 4))

data_test.hlNS <- window(logY_NS, start = c(2015, 1))
data_train.hlNS   <- window(logY_NS, end = c(2014, 4))

GDP_forecast.hlSA <- hw((data_train.hlSA), h=length(data_test.hlSA),alpha=0.3,beta=0.3,gamma=0.1,seasonal="additive")
GDP_forecast.hlNS <- hw((data_train.hlNS), h=length(data_test.hlNS),alpha=0.3,beta=0.3,gamma=0.1,seasonal="additive")


  ##nb: use HW for decomposition
  plot(GDP_forecast.hlNS$model$states)
  ###
  
  
autoplot(GDP_forecast.hlSA) +
  autolayer(fitted(GDP_forecast.hlSA), series="Fitted - SA") +
  autolayer(data_test.hlSA, series="Data - test") +
  ylab("log GDP") + xlab("Year")

autoplot(GDP_forecast.hlNS) +
  autolayer(fitted(GDP_forecast.hlNS), series="Fitted - NS") +
  autolayer(data_test.hlNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

## play a bit with parameters to see what happens!

GDP_forecast.a01b01g05 <- hw((data_train.hlNS), h=length(data_test.hlNS),alpha=0.1,beta=0.1,gamma=0.5,seasonal="additive")
GDP_forecast.a01b01g01 <- hw((data_train.hlNS), h=length(data_test.hlNS),alpha=0.1,beta=0.1,gamma=0.1,seasonal="additive")
GDP_forecast.a02b02g02 <- hw((data_train.hlNS), h=length(data_test.hlNS),alpha=0.2,beta=0.2,gamma=0.5,seasonal="additive")
GDP_forecast.opt <- hw((data_train.hlNS), h=length(data_test.hlNS),seasonal="additive")


autoplot(GDP_forecast.a01b01g05) +
  autolayer(fitted(GDP_forecast.a01b01g05), series="Fitted - NS") +
  autolayer(data_test.hlNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

autoplot(GDP_forecast.a01b01g01) +
  autolayer(fitted(GDP_forecast.a01b01g01), series="Fitted - NS") +
  autolayer(data_test.hlNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

autoplot(GDP_forecast.a02b02g02) +
  autolayer(fitted(GDP_forecast.a02b02g02), series="Fitted - NS") +
  autolayer(data_test.hlNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

autoplot(GDP_forecast.opt) +
  autolayer(fitted(GDP_forecast.opt), series="Fitted - NS") +
  autolayer(data_test.hlNS, series="Data - test") +
  ylab("log GDP") + xlab("Year")

