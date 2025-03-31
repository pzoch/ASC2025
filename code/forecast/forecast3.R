rm(list=ls())

require(zoo)
require(eurostat)
require(ggfortify)
require(forecast)
require(smooth)
require(pracma)
require(TSstudio)

# load data
input_dir  = "C:/users/Piotr/Dropbox/WNE/ASC_2023/data/"
data_us     = read.csv(file = paste0(input_dir, "CPI_US_NSA.csv"))
colnames(data_us) = c("date","Pi")

data_us$date = strptime(data_us$date, "%m/%d/%Y" )
data_us$Pi    = zoo(data_us$Pi,order.by = data_us$date)
data_us$date = as.POSIXct(data_us$date)

temp_Pi      =  ts(data_us$Pi,start = c(1960, 1,1), frequency =12)

# some plots
autoplot(temp_Pi)

Pi <- window(temp_Pi, start = c(2000,1))
autoplot(Pi)

# residual analysis
f_naive  <- naive(Pi)
f_snaive <- snaive(Pi)
f_meanf  <- meanf(Pi)
f_ses    <- ses(Pi,"alpha"=0.25)

autoplot(resid(f_naive)) + xlab("Time") + ylab("") +
  ggtitle("Residuals from naive method")

autoplot(resid(f_ses)) + xlab("Time") + ylab("") +
  ggtitle("Residuals from simple exponential smoothing")


# do everything simultaneously
checkresiduals(f_naive)
checkresiduals(f_snaive)
checkresiduals(f_meanf)
checkresiduals(f_ses)


# use holt-winters from now on
data_test <- window(Pi, start = c(2018, 1),end = c(2018, 12))
data_train   <- window(Pi, end = c(2017, 12))
f_hw1    <- hw((data_train), h=12,"alpha"=0.45,"beta"=0.1,"gamma"=0.50,seasonal="additive")
f_hw2    <- hw((data_train), h=12,"alpha"=0.25,"beta"=0.1,"gamma"=0.50,seasonal="additive")
f_hw3    <- hw((data_train), h=12,"alpha"=0.15,"beta"=0.1,"gamma"=0.50,seasonal="additive")
f_naive    <- naive((data_train), h=12,"alpha"=0.15,"beta"=0.1,"gamma"=0.50,seasonal="additive")
f_snaive    <- snaive((data_train), h=12,"alpha"=0.15,"beta"=0.1,"gamma"=0.50,seasonal="additive")

autoplot(f_hw1) +
  autolayer(fitted(f_hw1), series="Fitted - Holt-Winters") +
  autolayer(data_test, series="Data - test") +
  ylab("Inflation") + xlab("Year")

checkresiduals(f_hw1)
f_hw_opt    <- hw((data_train), h=12,seasonal="multiplicative")
checkresiduals(f_hw_opt)


accuracy(f_hw1, data_test)
accuracy(f_hw2, data_test)
accuracy(f_hw3, data_test)
accuracy(f_hw_opt, data_test)
accuracy(f_naive, data_test)
accuracy(f_snaive, data_test)

# time series cross-test

errors_tsCV <- tsCV(data_train, forecastfunction=hw, h=1,"alpha"=0.45,"beta"=0.1,"gamma"=0.50,seasonal="additive")


errors_base <- residuals(hw(data_train, h=1,"alpha"=0.45,"beta"=0.1,"gamma"=0.50,seasonal="additive"))

rmse_tsCV <- sqrt(mean(errors_tsCV^2, na.rm=TRUE))
rmse_base <- sqrt(mean(errors_base^2, na.rm=TRUE))

# which model is the best?
errors_tsCV1 <- tsCV(data_train, forecastfunction=hw, h=12,"alpha"=0.45,"beta"=0.1,"gamma"=0.50,seasonal="additive")
errors_tsCV2 <- tsCV(data_train, forecastfunction=hw, h=12,"alpha"=0.35,"beta"=0.1,"gamma"=0.10,seasonal="additive")
errors_tsCV3 <- tsCV(data_train, forecastfunction=hw, h=12,"alpha"=0.25,"beta"=0.05,"gamma"=0.10,seasonal="additive")
errors_tsCV4 <- tsCV(data_train, forecastfunction=snaive, h=12)

rmse_tsCV1 <- sqrt(mean(errors_tsCV1[,12]^2, na.rm=TRUE))
rmse_tsCV2 <- sqrt(mean(errors_tsCV2[,12]^2, na.rm=TRUE))
rmse_tsCV3 <- sqrt(mean(errors_tsCV3[,12]^2, na.rm=TRUE))
rmse_tsCV4 <- sqrt(mean(errors_tsCV4[,12]^2, na.rm=TRUE))

# which model performed best?

f_1  <- hw((data_train), h=12,"alpha"=0.45,"beta"=0.1,"gamma"=0.50,seasonal="additive")
f_2  <- hw((data_train), h=12,"alpha"=0.35,"beta"=0.1,"gamma"=0.10,seasonal="additive")
f_3  <- hw((data_train), h=12,"alpha"=0.25,"beta"=0.05,"gamma"=0.10,seasonal="additive")
f_4  <- snaive((data_train), h=12)

accuracy(f_1, data_test)
accuracy(f_2, data_test)
accuracy(f_3, data_test)
accuracy(f_4, data_test)

# bootstrap vs no bootstrap
f_snaive_nb   <-  snaive(data_train,bootstrap = FALSE)
f_snaive_b    <-  snaive(data_train,bootstrap = TRUE)


autoplot(f_snaive_nb) 
autoplot(f_snaive_b)


