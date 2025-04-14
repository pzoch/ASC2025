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
input_dir  = "D:/Dropbox (Personal)/WNE/ASC_2025/data/"
data_pl     = read.csv(file = paste0(input_dir, "GDP_POLAND.csv"))
colnames(data_pl) = c("date","Y")

data_pl$date = strptime(data_pl$date, "%m/%d/%Y" )
data_pl$Y    = zoo(data_pl$Y,order.by = data_pl$date)
data_pl$date = as.POSIXct(data_pl$date)
Y         =  data_pl$Y

# some plots
autoplot(Y)

# transform to logs

logY  <- ts(log(Y),start = c(1995, 1), frequency =4)
dlogY <- (diff(logY, lag = 1) * 100) # growth rate

autoplot(dlogY)

# forecasts from 2015 onwards
data_test <- window(dlogY, start = c(2015, 1))
data_train <- window(dlogY, end = c(2014, 4))

autoplot(data_test)
autoplot(data_train)

#### NAIVE FORECASTS
# naive forecast
GDP_forecast.naive  <- naive(data_train, h=length(data_test))
summary(GDP_forecast.naive)
autoplot(GDP_forecast.naive)

# plot without CIs
autoplot(data_train, series = "Train data") + 
  autolayer(GDP_forecast.naive$mean, series = "Point forecast") +
 autolayer(data_test, series = "Test data") + 
  xlab("Date") + 
  ylab("GDP growth")

#### MEAN FORECASTS
GDP_forecast.mean <-  meanf(data_train, h=length(data_test))
summary(GDP_forecast.mean)
autoplot(GDP_forecast.mean)

### Accuracy of these silly forecasts
accuracy(GDP_forecast.naive,data_test)
accuracy(GDP_forecast.mean,data_test)






####### MOVING AVERAGES

### SMOOTHING
# moving average - smoothing
# simple example
GDP_smooth.ma3 <- rollmean(dlogY, k = 3, fill = NA, align = "center")

# do plot many
# convert to dataframe
data_for_smoothing <- data.frame(dlogY = as.matrix(dlogY),date=as.Date(as.yearmon(time(dlogY))))

# mutate
smoothing <- data_for_smoothing %>% 
select(date, dlogY) %>%
  mutate(dlogY_ma02 = rollmean(dlogY, k = 3, fill = NA),
         dlogY_ma04 = rollmean(dlogY, k = 5, fill = NA),
         dlogY_ma06 = rollmean(dlogY, k = 7, fill = NA),
         dlogY_ma08 = rollmean(dlogY, k = 9, fill = NA),
         dlogY_ma10 = rollmean(dlogY, k = 11, fill = NA))

# reorganize...
smoothing = smoothing %>%
  gather(metric, value, dlogY:dlogY_ma10)

# ...and plot
  ggplot(data = smoothing,aes(date, value, color = metric)) +
  geom_line()


# can be done in one step!
  smoothing %>%
    gather(metric, value, dlogY:dlogY_ma10) %>%
    ggplot(aes(date, value, color = metric)) +
    geom_line() 
    
  
### forecast   
# moving average - forecast
# observe that k = 2 (use two quarters) and align = "right"
GDP_forecast.ma2 <- rollmean(dlogY, k = 2, fill = NA, align = "right")

# discuss what is going on here - this is not yet a forecast - this is just a trailing MA!


# convert to dataframe
data_for_forecasting <- data.frame(dlogY = as.matrix(dlogY),date=as.Date(as.yearmon(time(dlogY))))

# mutate - note there is "lag"!
forecast <- data_for_forecasting %>% 
  select(date, dlogY) %>%
  mutate(dlogY_MA_forecast = lag ( rollmean(dlogY, k = 4, align = "right",fill = NA), 1))

forecast %>%
  gather(metric, value, -date) %>%
  ggplot(aes(date, value, color = metric)) +
  geom_line()



# moving average - with pracma -- allows for more interesting weights

# compare with rollmean
forecast <- data_for_forecasting %>% 
  select(date, dlogY) %>%
  mutate(dlogY_MA_rollmean_forecast = lag ( rollmean(dlogY, k = 4, align = "right",fill = NA),1),
  dlogY_MA_pracma_forecast = lag ( movavg(dlogY, n = 4, type = 's'),1))


  
# compare weights
forecast <- data_for_forecasting %>% 
  select(date, dlogY) %>%
  mutate(dlogY_MA_s_forecast = lag ( movavg(dlogY, n = 4, type = 's'),1),
         dlogY_MA_t_forecast = lag ( movavg(dlogY, n = 4, type = 't'),1),
         dlogY_MA_w_forecast = lag ( movavg(dlogY, n = 4, type = 'w'),1),
         dlogY_MA_e_forecast = lag ( movavg(dlogY, n = 4, type = 'e'),1))

forecast %>%
  gather(metric, value, -date) %>%
  ggplot(aes(date, value, color = metric)) +
  geom_line()


### EXPONENTIAL SMOOTHING


GDP_forecast.exp001 <- ses((data_train), h=length(data_test),alpha=0.01)
GDP_forecast.exp03  <- ses((data_train), h=length(data_test),alpha=0.3)
GDP_forecast.exp05  <- ses((data_train), h=length(data_test),alpha=0.5)
GDP_forecast.exp07  <- ses((data_train), h=length(data_test),alpha=0.7)
GDP_forecast.exp099 <- ses((data_train), h=length(data_test),alpha=0.99)
GDP_forecast.expX   <- ses((data_train), h=length(data_test))


autoplot(GDP_forecast.exp001) +
  autolayer(fitted(GDP_forecast.exp001), series="Fitted") +
  ylab("GDP") + xlab("Year")

autoplot(GDP_forecast.exp03) +
  autolayer(fitted(GDP_forecast.exp03), series="Fitted") +
  ylab("GDP") + xlab("Year")

autoplot(GDP_forecast.exp05) +
  autolayer(fitted(GDP_forecast.exp05), series="Fitted") +
  ylab("GDP") + xlab("Year")

autoplot(GDP_forecast.exp07) +
  autolayer(fitted(GDP_forecast.exp07), series="Fitted") +
  ylab("GDP") + xlab("Year")

autoplot(GDP_forecast.exp099) +
  autolayer(fitted(GDP_forecast.exp099), series="Fitted") +
  ylab("GDP") + xlab("Year")

autoplot(GDP_forecast.expX) +
  autolayer(fitted(GDP_forecast.expX), series="Fitted") +
  ylab("GDP") + xlab("Year")

