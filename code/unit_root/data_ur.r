

# testing stationarity

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
# load data -- better to write it as a relative path...
input_dir  = "C:/users/Piotr/Dropbox/WNE/ASC_2025/data/"

data_pl     = read.csv(file = paste0(input_dir, "GDP_POLAND.csv"))
colnames(data_pl) = c("date","Y")

data_pl$date = strptime(data_pl$date, "%m/%d/%Y" )
data_pl$Y    = zoo(data_pl$Y,order.by = data_pl$date)
data_pl$date = as.POSIXct(data_pl$date)

Y         =  data_pl$Y
Y_short   =  window(data_pl$Y, end = parse_date_time("01/01/2015", orders = "%m-%d-%Y"))

# some plots
autoplot(Y)
autoplot(Y_short)


# transform to logs
y          = log(Y)
y_short    = log(Y_short)
autoplot(y)
autoplot(y_short)

# first log differences --- $ growth
dy   = diff(y)
dy_short   = diff(y_short)

autoplot(dy)
autoplot(dy_short)


# detrending
time = difftime(data_pl$date,min(data_pl$date), units = "weeks")
time_2 = as.numeric(time)^2

reg_detr = lm(y ~ time)
y_detr   = zoo(resid(reg_detr), data_pl$date)

autoplot(y_detr)

reg_detr_2 = lm(y ~ time +  time_2)
y_detr_2   = zoo(resid(reg_detr_2), data_pl$date)

autoplot(y_detr_2)

# plot acf and pacf

acf(coredata(y), lag.max = 12, plot = TRUE)
acf(coredata(y_short), lag.max = 12, plot = TRUE)

pacf(coredata(y), lag.max = 12, plot = TRUE)
pacf(coredata(y_short), lag.max = 12, plot = TRUE)

acf(coredata(dy), lag.max = 12, plot = TRUE)
acf(coredata(dy_short), lag.max = 12, plot = TRUE)

pacf(coredata(dy), lag.max = 12, plot = TRUE)
pacf(coredata(dy_short), lag.max = 12, plot = TRUE)

# Test ADF
########################################
summary(ur.df(y, type="none", lags=4))
summary(ur.df(y, type="drift", lags=4))
summary(ur.df(y, type="trend", lags=4))

summary(ur.df(dy, type="drift", lags=4))

summary(ur.df(y_short, type="none", lags=4))
summary(ur.df(y_short, type="drift", lags=4))
summary(ur.df(y_short, type="trend", lags=4))
summary(ur.df(dy_short, type="drift", lags=4))

# BG test - ADF
reg_df = dynlm(dy ~ lag(y,-1) - 1) 
summary(reg_df)
resid  = zoo(resid(reg_df), data_pl$date)

acf(coredata(resid), lag.max = 12, plot = TRUE)
pacf(coredata(resid), lag.max = 12, plot = TRUE)


bgtest(reg_df,order = 1)
bgtest(reg_df,order = 2)
bgtest(reg_df,order = 3)
bgtest(reg_df,order = 4)


# Test KPSS
########################################

summary(ur.kpss(y, type = "tau",lags = "short"))
summary(ur.kpss(dy, type = "mu",lags = "short"))

summary(ur.kpss(y_short, type = "tau",lags = "short"))
summary(ur.kpss(dy_short, type = "mu",lags = "short"))

