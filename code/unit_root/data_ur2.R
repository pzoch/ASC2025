
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
# load data

input_dir  = "C:/users/Piotr/Dropbox/WNE/ASC_2025/data/"
agg_data     = fread(paste0(input_dir,"/agg_data.csv"),header = TRUE)



forward_rate_list = c(paste0("TIPSF0",2:9),paste0("TIPSF",10:20))
hqm_TD_spread_list = paste0("hqm_TD_spread",c(2,5,10,20))


agg_data[, leverage := liquid_effective_leverage_agg]
agg_data = subset(agg_data,select = c("date","real_rate_1year",forward_rate_list,"CP_TB_spread",hqm_TD_spread_list,"leverage"))


agg_data      = na.omit(agg_data)
agg_data      = agg_data[order(date)]


# some plots
temp = melt(agg_data[,c("date","leverage")],id.vars="date")
ggplot()+
  geom_point(data = temp, aes(x=date, y= value, color=variable))+
  geom_line(data = temp, aes(x=date, y= value, color=variable))


y = agg_data$leverage
dy   = diff(y)

# plot acf and pacf
acf(coredata(y), lag.max = 12, plot = TRUE)
pacf(coredata(y), lag.max = 12, plot = TRUE)



# Test ADF
########################################
summary(ur.df(y, type="none", lags=4))
summary(ur.df(y, type="drift", lags=5))
summary(ur.df(y, type="trend", lags=5))
summary(ur.df(dy, type="drift", lags=4))


# BG test - ADF
# create "t" 
time = difftime(temp$date,min(temp$date), units = "weeks")

# need to specify regression manually - add as many lags of dy and as many deterministic regressors as in ur.df (line 53)
# line 63 runs the same regression as line 53 (we add time trend and an intercept (intercept is included automatically in R, if you want to get rid of it add "-1" in the formula))
reg_df = dynlm(zoo(dy) ~ L(zoo(y),1) + time + L(zoo(dy),1) +  L(zoo(dy),2) + L(zoo(dy),3) +  L(zoo(dy),4) + L(zoo(dy),5))
# line 65 verifies  if it is indeed the same regression - compare coefficients
summary(reg_df)
# now check for autocorrelation of residuals, first some plots
resid  = zoo(resid(reg_df), agg_data$date)
acf(coredata(resid), lag.max = 12, plot = TRUE)
pacf(coredata(resid), lag.max = 12, plot = TRUE)
# here Breusch-Godfrey test
bgtest(reg_df,order = 1)
bgtest(reg_df,order = 2)
bgtest(reg_df,order = 3)
bgtest(reg_df,order = 4)
bgtest(reg_df,order = 12)

# Test KPSS
########################################

summary(ur.kpss(y, type = "tau",lags = "short"))
summary(ur.kpss(dy, type = "mu",lags = "short"))


