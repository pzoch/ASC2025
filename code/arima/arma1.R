rm(list=ls())
set.seed(13) 

require(data.table)
require(ggplot2)
require(dyn)
require(zoo)
require(urca)
require(fpp2)
require(lmtest)
require(ggfortify)
require(lubridate)

# generate some MA
time = 1:1000
ma_vec_1 = c(1,1,1,1)
ma_vec_2 = c(-1,2,-3,4)
ma_vec_3 = c(2,1,0.5,0.2)
y1 <- data.table(time = 1:1000,y=arima.sim(model = list(ma = ma_vec_1), n = 1000,sd = 1))
y2 <- data.table(time = 1:1000,y=arima.sim(model = list(ma = ma_vec_2), n = 1000,sd = 1))
y3 <- data.table(time = 1:1000,y=arima.sim(model = list(ma = ma_vec_3), n = 1000,sd = 1))


ggplot(y1,aes(time,y))+geom_line(aes(color="MA VEC 1"))+
  geom_line(data=y2,aes(color="MA VEC 2"))+
  geom_line(data=y3,aes(color="MA VEC 3"))+
  labs(color="Legend")

par(mfrow=c(2,3))

acf(y1$y,lag.max=12)
acf(y2$y,lag.max=12)
acf(y3$y,lag.max=12)
pacf(y1$y,lag.max=12)
pacf(y2$y,lag.max=12)
pacf(y3$y,lag.max=12)

# generate some AR
time = 1:1000
ar_vec_1 = c(0.5,0,0,0)
ar_vec_2 = c(0.4,0.2,0.1,0.05)
ar_vec_3 = c(-0.3,0.2,-0.1,0.05)
z1 <- data.table(time = 1:1000,z=arima.sim(model = list(ar = ar_vec_1), n = 1000,sd = 1))
z2 <- data.table(time = 1:1000,z=arima.sim(model = list(ar = ar_vec_2), n = 1000,sd = 1))
z3 <- data.table(time = 1:1000,z=arima.sim(model = list(ar = ar_vec_3), n = 1000,sd = 1))


ggplot(z1,aes(time,z))+geom_line(aes(color="AR VEC 1"))+
  geom_line(data=z2,aes(color="AR VEC 2"))+
  geom_line(data=z3,aes(color="AR VEC 3"))+
  labs(color="Legend")

par(mfrow=c(2,3))

acf(z1$z,lag.max=12)
acf(z2$z,lag.max=12)
acf(z3$z,lag.max=12)
pacf(z1$z,lag.max=12)
pacf(z2$z,lag.max=12)
pacf(z3$z,lag.max=12)


# generate ARMA
# generate some MA
time = 1:100
ma_vec = c(0.3,3,0.1,0)
ar_vec = c(0.9,-0.3)
y_arma <-arima.sim(model = list(ma = ma_vec, ar = ar_vec),n=100,sd = 1)
  autoplot(y_arma)
  par(mfrow=c(1,2))
  acf(y_arma,lag.max=12) 
  pacf(y_arma,lag.max=12)     

model  <- arima(y_arma, order = c(2,0,3), method = "ML")

resids <- model$residuals
  par(mfrow=c(1,2))
  acf(resids,lag.max=12) 
  pacf(resids,lag.max=12)  
  coeftest(model)
  
# more specifications
  model32  <- arima(y_arma, order = c(3,0,2), method = "ML")
  model23  <- arima(y_arma, order = c(2,0,3), method = "ML")
  coeftest(model32)
  coeftest(model23)
# test autocorrelation
  
Box.test(resids,lag=12, type="Ljung-Box",fitdf = 2 + 3 )

checkresiduals(model)

# test normal
require(tseries)
JB <- jarque.bera.test(residuals(model))
JB



# roots 
require(fArma)
ar.coef <- model[["coef"]][1:2]
ma.coef <- model[["coef"]][3:5]
#ma.coef = NULL

roots <-armaRoots(ar.coef)
colnames(roots) <- c("re", "im", "abs")
options(digits=4)
roots

#MA representation
IRF <- ARMAtoMA(ar = ar.coef,ma = ma.coef,24)
plot.ts(IRF)

# information criteria

model_aic = auto.arima(y_arma,max.d=0,max.p=4,max.q=4,ic="aic")
model_bic = auto.arima(y_arma,max.d=0,max.p=4,max.q=4,ic="bic")
model32  <- arima(y_arma, order = c(3,0,2), method = "ML")  
model23  <- arima(y_arma, order = c(2,0,3), method = "ML")  

ar_aic.coef <- model_aic[["coef"]][1:1]
ma_aic.coef <- model_aic[["coef"]][2:3]

ar_bic.coef <- c(0)
ma_bic.coef <- model_bic[["coef"]][1:2]

ar_32.coef <- model32[["coef"]][1:3]
ma_32.coef <- model32[["coef"]][4:5]

ar_23.coef <- model23[["coef"]][1:2]
ma_23.coef <- model23[["coef"]][3:5]



resids  <- model12$residuals
coeftest(model12)
Box.test(resids,lag=12, type="Ljung-Box",fitdf = 1 + 2)
checkresiduals(model12)

#MA representation
IRF_aic <- ARMAtoMA(ar = ar_aic.coef,ma = ma_aic.coef,24)
IRF_bic <- ARMAtoMA(ar = ar_bic.coef,ma = ma_bic.coef,24)
IRF_23 <- ARMAtoMA(ar = ar_23.coef,ma = ma_23.coef,24)
IRF_true <- ARMAtoMA(ar = ar_vec,ma = ma_vec,24)


par(mfrow=c(1,4))
plot.ts(IRF_aic)
plot.ts(IRF_bic)
plot.ts(IRF_23)
plot.ts(IRF_true)

IRF_ic <- ARMAtoMA(ar = 0,ma = ma.coef,24)
plot.ts(IRF_ic)

#MA representation - true
IRF_true <- ARMAtoMA(ar = ar_vec,ma = ma_vec,24)
plot.ts(IRF_true)


IRF_true <- IRF_true / IRF_true[1]

