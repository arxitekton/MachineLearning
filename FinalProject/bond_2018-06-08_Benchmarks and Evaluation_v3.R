# BOND
library(Quandl)
library("forecast")
library("fpp")
library("fpp2")
library("GGally")
library("xts");
library("lubridate")
library("data.table")
library("ggplot2")


print(getwd())

# Read CSV into R
indices = read.csv(file="bond_indices.csv", as.is = 1)
num_of_indices = NROW(indices)

indices
indices$index[1]
#################################################################################

fh <- 10 #The forecasting horizon examined
frq <- 365 #The frequency of the data

start_date <- "1962-01-01"
end_date <- "2018-06-08"

Quandl.api_key("QWERTY")

test = Quandl(indices$index[1], start_date=start_date, end_date=end_date, collapse='daily')
typeof(test)
tail(test)
colnames(test)[colnames(test) == 'Date'] <- "dates"
colnames(test)[colnames(test) == 'Value'] <- "values"

alldates = seq(min(test$dates), max(test$dates), 1)

# Filter out timestamps that are already present in your `data.frame`:
# Construct a `data.frame` to append with missing values:
dates0 = alldates[!(alldates %in% test$dates)]
data0 = data.frame(dates = dates0, values = NA_real_)

# Append this `data.frame` and resort in time:
test = rbind(test, data0)
test = test[order(test$dates),]

# forward fill the values 
current = NA_real_
test$values = sapply(test$values, function(x) { 
  current <<- ifelse(is.na(x), current, x); current })



tail(test)
#test$date <- as.Date(test$date, format = "%Y-%m-%d")
#test = test[with(test, order(test$date, decreasing = FALSE)),]
tail(test,20)
head(test$value)
typeof(test$value)

tail(test,10)

data_train = data_test <- NULL #Train and test sample
for (i in 1:num_of_indices){
  quandl_data <- Quandl(indices$index[i], start_date=start_date, end_date=end_date, collapse='daily')
  colnames(quandl_data)[colnames(quandl_data) == 'Date'] <- "dates"
  colnames(quandl_data)[colnames(quandl_data) == 'Value'] <- "values"
  quandl_data$dates <- as.Date(quandl_data$dates, format = "%Y-%m-%d")
  #data <- data[with(data, order(data$dates, decreasing = FALSE)),]
  
  alldates = seq(min(quandl_data$dates), max(quandl_data$dates), 1)
  
  # Filter out timestamps that are already present in your `data.frame`:
  # Construct a `data.frame` to append with missing values:
  dates0 = alldates[!(alldates %in% quandl_data$dates)]
  data0 = data.frame(dates = dates0, values = NA_real_)
  
  # Append this `data.frame` and resort in time:
  data = rbind(quandl_data, data0)
  data = data[order(data$dates),]
  
  # forward fill the values 
  current = NA_real_
  data$values = sapply(data$values, function(x) { 
    current <<- ifelse(is.na(x), current, x); current })
  tail(data)
  
  
  data_all <- data$values
  data_train[length(data_train)+1] <- list(ts(head(data_all,length(data_all)-fh),frequency = frq))
  data_test[length(data_test)+1] <- list(tail(data_all,fh))
}
#################################################################################

smape_cal <- function(outsample, forecasts){
  #Used to estimate sMAPE
  outsample <- as.numeric(outsample) ; forecasts<-as.numeric(forecasts)
  smape <- (abs(outsample-forecasts)*200)/(abs(outsample)+abs(forecasts))
  return(smape)
}

mase_cal <- function(insample, outsample, forecasts){
  #Used to estimate MASE
  frq <- frequency(insample)
  forecastsNaiveSD <- rep(NA,frq)
  for (j in (frq+1):length(insample)){
    forecastsNaiveSD <- c(forecastsNaiveSD, insample[j-frq])
  }
  masep<-mean(abs(insample-forecastsNaiveSD),na.rm = TRUE)
  
  outsample <- as.numeric(outsample) ; forecasts <- as.numeric(forecasts)
  mase <- (abs(outsample-forecasts))/masep
  return(mase)
}

mape_cal <- function(x, x.hat){
  #Used to estimate MAPE
  return(mean(abs( (x-x.hat)/x )))
}

naive_seasonal <- function(input, fh){
  #Used to estimate Seasonal Naive
  frcy <- frequency(input)
  frcst <- naive(input, h=fh)$mean 
  if (frcy>1){ 
    frcst <- head(rep(as.numeric(tail(input,frcy)), fh), fh) + frcst - frcst
  }
  return(frcst)
}

Theta.classic <- function(input, fh){
  #Used to estimate Theta classic
  
  #Set parameters
  wses <- wlrl<-0.5 ; theta <- 2
  #Estimate theta line (0)
  observations <- length(input)
  xt <- c(1:observations)
  xf <- c((observations+1):(observations+fh))
  train <- data.frame(input=input, xt=xt)
  test <- data.frame(xt = xf)
  
  estimate <- lm(input ~ poly(xt, 1, raw=TRUE))
  thetaline0In <- as.numeric(predict(estimate))
  thetaline0Out <- as.numeric(predict(estimate,test))
  
  #Estimate theta line (2)
  thetalineT <- theta*input+(1-theta)*thetaline0In
  sesmodel <- ses(thetalineT, h=fh)
  thetaline2In <- sesmodel$fitted
  thetaline2Out <- sesmodel$mean
  
  #Theta forecasts
  forecastsIn <- (thetaline2In*wses)+(thetaline0In*wlrl)
  forecastsOut <- (thetaline2Out*wses)+(thetaline0Out*wlrl)
  
  #Zero forecasts become positive
  for (i in 1:length(forecastsOut)){
    if (forecastsOut[i]<0){ forecastsOut[i]<-0 }
  }
  
  output=list(fitted = forecastsIn, mean = forecastsOut,
              fitted0 = thetaline0In, mean0 = thetaline0Out,
              fitted2 = thetaline2In, mean2 = thetaline2Out)
  
  return(output)
}

SeasonalityTest <- function(input, ppy){
  #Used to determine whether a time series is seasonal
  tcrit <- 1.645
  if (length(input)<3*ppy){
    test_seasonal <- FALSE
  }else{
    xacf <- acf(input, plot = FALSE)$acf[-1, 1, 1]
    clim <- tcrit/sqrt(length(input)) * sqrt(cumsum(c(1, 2 * xacf^2)))
    test_seasonal <- ( abs(xacf[ppy]) > clim[ppy] )
    
    if (is.na(test_seasonal)==TRUE){ test_seasonal <- FALSE }
  }
  
  return(test_seasonal)
}

f_maf_3 <- function(input, fh){
  predict = tail(rollapply(input, width = 3, by = 1, FUN = mean, align = "right"), fh)
  return(predict)
}

f_maf_5 <- function(input, fh){
  predict = tail(rollapply(input, width = 5, by = 1, FUN = mean, align = "right"), fh)
  return(predict)
}

f_maf_7 <- function(input, fh){
  predict = tail(rollapply(input, width = 7, by = 1, FUN = mean, align = "right"), fh)
  return(predict)
}

f_ets <- function(input, fh){
  fit = ets(input);
  predict = forecast(fit, fh)
  return(predict)
}

f_EWMA <- function(input, fh){
  fit <- HoltWinters(input, beta=FALSE, gamma=FALSE);
  predict <- forecast(fit, fh)
  return(predict)
}

f_Holt <- function(input, fh){
  fit <- HoltWinters(input, gamma=FALSE);
  predict <- forecast(fit, fh)
  return(predict)
}

f_HoltWinters <- function(input, fh){
  fit <- HoltWinters(input);
  predict <- forecast(fit, fh)
  return(predict)
}

f_tbats <- function(input, fh){
  fit <- tbats(input);
  predict <- forecast(fit, fh)
  return(predict)
}

f_autoArima <- function(input, fh){
  fit <- auto.arima(input);
  predict <- forecast(fit, fh)
  return(predict)
}

Benchmarks <- function(input, fh){
  #Used to estimate the statistical benchmarks of the M4 competition
  
  #Estimate seasonaly adjusted time series
  ppy <- frequency(input) ; ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input,ppy) }
  if (ST==T){
    Dec <- decompose(input,type="multiplicative")
    des_input <- input/Dec$seasonal
    SIout <- head(rep(Dec$seasonal[(length(Dec$seasonal)-ppy+1):length(Dec$seasonal)], fh), fh)
  }else{
    des_input <- input ; SIout <- rep(1, fh)
  }
  
  f1 <- naive(input, h=fh)$mean # Naive
  f1_0 <- rwf(input, h=fh)$mean # Alternative Naive
  f1_1 <- snaive(input, h=fh)$mean # Seasonal naïve
  f2 <- naive_seasonal(input, fh=fh) # Seasonal Naive2
  f3 <- naive(des_input, h=fh)$mean*SIout # Naive2
  f3_1 <- rwf(input, h=fh, drift=TRUE)$mean # Drift
  f4_0 <- meanf(input, h=fh)$mean # Average
  f4 <- f_maf_3(input, fh=fh) # maf3
  f5 <- f_maf_5(input, fh=fh) # maf5
  f6 <- f_maf_7(input, fh=fh) # maf7
  f7 <- ses(des_input, h=fh)$mean*SIout # Ses
  f8 <- holt(des_input, h=fh, damped=F)$mean*SIout # Holt
  f9 <- holt(des_input, h=fh, damped=T)$mean*SIout # Damped
  f10 <- Theta.classic(input=des_input, fh=fh)$mean*SIout # Theta
  f11 <- f_EWMA(input, fh=fh)$mean # EWMA
  f11_1 <- f_Holt(input, fh=fh)$mean # Holt2
  f12 <- f_HoltWinters(input, fh=fh)$mean # HoltWinters
  f13 <- f_tbats(input, fh=fh)$mean # tbats
  f14 <- f_autoArima(input, fh=fh)$mean # auto.arima
  f15 <- f_ets(input, fh=fh)$mean # ets
  
  return(list(f1,f1_0,f1_1,f2,f3,f3_1,f4_0,f4,f5,f6,f7,f8,f9,f10,f11,f11_1,f12,f13,f14,f15))
}

Names_benchmarks <- c("Naive","Alternative Naive","Seasonal naïve", "sNaive2", "Naive2", "Drift","Average",'maf3','maf5','maf7', "SES", "Holt", "Damped", "Theta", "EWMA", "Holt2", 'HoltWinters','tbats','auto.arima',"ets")
Total_smape=Total_mase=Total_mape <- array(NA,dim = c(length(Names_benchmarks), fh, length(data_train)))
#Methods, Horizon, time-series
for (i in 1:length(data_train)){
  
  insample <- data_train[[i]]
  outsample <- data_test[[i]]
  forecasts <- Benchmarks(input=insample, fh=fh)
  
  #sMAPE
  for (j in 1:length(Names_benchmarks)){
    Total_smape[j,,i] <- smape_cal(outsample, forecasts[[j]]) #j the # of the benchmark
  }
  #MASE
  for (j in 1:length(Names_benchmarks)){
    Total_mase[j,,i] <- mase_cal(insample, outsample, forecasts[[j]]) #j the # of the benchmark
  }
  #MAPE
  for (j in 1:length(Names_benchmarks)){
    Total_mape[j,,i] <- mape_cal(outsample, forecasts[[j]]) #j the # of the benchmark
  }
  
}

print("########### sMAPE ###############")
for (i in 1:length(Names_benchmarks)){
  print(paste(Names_benchmarks[i], round(mean(Total_smape[i,,]), 5)))
}
print("########### MASE ################")
for (i in 1:length(Names_benchmarks)){
  print(paste(Names_benchmarks[i], round(mean(Total_mase[i,,]), 5)))
}
print("########### MAPE ################")
for (i in 1:length(Names_benchmarks)){
  print(paste(Names_benchmarks[i], round(mean(Total_mape[i,,]), 5)))
}

print("########### OWA ################")
for (i in 1:length(Names_benchmarks)){
  print(paste(Names_benchmarks[i],
              round(((mean(Total_mase[i,,])/mean(Total_mase[3,,]))+(mean(Total_smape[i,,])/mean(Total_smape[3,,])))/2, 3)))
}


