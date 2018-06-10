# Bond US

library("forecast")
library("fpp")
library("fpp2")
library("GGally")
library("xts");
library("lubridate")
library("data.table")


print(getwd())

# Read CSV into R
data = read.csv(file="test_bond.csv")

# asDate
data$date <- as.Date(data$date, format = "%Y-%m-%d")

data = data[with(data, order(data$date, decreasing = FALSE)),]

head(data)
tail(data)

nByYear = table(year(data$date))
nByYear

nByYear = data.frame(table(year(data$date)))
columns <- c("year", "freq")
colnames(nByYear) <- columns
nByYear$year = as.numeric(as.character(nByYear$year))
nByYear


sapply(nByYear, class)

nByYear <- subset(nByYear, year > 1962 & year < 2018)

nByYear

mean_Yfreq = mean(nByYear$freq, na.rm = TRUE)
mean_freq = round(mean_Yfreq)
mean_freq

data <- subset(data, date > as.Date("1962-01-01"))

head(data)
tail(data)

# get last date
lastdate = max(data$date)

# get startdate
startdate = min(data$date)

print(lastdate)
print(startdate)


# make index
rownames(data) <- data$date

data$date <- NULL

# check
head(data)
tail(data)


daily_ts = ts(data, frequency = 253)
daily_msts <- msts(data, seasonal.periods=c(5, 253))

head(daily_ts)
tail(daily_ts)

plot(daily_ts)

# optional
monthly_ts = ts(data, start=c(year(startdate),month(startdate)),frequency=12)
quart_ts = ts(data, start=c(year(startdate),quarter(startdate)),frequency=4)

ggseasonplot(monthly_ts)
ggmonthplot(monthly_ts)

# step a head
tau=10;

n=length(daily_ts);
print(n)


loss.all=list();

loss.functions = function(x.hat, x)
{
  return(c(mean((x-x.hat)^2), mean(abs(x-x.hat)), mean(abs( (x-x.hat)/x )) ));
}


my.forecast = function (which.model)
{
  my.mod = NULL;
    # MAFS
  if (which.model==1) {ee.for = tail(rollapply(daily_ts, width = 3, by = 1, FUN = mean, align = "right"), tau);}
  if (which.model==2) {ee.for = tail(rollapply(daily_ts, width = 5, by = 1, FUN = mean, align = "right"), tau);}
  if (which.model==3) {ee.for = tail(rollapply(daily_ts, width = 7, by = 1, FUN = mean, align = "right"), tau);}

  if (which.model==4) {
    # EWMA
    my.mod = HoltWinters(daily_ts, beta=FALSE, gamma=FALSE);
    ee.for = tail(my.mod$fitted[,1], tau);
  }
  if (which.model==5) {
    # Holt
    my.mod = HoltWinters(daily_ts, gamma=FALSE);
    ee.for = tail(my.mod$fitted[,1], tau);
  }
  if (which.model==6) {
    # Holt-Winters
    my.mod = HoltWinters(daily_ts);
    ee.for = tail(my.mod$fitted[,1], tau);
  }
  if (which.model==7) {
    # TBATS
    my.mod = tbats(daily_msts);
    ee.for = tail(my.mod$fitted[,1], tau);
  }
  return(list(my.mod, ee.for))
}

# gathering all forecasts
all.forec = matrix(0, ncol=7, nrow=tau);
for(i in 1:7)
{
  all.forec[,i] = my.forecast(i)[[2]];
}

#  compute the corresponding losses.
losses = apply(all.forec, 2, function(x) loss.functions(x, tail(daily_ts,tau)))
losses

# plot
model.names = c("maf 3", "maf 5", "maf 7", "emwa", "holt", "holtwinters", "tbats")
plot.names = c("true", model.names)
all = cbind(tail(daily_ts,tau),all.forec);
names(all) <- plot.names;
plot.ts(all, plot.type="single", col=c(1:6,1:5), lty=c(rep(1,6),rep(2,6)), ylab="forecasts")
legend(x = "bottomleft", legend=model.names, ncol=3, bty="n", col=c(1:6,1:5), lty=c(rep(1,6),rep(2,6)))

df_losses = data.frame(losses)
colnames(df_losses) <- model.names
rownames(df_losses) <- c("mse", "mae", "mape")
df_losses
write.csv(df_losses, file = "R_losses.csv")

