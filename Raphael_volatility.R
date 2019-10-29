rm(list = ls())

library(readr)
library(ggplot2)
library(scales)
library(dplyr)
library(lubridate)
library(tidyverse)
library(zoo)
library(xts)
library(forecast)
library(tseries)
library(urca)
library(sweep)

# ***** Import Data *****
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)

FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)

last_day <- 20150430
last_month <- 201504

last_entry_day <- which(FF_daily[,1] == last_day)
last_entry_month <- which(FF_monthly[,1] == last_month)
FF_daily <- head(FF_daily, last_entry_day)
FF_monthly <- head(FF_monthly, last_entry_month)

# ***** Manipulate Data *****
FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")

# ***** Calculate Monthly Variances *****
trading_days <- 22
trading_months <- 12

monthly_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(Mkt) * trading_days)

monthly_vars <- monthly_vars %>% mutate(volatility = sqrt(variance))
plot(monthly_vars$volatility * sqrt(trading_months), type = "l")

# ******** Scenario 1: ARIMA model**

variance_ts <- xts(monthly_vars$variance, order.by = FF_monthly$Date)

ARMA_var <- c(1:last_entry_month)
for (i in 1:last_entry_month) {
  if (i <= 12) ARMA_var[i] <- variance_ts[i]
  else {
    model <- auto.arima(variance_ts[1:i], d = 0)
    ARMA_var[i] <- as.numeric(forecast(model, h = 1)$mean)
  }
}

monthly_vars$ARMA_var <- ARMA_var

# ***** Scenario 2: EWMA ****

n <- 0.94
EWMA_var <- c(1:last_entry_month)
for (i in 1:last_entry_month) {
  if (i == 1) EWMA_var[i] <- monthly_vars$variance[i]
  else {
    EWMA_var[i] <- n * EWMA_var[i - 1] + (1 - n) * monthly_vars$variance[i]
  }
}

monthly_vars$EWMA_var <- EWMA_var

# ***** Calculate c *****
c0 = sqrt(var(FF_monthly$Mkt[-1])/
           var(1/monthly_vars$variance[-last_entry_month]*FF_monthly$Mkt[-1]))

c1 = sqrt(var(FF_monthly$Mkt[-1])/
           var(1/monthly_vars$ARMA_var[-last_entry_month]*FF_monthly$Mkt[-1]))

c2 = sqrt(var(FF_monthly$Mkt[-1])/
            var(1/monthly_vars$EWMA_var[-last_entry_month]*FF_monthly$Mkt[-1]))

# ***** Calculate weights and volatility managed returns *****
weights0 <- c(1:(last_entry_month-1))
var_returns0 <- c(1:(last_entry_month-1))

weights1 <- c(1:(last_entry_month-1))
var_returns1 <- c(1:(last_entry_month-1))

weights2 <- c(1:(last_entry_month-1))
var_returns2 <- c(1:(last_entry_month-1))

for (month in 1:(last_entry_month-1)) {
  weights0[month] <- c0/monthly_vars$variance[month]
  weights1[month] <- c1/monthly_vars$ARMA_var[month]
  weights2[month] <- c2/monthly_vars$EWMA_var[month]
  var_returns0[month] <- weights0[month]*FF_monthly$Mkt[month+1]
  var_returns1[month] <- weights1[month]*FF_monthly$Mkt[month+1]
  var_returns2[month] <- weights2[month]*FF_monthly$Mkt[month+1]
}
returns <- data.frame(FF_monthly$Date[-1],FF_monthly$Mkt[-1], 
                      var_returns0, var_returns1, var_returns2)
colnames(returns) <- c("Date", "Mkt", "vol_managed", "ARMA_vol_managed", "EWMA_vol_managed")


# ***** Some descriptive statistics *****
print(var(returns$Mkt))
print(var(returns$vol_managed))
print(var(returns$ARMA_vol_managed))
print(var(returns$EWMA_vol_managed))

print(quantile(weights0, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
print(quantile(weights1, probs = c(0.5, 0.75, 0.9, 0.99)))
print(quantile(weights2, probs = c(0.5, 0.75, 0.9, 0.99)))

# ***** Calculate performance / total returns *****
tot_ret_Mkt = c(1:last_entry_month)
tot_ret0 = c(1:last_entry_month)
tot_ret1 = c(1:last_entry_month)
tot_ret2 = c(1:last_entry_month)

for (month in 2:last_entry_month) {
  tot_ret_Mkt[month] = tot_ret_Mkt[month - 1]*(1 + returns$Mkt[month-1]/100)
  tot_ret0[month] = tot_ret0[month - 1]*(1 + returns$vol_managed[month - 1]/100)
  tot_ret1[month] = tot_ret1[month - 1]*(1 + returns$ARMA_vol_managed[month - 1]/100)
  tot_ret2[month] = tot_ret2[month - 1]*(1 + returns$EWMA_vol_managed[month - 1]/100)
}

performance <- data.frame(FF_monthly$Date, tot_ret_Mkt, tot_ret0, tot_ret1, tot_ret2)
colnames(performance) <- c("Date", "Mkt", "vol_managed", "ARMA_vol_managed", "EWMA_vol_managed")
performance[last_entry_month,]

# ***** Plot market and VM returns on log scale *****
months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(performance, aes(months)) +
  geom_line(aes(y=Mkt)) +
  geom_line(aes(y=vol_managed)) +
  geom_line(aes(y=ARMA_vol_managed)) +
  geom_line(aes(y=EWMA_vol_managed)) +
  scale_x_date(limits = as.Date(c("1926-7-1", "2015-4-1")),
               expand = c(0,0),
               breaks = dates,
               labels = format(dates, "%Y"),
               minor_breaks = NULL) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks('log10', function(x) 10^x),
                     minor_breaks = scale,
                     labels = trans_format('log10', math_format(10^.x)),
                     limits = c(0.1,100000),
                     expand = c(0,0)) +
  ggtitle("Cumulative Performance") + xlab("") + ylab("")


# ***** Perform regressions to obtain alpha *****
reg <- lm(returns$vol_managed ~ returns$Mkt)
reg2 <- lm(returns$vol_managed ~ FF_monthly$Mkt[2:1066] + 
             FF_monthly$SMB[2:1066] + FF_monthly$HML[2:1066])
print(alpha <- reg$coefficients[1]*12) # paper: 4.86
print(alpha <- reg2$coefficients[1]*12) # paper: 5.45

