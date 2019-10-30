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
  summarise(variance = var(`Mkt-RF`) * 21)

monthly_vars <- monthly_vars %>% mutate(volatility = sqrt(variance))
plot(monthly_vars$volatility * sqrt(trading_months), type = "l")

# ******** Scenario 2: ARIMA model**

variance_ts <- xts(monthly_vars$variance, order.by = FF_monthly$Date)
min_obs <- 12
max_obs <- 100
model <- vector(mode = "list", length = last_entry_month - min_obs)

monthly_vars$ARMA_var <- c(1:last_entry_month)
for (i in 1:last_entry_month) {
  if (i <= min_obs) monthly_vars$ARMA_var[i] <- variance_ts[i]
  else {
    model[[i-min_obs]] <- auto.arima(variance_ts[max(1,(i-max_obs)):i], d = 0)
    monthly_vars$ARMA_var[i] <- as.numeric(forecast(model[[i-min_obs]], h = 1)$mean)
  }
}

p <- vector(length = last_entry_month - min_obs)
q <- vector(length = last_entry_month - min_obs)
for (i in 1:(last_entry_month - min_obs)) {
  p[i] <- arimaorder(model[[i]])[1]
  q[i] <- arimaorder(model[[i]])[3]
}
plot(p, type = "l")
plot(q, type = "l")

# ***** Scenario 3: EWMA ****

n <- 0.94
monthly_vars$EWMA_var <- c(1:last_entry_month)
for (i in 1:last_entry_month) {
  if (i == 1) monthly_vars$EWMA_var[i] <- monthly_vars$variance[i]
  else {
    monthly_vars$EWMA_var[i] <- 
      n * monthly_vars$EWMA_var[i - 1] + (1 - n) * monthly_vars$variance[i]
  }
}

# ***** 4: Var of Var****
var_of_var_periods <- 12

monthly_vars$var_of_var <- c(1:last_entry_month)
for (i in 1:last_entry_month) {
  if (i <= var_of_var_periods) {
    monthly_vars$var_of_var[i] <- 
      var(monthly_vars$variance[1:var_of_var_periods])
  }
  else {
    monthly_vars$var_of_var[i] <- 
      var(monthly_vars$variance[(i-var_of_var_periods + 1):i])
  }
}

# **** 5: Vol of vol
vol_of_vol_periods <- 12

monthly_vars$vol_of_vol <- c(1:last_entry_month)
for (i in 1:last_entry_month) {
  if (i <= vol_of_vol_periods) {
    monthly_vars$vol_of_vol[i] <- sd(monthly_vars$volatility[1:vol_of_vol_periods])
  }
  else {
    monthly_vars$vol_of_vol[i] <- sd(monthly_vars$volatility[(i-vol_of_vol_periods + 1):i])
  }
}

#*** 6: recent emphasis
emphasized_days <- 10

recent_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(tail(Mkt, emphasized_days)) * trading_days)

monthly_vars$recent_var <- recent_vars$variance
monthly_vars <- monthly_vars %>% mutate(recent_vol = sqrt(recent_var))

# initiate strategy name vector
names <- c("var_managed", "vol_managed", "ARMA_vol_managed", "EWMA_vol_managed", 
           "var_of_var", "vol_of_vol", "recent_emphasize")

# decide denominator
denom <- data.frame(matrix(ncol = length(names), nrow = last_entry_month - 1))
colnames(denom) <- names

denom[,1] <- monthly_vars$variance[-last_entry_month]
denom[,2] <- monthly_vars$volatility[-last_entry_month]
denom[,3] <- monthly_vars$ARMA_var[-last_entry_month]
denom[,4] <- monthly_vars$EWMA_var[-last_entry_month]
denom[,5] <- monthly_vars$variance[-last_entry_month] * 
  monthly_vars$var_of_var[-last_entry_month]
denom[,6] <- monthly_vars$volatility[-last_entry_month] * 
  monthly_vars$vol_of_vol[-last_entry_month]
denom[,7] <- monthly_vars$volatility[-last_entry_month] *
  monthly_vars$recent_vol[-last_entry_month]

# ***** Calculate c *****
c <- data.frame(matrix(ncol = length(names)))
colnames(c) <- names

#************************************************************************************************************
# new c test (2 solutions due to quadratic equation)
c1 = 1/(2*var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month]))*
  (-2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1])+
     sqrt((2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1]))^2-4*
            var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month])*
            (var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1]))))
c2 = 1/(2*var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month]))*
  (-2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1])-
     sqrt((2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1]))^2-4*
            var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month])*
            (var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1]))))

weights <- c(1:1065)
vola_managed_returns <- c(1:1065)

for (month in 2:1066) {
  weights[month-1] <- c1/monthly_vars$variance[month-1]
  vola_managed_returns[month-1] <- weights[month-1]*(FF_monthly$Mkt[month]-FF_monthly$RF[month])+FF_monthly$RF[month]
}
returns <- data.frame(FF_monthly$Mkt[2:1066], vola_managed_returns)
colnames(returns) <- c("Market Returns", "Vola Managed Returns")


print(var(vola_managed_returns))
print(var(FF_monthly$Mkt[2:1066]))
print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39

tot_ret = c(1:1066)
tot_ret_VM = c(1:1066)

for (month in 2:1066) {
  tot_ret[month] = tot_ret[month - 1]*(1 + FF_monthly$Mkt[month]/100)
  tot_ret_VM[month] = tot_ret_VM[month - 1]*(1 + vola_managed_returns[month - 1]/100)
}
tot_ret[1066]
tot_ret_VM[1066]

#********************************************************************************************




for (i in 1:length(names)) {
  c[i] <- sqrt(var(FF_monthly$`Mkt-RF`[-1]) / var(FF_monthly$`Mkt-RF`[-1] / denom[,i]))
}

# ***** Calculate weights and volatility managed returns *****
weights <- data.frame(matrix(ncol = length(names), nrow = last_entry_month - 1))
colnames(weights) <- names

returns <- data.frame(matrix(ncol = length(names) + 3, nrow = last_entry_month - 1))
colnames(returns) <- c("Date", "Mkt", "rf", names)
returns <- returns %>% mutate(Date = FF_monthly$Date[-1], Mkt = FF_monthly$Mkt[-1], 
                              rf = FF_monthly$RF[-1])

for (month in 1:(last_entry_month-1)) {
  for (j in 1:length(names)) {
    weights[month,j] <- c[1,j] / denom[month,j]
    returns[month,names[j]] <- weights[month, j] * (returns$Mkt[month] - returns$rf[month]) +
      returns$rf[month]
  }
}

# ***** Some descriptive statistics *****
print(apply(returns[,-c(1,3)], 2, var))
print(apply(weights, 2, quantile, probs = c(0.5, 0.75, 0.9, 0.99)))
# paper: 0.93 1.59 2.64 6.39

# ***** Calculate performance / total returns *****
tot_ret <- data.frame(matrix(ncol = length(names) + 2, nrow = last_entry_month))
colnames(tot_ret) <- c("Date", "Mkt", names)
tot_ret <- tot_ret %>% mutate(Date = FF_monthly$Date)

tot_ret[1,-1] <- 1

for (month in 2:last_entry_month) {
  tot_ret$Mkt[month] <- tot_ret$Mkt[month-1] * 
    (1 + (returns$Mkt[month-1]/100))
  for (j in 1:length(names)) {
    tot_ret[month,names[j]] <- tot_ret[month-1, names[j]] * 
      (1 + (returns[month-1, names[j]]/100))
  }
}

tot_ret[last_entry_month,]

# ***** Plot market and VM returns on log scale *****
months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(tot_ret, aes(months)) +
  geom_line(aes(y=Mkt)) +
  geom_line(aes(y=var_managed)) +
  geom_line(aes(y=ARMA_vol_managed)) +
  geom_line(aes(y=EWMA_vol_managed)) +
  geom_line(aes(y=var_of_var)) +
  geom_line(aes(y=vol_of_vol)) +
  geom_line(aes(y=recent_emphasize)) +
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

# regressions
a <- 12 * (returns$var_managed - returns$rf)
b <- 12 * (returns$Mkt - returns$rf)
lm(a ~ b)

a1 <- 12 * (returns$vol_managed - returns$rf)
a2 <- 12 * (returns$ARMA_vol_managed - returns$rf)
a3 <- 12 * (returns$EWMA_vol_managed - returns$rf)
a4 <- 12 * (returns$var_of_var - returns$rf)
a5 <- 12 * (returns$vol_of_vol - returns$rf)
a6 <- 12 * (returns$recent_emphasize - returns$rf)

lm(a~b)
lm(a1~b)
lm(a2~b)
lm(a3~b)
lm(a4~b)
lm(a5~b)
lm(a6~b)