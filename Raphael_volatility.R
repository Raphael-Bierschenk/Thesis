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
VIX_daily <- read_csv("VIXCLS.csv", col_names = TRUE, 
                      col_types = cols(Date = "D", VIXCLS = "n"))

FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)
VIX_daily <- VIX_daily %>% rename(Date = DATE, VIX = VIXCLS)

first_day <- 19900101
first_month <- 199001
last_day <- 20190630
last_month <- 201906

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)
VIX_daily <- VIX_daily %>% subset(subset = Date <= ymd(last_day) & Date >= ymd(first_day))

n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))

VIX_monthly <- data.frame(matrix(ncol = ncol(VIX_daily), nrow = n_months))
colnames(VIX_monthly) <- colnames(VIX_daily)
VIX_monthly <- VIX_monthly %>% 
for (i in 2:as.integer(count(VIX_daily))) {
  if (month(VIX_daily$Date[i]) != month(VIX_daily$Date[i-1])) {
    VIX_monthly$VIX[i] <- VIX_daily$VIX[i]
  }
}

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
model <- vector(mode = "list", length = n_months - min_obs)

monthly_vars$ARMA_var <- c(1:n_months)
for (i in 1:n_months) {
  if (i <= min_obs) monthly_vars$ARMA_var[i] <- variance_ts[i]
  else {
    model[[i-min_obs]] <- auto.arima(variance_ts[max(1,(i-max_obs)):i], d = 0)
    monthly_vars$ARMA_var[i] <- as.numeric(forecast(model[[i-min_obs]], h = 1)$mean)
  }
}

p <- vector(length = n_months - min_obs)
q <- vector(length = n_months - min_obs)
for (i in 1:(n_months - min_obs)) {
  p[i] <- arimaorder(model[[i]])[1]
  q[i] <- arimaorder(model[[i]])[3]
}
plot(p, type = "l")
plot(q, type = "l")

# ***** Scenario 3: EWMA ****

n <- 0.94
monthly_vars$EWMA_var <- c(1:n_months)
for (i in 1:n_months) {
  if (i == 1) monthly_vars$EWMA_var[i] <- monthly_vars$variance[i]
  else {
    monthly_vars$EWMA_var[i] <- 
      n * monthly_vars$EWMA_var[i - 1] + (1 - n) * monthly_vars$variance[i]
  }
}

# ***** 4: Var of Var****
var_of_var_periods <- 12

monthly_vars$var_of_var <- c(1:n_months)
for (i in 1:n_months) {
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

monthly_vars$vol_of_vol <- c(1:n_months)
for (i in 1:n_months) {
  if (i <= vol_of_vol_periods) {
    monthly_vars$vol_of_vol[i] <- sd(monthly_vars$volatility[1:vol_of_vol_periods])
  }
  else {
    monthly_vars$vol_of_vol[i] <- sd(monthly_vars$volatility[(i-vol_of_vol_periods + 1):i])
  }
}

#*** 6: recent emphasis
emphasized_days <- 5

recent_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(tail(Mkt, emphasized_days)) * trading_days)

monthly_vars$recent_var <- recent_vars$variance
monthly_vars <- monthly_vars %>% mutate(recent_vol = sqrt(recent_var))

# **** additional daily vol of vol factor
monthly_vars$vol_of_vol_d <- c(1:n_months)
j <- 1

for (i in 2:n_days
     ) {
  if (month(FF_daily$Date[i]) != month(FF_daily$Date[i-1])) {
    vol_1 <- sd(FF_daily$`Mkt-RF`[(i - 5):(i - 1)])
    vol_2 <- sd(FF_daily$`Mkt-RF`[(i - 10):(i - 6)])
    vol_3 <- sd(FF_daily$`Mkt-RF`[(i - 15):(i - 11)])
    vol_4 <- sd(FF_daily$`Mkt-RF`[(i - 20):(i - 16)])
    monthly_vars$vol_of_vol_d[j] <- sd(c(vol_1, vol_2, vol_3, vol_4))
    j <- j + 1
  }
  else if (i == n_days
           ) {
    vol_1 <- sd(FF_daily$`Mkt-RF`[(i - 4):i])
    vol_2 <- sd(FF_daily$`Mkt-RF`[(i - 9):(i - 5)])
    vol_3 <- sd(FF_daily$`Mkt-RF`[(i - 14):(i - 10)])
    vol_4 <- sd(FF_daily$`Mkt-RF`[(i - 19):(i - 15)])
    monthly_vars$vol_of_vol_d[j] <- sd(c(vol_1, vol_2, vol_3, vol_4))
  }
  else {}
}

vol_of_vol_weight <- 1.2

vol_of_vol_factor <- c(1:n_months)
for (i in 1:n_months) {
  if (i <= vol_of_vol_periods) vol_of_vol_factor[i] <- 1
  else {
    if (monthly_vars$vol_of_vol_d[i] <= 
        mean(monthly_vars$vol_of_vol_d[(i - vol_of_vol_periods + 1):i])) {
      vol_of_vol_factor[i] <- vol_of_vol_weight
    }
    else {vol_of_vol_factor[i] <- 2 - vol_of_vol_weight}
  }
}

# initiate strategy name vector
names <- c("var_managed", "vol_managed", "ARMA_vol_managed", "EWMA_vol_managed", 
           "var_of_var", "vol_of_vol", "recent_emphasize", "vol_of_vol_d", "test", "test1")

# decide denominator
denom <- data.frame(matrix(ncol = length(names), nrow = n_months - 1))
colnames(denom) <- names

denom[,1] <- monthly_vars$variance[-n_months]
denom[,2] <- monthly_vars$volatility[-n_months]
denom[,3] <- monthly_vars$ARMA_var[-n_months]
denom[,4] <- monthly_vars$EWMA_var[-n_months]
denom[,5] <- monthly_vars$variance[-n_months] * 
  monthly_vars$var_of_var[-n_months]
denom[,6] <- monthly_vars$volatility[-n_months] * 
  monthly_vars$vol_of_vol[-n_months]
denom[,7] <- monthly_vars$volatility[-n_months] *
  monthly_vars$recent_vol[-n_months]
denom[,8] <- monthly_vars$variance[-n_months] *
  vol_of_vol_factor[-n_months]
denom[,9] <- monthly_vars$variance[-n_months] / 
  monthly_vars$vol_of_vol_d[-n_months]
denom[,10] <- monthly_vars$recent_var[-n_months]


# ***** Calculate c *****
c <- data.frame(matrix(ncol = length(names)))
colnames(c) <- names

for (i in 1:length(names)) {
  c[i] <- sqrt(var(FF_monthly$`Mkt-RF`[-1]) / var(FF_monthly$`Mkt-RF`[-1] / denom[,i]))
}

a_qe <- c(1:length(names))
b_qe <- c(1:length(names))
c_qe <- c(1:length(names))

for (i in 1:length(names)) {
  # determine a,b,c of midnight formula
  a_qe[i] <- var(FF_monthly$`Mkt-RF`[-1]/denom[,i])
  b_qe[i] <- 2*cov(FF_monthly$`Mkt-RF`[-1]/denom[,i], FF_monthly$RF[-1])
  c_qe[i] <- var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1])
  
  # apply midnight formula (we always take the x1 solution)
  c[i] <- 1/(2*a_qe[i])*(-b_qe[i]+sqrt((b_qe[i])^2-4*a_qe[i]*c_qe[i]))
}





#************************************************************************************************************
# new c test (2 solutions due to quadratic equation)
#c1 = 1/(2*var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month]))*
#  (-2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1])+
#     sqrt((2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1]))^2-4*
#            var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month])*
#            (var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1]))))
#c2 = 1/(2*var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month]))*
#  (-2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1])-
#     sqrt((2*cov(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month], FF_monthly$RF[-1]))^2-4*
#            var(FF_monthly$`Mkt-RF`[-1]/monthly_vars$variance[-last_entry_month])*
#            (var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1]))))
# 
# weights <- c(1:1065)
# vola_managed_returns <- c(1:1065)
# 
# for (month in 2:1066) {
#   weights[month-1] <- c1/monthly_vars$variance[month-1]
#   vola_managed_returns[month-1] <- weights[month-1]*(FF_monthly$Mkt[month]-FF_monthly$RF[month])+FF_monthly$RF[month]
# }
# returns <- data.frame(FF_monthly$Mkt[2:1066], vola_managed_returns)
# colnames(returns) <- c("Market Returns", "Vola Managed Returns")
# 
# 
# print(var(vola_managed_returns))
# print(var(FF_monthly$Mkt[2:1066]))
# print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
# 
# tot_ret = c(1:1066)
# tot_ret_VM = c(1:1066)
# 
# for (month in 2:1066) {
#   tot_ret[month] = tot_ret[month - 1]*(1 + FF_monthly$Mkt[month]/100)
#   tot_ret_VM[month] = tot_ret_VM[month - 1]*(1 + vola_managed_returns[month - 1]/100)
# }
# tot_ret[1066]
# tot_ret_VM[1066]
# 
#********************************************************************************************





# ***** Calculate weights and volatility managed returns *****
weights <- data.frame(matrix(ncol = length(names), nrow = n_months - 1))
colnames(weights) <- names

returns <- data.frame(matrix(ncol = length(names) + 3, nrow = n_months - 1))
colnames(returns) <- c("Date", "Mkt", "rf", names)
returns <- returns %>% mutate(Date = FF_monthly$Date[-1], Mkt = FF_monthly$Mkt[-1], 
                              rf = FF_monthly$RF[-1])

for (month in 1:(n_months-1)) {
  for (j in 1:length(names)) {
    weights[month,j] <- c[1,j] / denom[month,j]
    returns[month,names[j]] <- weights[month, j] * (returns$Mkt[month] - returns$rf[month]) +
      returns$rf[month]
  }
}

# ***** Some descriptive statistics *****
print(apply(returns[,-c(1,3)], 2, var))
print(apply(weights, 2, quantile, probs = c(0.5, 0.75, 0.9, 0.99)))
# paper: 0.93 1.59 2.64 6.39 for var managed

# ***** Calculate performance / total returns *****
tot_ret <- data.frame(matrix(ncol = length(names) + 2, nrow = n_months))
colnames(tot_ret) <- c("Date", "Mkt", names)
tot_ret <- tot_ret %>% mutate(Date = FF_monthly$Date)

tot_ret[1,-1] <- 1

for (month in 2:n_months) {
  tot_ret$Mkt[month] <- tot_ret$Mkt[month-1] * 
    (1 + (returns$Mkt[month-1]/100))
  for (j in 1:length(names)) {
    tot_ret[month,names[j]] <- tot_ret[month-1, names[j]] * 
      (1 + (returns[month-1, names[j]]/100))
  }
}

tot_ret[n_months,]

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
reg_mkt <- vector(mode = "list", length = length(names))
reg_FF3 <- vector(mode = "list", length = length(names))
b <- 12 * (returns$Mkt - returns$rf)
b1 <- 12 * (FF_monthly$SMB[-1])
b2 <- 12 * (FF_monthly$HML[-1])

for (i in 1:length(names)) {
  a <- 12 * (returns[, names[i]] - returns$rf)
  reg_mkt[[i]] <- lm(a ~ b)
}

for (i in 1:length(names)) {
  a <- 12 * (returns[, names[i]] - returns$rf)
  reg_FF3[[i]] <- lm(a ~ b + b1 + b2)
}


output_names <- c("alpha_mkt", "R^2_mkt", "RMSE", "SR", "Appr_Ratio", "alpha_FF3")
reg_output <- data.frame(matrix(ncol = length(names), nrow = length(output_names)))
colnames(reg_output) <- names
rownames(reg_output) <- output_names

for (i in 1:length(names)) {
  reg_output["alpha_mkt", i] <- reg_mkt[[i]]$coefficients[1]
  reg_output["R^2_mkt", i] <- summary(reg_mkt[[i]])$r.squared
  reg_output["RMSE", i] <- sigma(reg_mkt[[i]])
  reg_output["SR", i] <- 12 * (mean(returns[,names[i]] - returns$rf)) / 
    (sqrt(trading_months) * sd(returns[,names[i]]))
  reg_output["Appr_Ratio", i] <- sqrt(trading_months) * reg_output["alpha_mkt", i] /
    reg_output["RMSE", i]
  reg_output["alpha_FF3", i] <- reg_FF3[[i]]$coefficients[1]
}

round(reg_output, 2)
