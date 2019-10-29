library(readr)
library(ggplot2)
library(scales)
library(dplyr)
library(lubridate)
library(tidyverse)
library(glue)

# ***** Import Data *****
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)

FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)

# Adjust end date --> end of 2015
# @Stefan: why April? only says "2015" in paper

last_day <- 20150330
last_month <- 201504

last_entry_day <- which(FF_daily[,1] == last_day,)
last_entry_month <- which(FF_monthly[,1] == last_month,)
FF_daily <- head(FF_daily, last_entry_day)
FF_monthly <- head(FF_monthly, last_entry_month)

# ***** Manipulate Data *****
FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")


# ***** Volatility Estimating *****
FF_daily$u <- log(1+FF_daily$Mkt/100)
FF_daily$u_sq <- FF_daily$u^2
plot(FF_daily$u, type = "l")
mean(FF_daily$u) # mean of log changes basically zero as required




# ****************************************
# ***** Volatility Estimating - EWMA *****
# ****************************************

lambda = 0.94
FF_daily$EWMA_vars <- c(1:nrow(FF_daily))
FF_daily$EWMA_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$EWMA_vars[i] <- lambda*FF_daily$EWMA_vars[i-1]+(1-lambda)*FF_daily$u_sq[i]
}
plot(FF_daily$EWMA_vars, type = "l")


trading_days <- 250
monthly_vars_EWMA <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = mean(EWMA_vars) * trading_days * 10000)

monthly_vars_EWMA <- monthly_vars_EWMA %>% mutate(volatility = sqrt(variance))
plot(monthly_vars_EWMA$volatility, type = "l")

c_EWMA = sqrt(var(FF_monthly$Mkt[2:1066])/var(1/monthly_vars_EWMA$variance[1:1065]*FF_monthly$Mkt[2:1066]))

weights_EWMA <- c(1:1065)
vola_managed_returns_EWMA <- c(1:1065)

for (month in 2:1066) {
  weights_EWMA[month-1] <- c_EWMA/monthly_vars_EWMA$variance[month-1]
  vola_managed_returns_EWMA[month-1] <- weights_EWMA[month-1]*FF_monthly$Mkt[month]
}
returns <- data.frame(FF_monthly$Mkt[2:1066], vola_managed_returns_EWMA)
colnames(returns) <- c("Market Returns", "Vola Managed Returns EWMA")

print(var(vola_managed_returns_EWMA))
print(var(FF_monthly$Mkt[2:1066]))
print(quantile(weights_EWMA, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
plot(weights_EWMA)
summary(FF_monthly$Mkt)
summary(vola_managed_returns_EWMA)

tot_ret_EWMA = c(1:1066)
tot_ret_VM_EWMA = c(1:1066)

for (month in 2:1066) {
  tot_ret_EWMA[month] = tot_ret_EWMA[month - 1]*(1 + FF_monthly$Mkt[month]/100)
  tot_ret_VM_EWMA[month] = tot_ret_VM_EWMA[month - 1]*(1 + vola_managed_returns_EWMA[month - 1]/100)
}
tot_ret_EWMA[1066]
tot_ret_VM_EWMA[1066]

months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
performance <- data.frame(months, tot_ret_EWMA, tot_ret_VM_EWMA)

scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(performance, aes(months)) +
  geom_line(aes(y=tot_ret_EWMA)) +
  geom_line(aes(y=tot_ret_VM_EWMA)) +
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

reg_EWMA <- lm(vola_managed_returns_EWMA ~ FF_monthly$Mkt[2:1066])
reg2_EWMA <- lm(vola_managed_returns_EWMA ~ FF_monthly$Mkt[2:1066] + FF_monthly$SMB[2:1066] + FF_monthly$HML[2:1066])
print(alpha <- reg_EWMA$coefficients[1]*12) # paper: 4.86
print(alpha <- reg2_EWMA$coefficients[1]*12) # paper: 5.45




# *****************************************
# ***** Volatility Estimating - GARCH *****
# *****************************************

omega = 0.000001
alpha = 0.074715
beta = 0.9175
FF_daily$GARCH_vars <- c(1:nrow(FF_daily))
FF_daily$GARCH_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$GARCH_vars[i] <- omega+alpha*FF_daily$u_sq[i]+beta*FF_daily$GARCH_vars[i-1]
}
plot(FF_daily$GARCH_vars, type = "l")


trading_days <- 250
monthly_vars_GARCH <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = mean(GARCH_vars) * trading_days * 10000)

monthly_vars_GARCH <- monthly_vars_GARCH %>% mutate(volatility = sqrt(variance))
plot(monthly_vars_GARCH$volatility, type = "l")

c_GARCH = sqrt(var(FF_monthly$Mkt[2:1066])/var(1/monthly_vars_GARCH$variance[1:1065]*FF_monthly$Mkt[2:1066]))

weights_GARCH <- c(1:1065)
vola_managed_returns_GARCH <- c(1:1065)

for (month in 2:1066) {
  weights_GARCH[month-1] <- c_GARCH/monthly_vars_GARCH$variance[month-1]
  vola_managed_returns_GARCH[month-1] <- weights_GARCH[month-1]*FF_monthly$Mkt[month]
}
returns <- data.frame(FF_monthly$Mkt[2:1066], vola_managed_returns_GARCH)
colnames(returns) <- c("Market Returns", "Vola Managed Returns GARCH")

print(var(vola_managed_returns_GARCH))
print(var(FF_monthly$Mkt[2:1066]))
print(quantile(weights_GARCH, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
plot(weights_GARCH)
summary(FF_monthly$Mkt)
summary(vola_managed_returns_GARCH)

tot_ret_GARCH = c(1:1066)
tot_ret_VM_GARCH = c(1:1066)

for (month in 2:1066) {
  tot_ret_GARCH[month] = tot_ret_GARCH[month - 1]*(1 + FF_monthly$Mkt[month]/100)
  tot_ret_VM_GARCH[month] = tot_ret_VM_GARCH[month - 1]*(1 + vola_managed_returns_GARCH[month - 1]/100)
}
tot_ret_GARCH[1066]
tot_ret_VM_GARCH[1066]

months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
performance <- data.frame(months, tot_ret_GARCH, tot_ret_VM_GARCH)

scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(performance, aes(months)) +
  geom_line(aes(y=tot_ret_GARCH)) +
  geom_line(aes(y=tot_ret_VM_GARCH)) +
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

reg_GARCH <- lm(vola_managed_returns_GARCH ~ FF_monthly$Mkt[2:1066])
reg2_GARCH <- lm(vola_managed_returns_GARCH ~ FF_monthly$Mkt[2:1066] + FF_monthly$SMB[2:1066] + FF_monthly$HML[2:1066])
print(alpha <- reg_GARCH$coefficients[1]*12) # paper: 4.86
print(alpha <- reg2_GARCH$coefficients[1]*12) # paper: 5.45



# ***** Calculate Monthly Variances *****
trading_days <- 250
monthly_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(Mkt) * trading_days)

monthly_vars <- monthly_vars %>% mutate(volatility = sqrt(variance))
plot(monthly_vars$volatility, type = "l")


# *******************************************************
# ***** Historical Simulation - Volatility Updating *****
# *******************************************************


# Default + EWMA + GARCH
FF_monthly$MktVUp <- c(1:nrow(FF_monthly)-1)
FF_monthly$MktVUpEWMA <- c(1:nrow(FF_monthly)-1)
FF_monthly$MktVUpGARCH <- c(1:nrow(FF_monthly)-1)

for (i in 2:nrow(FF_monthly)) {
  FF_monthly$MktVUp[i] <- FF_monthly$Mkt[i]*
    sqrt(monthly_vars$variance[nrow(FF_monthly)-1]/monthly_vars$variance[i-1])
  FF_monthly$MktVUpEWMA[i] <- FF_monthly$Mkt[i]*
    sqrt(monthly_vars_EWMA$variance[nrow(FF_monthly)-1]/monthly_vars_EWMA$variance[i-1])
  FF_monthly$MktVUpGARCH[i] <- FF_monthly$Mkt[i]*
    sqrt(monthly_vars_GARCH$variance[nrow(FF_monthly)-1]/monthly_vars_GARCH$variance[i-1])
}

# ***** Calculate c *****
c_VUp = sqrt(var(FF_monthly$Mkt[2:1066])/var(1/monthly_vars$variance[1:1065]*FF_monthly$Mkt[2:1066]))
c_VUp_EWMA = sqrt(var(FF_monthly$Mkt[2:1066])/var(1/monthly_vars_EWMA$variance[1:1065]*FF_monthly$Mkt[2:1066]))
c_VUp_GARCH = sqrt(var(FF_monthly$Mkt[2:1066])/var(1/monthly_vars_GARCH$variance[1:1065]*FF_monthly$Mkt[2:1066]))


# ***** Calculate weights and volatility managed returns *****
weights_VUp <- c(1:1065)
weights_VUp_EWMA <- c(1:1065)
weights_VUp_GARCH <- c(1:1065)
vola_managed_returns_VUp <- c(1:1065)
vola_managed_returns_VUp_EWMA <- c(1:1065)
vola_managed_returns_VUp_GARCH <- c(1:1065)

for (month in 2:1066) {
  weights_VUp[month-1] <- c_VUp/monthly_vars$variance[month-1]
  weights_VUp_EWMA[month-1] <- c_VUp_EWMA/monthly_vars_EWMA$variance[month-1]
  weights_VUp_GARCH[month-1] <- c_VUp_GARCH/monthly_vars_GARCH$variance[month-1]
  vola_managed_returns_VUp[month-1] <- weights_VUp[month-1]*FF_monthly$Mkt[month]
  vola_managed_returns_VUp_EWMA[month-1] <- weights_VUp_EWMA[month-1]*FF_monthly$Mkt[month]
  vola_managed_returns_VUp_GARCH[month-1] <- weights_VUp_GARCH[month-1]*FF_monthly$Mkt[month]
}
returns_VUp <- data.frame(FF_monthly$Mkt[2:1066], vola_managed_returns_VUp, vola_managed_returns_VUp_EWMA, vola_managed_returns_VUp_GARCH)
colnames(returns_VUp) <- c("Market Returns", "VM Returns VUp", "VM Returns VUp + EWMA", "VM Returns VUp + GARCH")


# ***** Some descriptive statistics *****
print(var(FF_monthly$Mkt[2:1066]))
print(var(vola_managed_returns_VUp))
print(var(vola_managed_returns_VUp_EWMA))
print(var(vola_managed_returns_VUp_GARCH))
print(quantile(weights_VUp, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
print(quantile(weights_VUp_EWMA, probs = c(0.5, 0.75, 0.9, 0.99)))
print(quantile(weights_VUp_GARCH, probs = c(0.5, 0.75, 0.9, 0.99)))


# ***** Calculate performance / total returns *****
tot_ret = c(1:1066)
tot_ret_VM_VUp = c(1:1066)
tot_ret_VM_VUp_EWMA = c(1:1066)
tot_ret_VM_VUp_GARCH = c(1:1066)

for (month in 2:1066) {
  tot_ret[month] = tot_ret[month - 1]*(1 + FF_monthly$Mkt[month]/100)
  tot_ret_VM_VUp[month] = tot_ret_VM_VUp[month - 1]*(1 + vola_managed_returns_VUp[month - 1]/100)
  tot_ret_VM_VUp_EWMA[month] = tot_ret_VM_VUp_EWMA[month - 1]*(1 + vola_managed_returns_VUp_EWMA[month - 1]/100)
  tot_ret_VM_VUp_GARCH[month] = tot_ret_VM_VUp_GARCH[month - 1]*(1 + vola_managed_returns_VUp_GARCH[month - 1]/100)
}

months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
performance_VUp <- data.frame(months, tot_ret, tot_ret_VM_VUp, tot_ret_VM_VUp_EWMA, tot_ret_VM_VUp_GARCH)
performance_VUp[last_entry_month,]


# ***** Plot market and VM returns on log scale *****
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(performance, aes(months)) +
  geom_line(aes(y=tot_ret)) +
  geom_line(aes(y=tot_ret_VM_VUp)) +
  geom_line(aes(y=tot_ret_VM_VUp_EWMA)) +
  geom_line(aes(y=tot_ret_VM_VUp_GARCH)) +
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







# ************************** default *****************************************

# ***** Calculate c *****
c = sqrt(var(FF_monthly$Mkt[2:1066])/var(1/monthly_vars$variance[1:1065]*FF_monthly$Mkt[2:1066]))
# Tried to hardcode c to obtain weight quantiles from paper but distribution itself seems to be unequal
# -> Even if c is wrong, the problem occurs earlier


# ***** Calculate weights and volatility managed returns *****
weights <- c(1:1065)
vola_managed_returns <- c(1:1065)

for (month in 2:1066) {
  weights[month-1] <- c/monthly_vars$variance[month-1]
  vola_managed_returns[month-1] <- weights[month-1]*FF_monthly$Mkt[month]
}
returns <- data.frame(FF_monthly$Mkt[2:1066], vola_managed_returns)
colnames(returns) <- c("Market Returns", "Vola Managed Returns")


# ***** Some descriptive statistics *****
print(var(vola_managed_returns))
print(var(FF_monthly$Mkt[2:1066]))
print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
plot(weights)
summary(FF_monthly$Mkt)
summary(vola_managed_returns)
sort(weights, decreasing = TRUE)[1:20]


# ***** Calculate performance / total returns *****
tot_ret = c(1:1066)
tot_ret_VM = c(1:1066)

for (month in 2:1066) {
  tot_ret[month] = tot_ret[month - 1]*(1 + FF_monthly$Mkt[month]/100)
  tot_ret_VM[month] = tot_ret_VM[month - 1]*(1 + vola_managed_returns[month - 1]/100)
}
tot_ret[1066]
tot_ret_VM[1066]

months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
performance <- data.frame(months, tot_ret, tot_ret_VM)


# ***** Plot market and VM returns on log scale *****
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(performance, aes(months)) +
  geom_line(aes(y=tot_ret)) +
  geom_line(aes(y=tot_ret_VM)) +
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
reg <- lm(vola_managed_returns ~ FF_monthly$Mkt[2:1066])
reg2 <- lm(vola_managed_returns ~ FF_monthly$Mkt[2:1066] + FF_monthly$SMB[2:1066] + FF_monthly$HML[2:1066])
print(alpha <- reg$coefficients[1]*12) # paper: 4.86
print(alpha <- reg2$coefficients[1]*12) # paper: 5.45

