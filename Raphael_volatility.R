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

first_day <- 19260701
first_month <- 192607
last_day <- 20150430
last_month <- 201504

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)
VIX_daily <- VIX_daily %>% subset(subset = Date <= ymd(last_day) & Date >= ymd(first_day))

n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))



# VIX_monthly <- data.frame(matrix(ncol = ncol(VIX_daily), nrow = n_months))
# colnames(VIX_monthly) <- colnames(VIX_daily)
# VIX_monthly <- VIX_monthly %>% 
# for (i in 2:as.integer(count(VIX_daily))) {
#   if (month(VIX_daily$Date[i]) != month(VIX_daily$Date[i-1])) {
#     VIX_monthly$VIX[i] <- VIX_daily$VIX[i]
#   }
# }





# 
# # ***** Volatility Estimating *****
# FF_daily$u <- log(1+FF_daily$Mkt/100)
# FF_daily$u_sq <- FF_daily$u^2
# plot(FF_daily$u, type = "l")
# mean(FF_daily$u) # mean of log changes basically zero as required
# 
# # ***** Volatility Estimating - EWMA *****
# lambda = 0.94
# FF_daily$EWMA_vars <- c(1:nrow(FF_daily))
# FF_daily$EWMA_vars[1] <- FF_daily$u_sq[1]
# for (i in 2:nrow(FF_daily)) {
#   FF_daily$EWMA_vars[i] <- lambda*FF_daily$EWMA_vars[i-1]+(1-lambda)*FF_daily$u_sq[i]
# }
# 
# # ***** Volatility Estimating - GARCH *****
# omega = 0.000001
# alpha = 0.074715
# beta = 0.9175
# FF_daily$GARCH_vars <- c(1:nrow(FF_daily))
# FF_daily$GARCH_vars[1] <- FF_daily$u_sq[1]
# for (i in 2:nrow(FF_daily)) {
#   FF_daily$GARCH_vars[i] <- omega+alpha*FF_daily$u_sq[i]+beta*FF_daily$GARCH_vars[i-1]
# }
# 
# test = 1
# for (i in c(1:nrow(FF_daily))) {
#   test = (1+FF_daily$Mkt[i]/100)*test
# }
# test
# 
# vars_flexible <- data.frame(0, 0, 0, 0, 0)
# colnames(vars_flexible) <- c("Index", "Variance", "Volatility", "Mkt", "RF")
# min_days_for_var = 10
# max_days_for_var = 22
# deviation = 0.6
# i = 1
# last_i = 0
# while (i <= nrow(FF_daily)) { # nrow(FF_daily)
#   if (i - last_i < min_days_for_var) {
#     i = i + 1
#   }
#   else {
#     if (vars_flexible$Index[1] == 0) {
#       ret_temp = 0
#       rf_temp = 0
#       for (j in c((last_i+1):i)) {
#         ret_temp <- ((1 + FF_daily$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
#         rf_temp <- ((1 + FF_daily$RF[j]/100)*(1 + rf_temp/100) - 1)*100
#       }
#       df_temp <- data.frame(i, var(FF_daily$Mkt[(last_i+1):i]), sd(FF_daily$Mkt[(last_i+1):i]), ret_temp, rf_temp)
#       vars_flexible[nrow(vars_flexible),] <- df_temp
#       last_i = i
#     }
#     else {
#       last_vola = vars_flexible$Volatility[nrow(vars_flexible)]
#       if (abs((sd(FF_daily$Mkt[(last_i+1):i]) - last_vola)/last_vola) >= deviation || i - last_i >= max_days_for_var) {
#         ret_temp = 0
#         rf_temp = 0
#         for (j in c((last_i+1):i)) {
#           ret_temp <- ((1 + FF_daily$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
#           rf_temp <- ((1 + FF_daily$RF[j]/100)*(1 + rf_temp/100) - 1)*100
#         }
#         df_temp <- data.frame(i, var(FF_daily$Mkt[(last_i+1):i]), sd(FF_daily$Mkt[(last_i+1):i]), ret_temp, rf_temp)
#         vars_flexible[nrow(vars_flexible) + 1,] <- df_temp
#         last_i = i
#       }
#     }
#     i = i + 1
#   }
# }
# ret_temp = 0
# for (j in c((last_i+1):(i-1))) {
#   ret_temp <- ((1 + FF_daily$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
#   rf_temp <- ((1 + FF_daily$RF[j]/100)*(1 + rf_temp/100) - 1)*100
# }
# df_temp <- data.frame(i, var(FF_daily$Mkt[(last_i+1):i]), sd(FF_daily$Mkt[(last_i+1):i]), ret_temp, rf_temp)
# vars_flexible[nrow(vars_flexible) + 1,] <- df_temp
# 
# diff_test <- diff(vars_flexible$Index)
# summary(diff_test)
# sum(diff_test==10)
# 
# vars_flexible$`Mkt-RF` <- vars_flexible$Mkt - vars_flexible$RF
# 
# n_periods <- nrow(vars_flexible)
# c_flex <- sqrt(var(vars_flexible$`Mkt-RF`[-1]) / var(vars_flexible$`Mkt-RF`[-1] / vars_flexible$Variance[-n_periods]))
# 
# weights <- c(1:(n_periods-1))
# vola_managed_returns <- c(1:(n_periods-1))
# 
# for (period in 2:n_periods) {
#   weights[period-1] <- c_flex/vars_flexible$Variance[period-1]
#   vola_managed_returns[period-1] <- weights[period-1]*
#     (vars_flexible$Mkt[period]-vars_flexible$RF[period])+vars_flexible$RF[period]
# }
# returns <- data.frame(vars_flexible$Mkt[2:n_periods], vola_managed_returns)
# colnames(returns) <- c("Market Returns", "Vola Managed Returns")
# 
# plot(vars_flexible$Variance, type = "l")
# plot(weights)
# print(var(vola_managed_returns))
# print(var(vars_flexible$Mkt[2:n_periods]))
# print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
# 
# tot_ret = c(1:n_periods)
# tot_ret_VM = c(1:n_periods)
# 
# for (period in 2:n_periods) {
#   tot_ret[period] = tot_ret[period - 1]*(1 + vars_flexible$Mkt[period]/100)
#   tot_ret_VM[period] = tot_ret_VM[period - 1]*(1 + vola_managed_returns[period - 1]/100)
# }
# tot_ret[n_periods]
# tot_ret_VM[n_periods]
# 
# 
# # ***** Plot market and VM returns on log scale *****
# performance <- data.frame(vars_flexible$Index, tot_ret, tot_ret_VM)
# scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
#            200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
#            20000,30000,40000,50000,60000,70000,80000,90000,100000)
# ggplot(performance, aes(vars_flexible$Index)) +
#   geom_line(aes(y=tot_ret)) +
#   geom_line(aes(y=tot_ret_VM)) +
#   scale_y_continuous(trans = "log10",
#                      breaks = trans_breaks('log10', function(x) 10^x),
#                      minor_breaks = scale,
#                      labels = trans_format('log10', math_format(10^.x)),
#                      limits = c(0.1,100000),
#                      expand = c(0,0)) +
#   ggtitle("Cumulative Performance") + xlab("") + ylab("")
# 
# reg_flex <- lm(vola_managed_returns ~ vars_flexible$`Mkt-RF`[-1])
# alpha <- reg_flex$coefficients[1]*12*(nrow(vars_flexible)-1)/1065
# rmse <- sigma(reg_flex)
# SR <- 12 *(nrow(vars_flexible)-1)/1065 * (mean(vola_managed_returns - vars_flexible$RF[-1])) / 
#   (sqrt(12 *(nrow(vars_flexible)-1)/1065) * sd(vola_managed_returns))
# appr <- sqrt(12 *(nrow(vars_flexible)-1)/1065) * alpha / (rmse*22)
# alpha
# rmse
# SR
# appr









# ***** Manipulate Data *****
FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")

VIX_daily <- VIX_daily[VIX_daily$Date %in% FF_daily$Date,]

VIX_monthly <- data.frame(matrix(ncol = ncol(VIX_daily), nrow = n_months))
colnames(VIX_monthly) <- colnames(VIX_daily)
j <- 1
VIX_monthly <- VIX_monthly %>% mutate(Date = FF_monthly$Date)

for (i in 1:as.integer(count(VIX_daily))) {
  if (i == 1) {
    VIX_monthly$VIX[j] <- VIX_daily$VIX[j]
    j <- j + 1
    }
  else if (month(VIX_daily$Date[i]) != month(VIX_daily$Date[i-1])) {
    VIX_monthly$VIX[j] <- ifelse(!is.na(VIX_daily$VIX[i-1]), 
                                 VIX_daily$VIX[i-1], VIX_daily$VIX[i])
    j <- j + 1
  }
}



# test with higher frequency
frequ <- 5
j <- 1
n <- as.integer(n_days/frequ)
daily <- data.frame(matrix(ncol = 4, nrow = n))
colnames(daily) <- c("Date", "Var", "Mkt", "rf")
daily <- daily %>% mutate(Date = FF_daily$Date[1:n])
# workaround to set format
for (i in 1:n_days) {
  if (i %% frequ == 0) {
    daily$Date[j] <- FF_daily$Date[i]
    daily$Var[j] <- max(var(FF_daily$`Mkt-RF`[(i+1-frequ):i]) * 22, 1)
    ret_temp <- 0
    rf_temp <- 0
    for (a in 1:min(frequ, n_days - i)) {
      ret_temp <- ((1 + FF_daily$Mkt[i+a]/100)*(1 + ret_temp/100) - 1)*100
      rf_temp <- ((1 + FF_daily$RF[i+a]/100)*(1 + rf_temp/100) - 1)*100
    }
    daily$Mkt[j] <- ret_temp
    daily$rf[j] <- rf_temp
    j <- j + 1
  }
}
# 
# long <- 5
# short <- 5
# daily1 <- data.frame(matrix(ncol = 4, nrow = n_days - 1))
# colnames(daily1) <- c("Date", "Var", "Mkt", "rf")
# daily1 <- daily1 %>% mutate(Date = FF_daily$Date[-1], Mkt = FF_daily$Mkt[-1],
#                             rf = FF_daily$RF[-1], "Mkt-rf" = Mkt - rf)
# 
# for (i in 1:(n_days-1)) {
#   if (i < long) {
#     daily1$Var[i] <- var(daily1$`Mkt-rf`)
#   }
#   else {
#     if (var(daily1$`Mkt-rf`[(i-short):(i-1)]) > var(daily1$`Mkt-rf`[(i-long):(i-1)])) {
#       daily1$Var[i] <- max(var(daily1$`Mkt-rf`[(i-short):(i-1)]), 0.5)
#     }
#     else {
#       daily1$Var[i] <- max(var(daily1$`Mkt-rf`[(i-long):(i-1)]), 0.5)
#     }
#   }
# } 
# 
# a1 <- var(daily1$`Mkt-rf`/daily1$Var)
# b1 <- 2*cov(daily1$`Mkt-rf`/daily1$Var, daily1$rf)
# c1 <- var(daily1$rf)-var(daily1$Mkt)
# 
# c1 <- 1/(2*a1)*(-b1+sqrt((b1)^2-4*a1*c1))
# 
# weights1 <- c(1:(n_days-1))
# returns1 <- c(1:(n_days-1))
# 
# for (month in 1:(n_days-1)) {
#   weights1[month] <- c1 / daily1$Var[month]
#   returns1[month] <- weights1[month] * (daily1$Mkt[month] - daily1$rf[month]) + 
#     daily1$rf[month]
# }
# 
# tot_ret1 <- data.frame(matrix(ncol = 2, nrow = (n_days-1)))
# colnames(tot_ret1) <- c("Mkt", "Strategy")
# 
# tot_ret1[1,] <- 1
# 
# for (month in 2:(n_days-1)) {
#   tot_ret1$Mkt[month] <- tot_ret1$Mkt[month-1] * 
#     (1 + (daily1$Mkt[month-1]/100))
#   tot_ret1$Strategy[month] <- tot_ret1$Strategy[month-1] * 
#     (1 + (returns1[month-1]/100))
# }
# 
# returns1_reg <- (returns1 - daily1$rf) * 252
# market1_reg <- daily1$`Mkt-rf` * 252
# summary(lm(returns1_reg ~ market1_reg))
# summary(lm(returns1_reg ~ market1_reg))$coefficient[1]/ 
#   sigma(lm(returns1_reg ~ market1_reg)) * sqrt(252)




daily <- daily %>% mutate(Vol = sqrt(Var), "Mkt-rf" = Mkt - rf)

a <- var(daily$`Mkt-rf`/daily$Var)
b <- 2*cov(daily$`Mkt-rf`/daily$Var, daily$rf)
c <- var(daily$rf)-var(daily$Mkt)

c <- 1/(2*a)*(-b+sqrt((b)^2-4*a*c))

weights <- c(1:n)
returns <- c(1:n)

for (month in 1:n) {
  weights[month] <- c / daily$Var[month]
  returns[month] <- weights[month] * (daily$Mkt[month] - daily$rf[month]) + 
    daily$rf[month]
}

tot_ret <- data.frame(matrix(ncol = 2, nrow = n))
colnames(tot_ret) <- c("Mkt", "Strategy")

tot_ret[1,] <- 1

for (month in 2:n) {
  tot_ret$Mkt[month] <- tot_ret$Mkt[month-1] * 
    (1 + (daily$Mkt[month-1]/100))
  tot_ret$Strategy[month] <- tot_ret$Strategy[month-1] * 
    (1 + (returns[month-1]/100))
}
tail(tot_ret)
plot(tot_ret$Strategy, type = "l", col = "red")
lines(tot_ret$Mkt, type = "l")
returns1 <- (returns - daily$rf) * 252 / frequ
market <- daily$`Mkt-rf` * 252 / frequ
summary(lm(returns1 ~ market))
summary(lm(returns1 ~ market))$coefficient[1]/ sigma(lm(returns1 ~ market)) * sqrt(252 / frequ)


# ***** Calculate Monthly Variances *****
trading_days <- 22
trading_months <- 12

monthly_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(`Mkt-RF`) * 21)

monthly_vars <- monthly_vars %>% mutate(volatility = sqrt(variance), 
                                        log_var = log(variance), log_vol = log(volatility))
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
monthly_vars$var_of_var_log <- c(1:n_months)
for (i in 1:n_months-1) {
  if (i <= var_of_var_periods) {
    monthly_vars$var_of_var[i] <- 
      var(monthly_vars$variance[1:var_of_var_periods])
    monthly_vars$var_of_var_log[i] <- 
      var(monthly_vars$log_var[1:var_of_var_periods])
  }
  else {
    monthly_vars$var_of_var[i] <- 
      var(monthly_vars$variance[(i-var_of_var_periods + 1):i])
    monthly_vars$var_of_var_log[i] <- 
      var(monthly_vars$log_var[(i-var_of_var_periods + 1):i])
  }
}

# **** 5: Vol of vol
vol_of_vol_periods <- 12

monthly_vars$vol_of_vol <- c(1:n_months)
monthly_vars$vol_of_vol_log <- c(1:n_months)
for (i in 1:n_months) {
  if (i <= vol_of_vol_periods) {
    monthly_vars$vol_of_vol[i] <- sd(monthly_vars$volatility[1:vol_of_vol_periods])
    monthly_vars$vol_of_vol_log[i] <- sd(monthly_vars$log_vol[1:vol_of_vol_periods])
  }
  else {
    monthly_vars$vol_of_vol[i] <- sd(monthly_vars$volatility[(i-vol_of_vol_periods + 1):i])
    monthly_vars$vol_of_vol_log[i] <- sd(monthly_vars$log_vol[(i-vol_of_vol_periods + 1):i])
  }
}

#*** 6: recent emphasis
emphasized_days <- 5

recent_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(tail(`Mkt-RF`, emphasized_days)) * trading_days)

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

vol_of_vol_factor <- c(1:(n_months))
for (i in 1:n_months) {
  if (i <= vol_of_vol_periods) vol_of_vol_factor[i] <- 1
  else {
    if (monthly_vars$vol_of_vol_d[i] <= 
        mean(monthly_vars$vol_of_vol_d[(i - vol_of_vol_periods):(i - 1)])) {
      vol_of_vol_factor[i] <- vol_of_vol_weight
    }
    else {vol_of_vol_factor[i] <- 2 - vol_of_vol_weight}
  }
}

var_factor <- c(1:(n_months))
var_diff <- diff(monthly_vars$variance)
var_factor[i] <- 1
var_periods <- 12
var_weight <- 1.2
for (i in 1:n_months-1) {
  if (i <= var_periods) var_factor[i+1] <- 1
  else {
    if (var_diff[i] <= 
        mean(var_diff[(i - var_periods):(i - 1)])) {
      var_factor[i+1] <- 2 - var_weight
    }
    else {var_factor[i+1] <- var_weight}
  }
}

var_factor2 <- c(1:(n_months))
var_factor2[i] <- 1
for (i in 1:n_months-1) {
  if (i <= var_periods) var_factor[i+1] <- 1
  else {
    if (var_diff[i] >= 
        mean(var_diff[(i - var_periods):(i - 1)]) + sd(var_diff[(i - var_periods):(i - 1)])) {
      var_factor2[i+1] <- 0.8
    }
    else {var_factor2[i+1] <- 1.075}
  }
}

var_factor1 <- c(1:(n_months))
var_periods <- 12
var_weight <- 0.8
for (i in 1:n_months) {
  if (i <= var_periods) var_factor1[i] <- 1
  else {
    if (monthly_vars$variance[i] > 
        (mean(monthly_vars$variance[(i - var_periods):(i - 1)]) + 
        sd(monthly_vars$variance[(i - var_periods):(i - 1)]))) {
      var_factor1[i] <- var_weight
    }
    else if (monthly_vars$variance[i] <
             (mean(monthly_vars$variance[(i - var_periods):(i - 1)]) - 
              sd(monthly_vars$variance[(i - var_periods):(i - 1)]))) {
      var_factor1[i] <- var_weight
    }
    else {var_factor1[i] <- 2 - var_weight}
  }
}
# initiate strategy name vector
names <- c("var_managed", "vol_managed", "ARMA_vol_managed", "EWMA_vol_managed", 
           "var_of_var", "vol_of_vol", "recent_emphasize", "vol_of_vol_d", 
           "Vol_of_vol_inv", "log_var", "log_vol", "log_vol_inv", "var_factor", "test", "test1")

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
denom[,10] <- monthly_vars$variance[-n_months] * 
  monthly_vars$var_of_var_log[-n_months]
denom[,11] <- monthly_vars$volatility[-n_months] * 
  monthly_vars$vol_of_vol_log[-n_months]
denom[,12] <- monthly_vars$variance[-n_months] / 
  monthly_vars$vol_of_vol_log[-n_months]
denom[,13] <- monthly_vars$variance[-n_months] *
  var_factor[-n_months]
denom[,14] <- monthly_vars$variance[-n_months] *
  var_factor1[-n_months]
denom[,15] <- monthly_vars$variance[-n_months] *
  var_factor2[-n_months]

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
  geom_line(aes(y=var_factor)) +
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
