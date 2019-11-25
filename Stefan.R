rm(list = ls())
options(warn=-1)

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

first_day <- 19260701
first_month <- 192607
last_day <- 20190830
last_month <- 201908

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)

n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))


# ***** Manipulate Data *****
FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")
FF_monthly$Date <- as.Date(FF_monthly$Date)

FF_daily$u <- log(1+FF_daily$Mkt/100)
FF_daily$u_sq <- FF_daily$u^2
mean(FF_daily$u) # mean of log changes basically zero as required

# Calculate EWMA Variances
lambda = 0.937456814037826
FF_daily$EWMA_vars <- c(1:nrow(FF_daily))
FF_daily$EWMA_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$EWMA_vars[i] <- lambda*FF_daily$EWMA_vars[i-1]+(1-lambda)*FF_daily$u_sq[i]
}

# Calculate GARCH Variances
omega = 0.00000125947711345891
alpha = 0.0986757938205847
beta = 0.890202868683846
FF_daily$GARCH_vars <- c(1:nrow(FF_daily))
FF_daily$GARCH_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$GARCH_vars[i] <- omega+alpha*FF_daily$u_sq[i]+beta*FF_daily$GARCH_vars[i-1]
}

#########################################################################################################
# Strategy: Flexible time periods determines by fixed deviation of variance, including min and max length
vars_flexible <- data.frame(0, 0, 0, 0, 0)
colnames(vars_flexible) <- c("Index", "Variance", "Volatility", "Mkt", "RF")
min_days_for_var = 5
max_days_for_var = 22
deviation = 0.6
i = 1
last_i = 0
while (i <= nrow(FF_daily)) {
  if (i - last_i < min_days_for_var) {
    i = i + 1
  }
  else {
    if (vars_flexible$Index[1] == 0) {
      ret_temp = 0
      rf_temp = 0
      for (j in c((last_i+1):i)) {
        ret_temp <- ((1 + FF_daily$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
        rf_temp <- ((1 + FF_daily$RF[j]/100)*(1 + rf_temp/100) - 1)*100
      }
      df_temp <- data.frame(i, var(FF_daily$Mkt[(last_i+1):i]), sd(FF_daily$Mkt[(last_i+1):i]), ret_temp, rf_temp)
      vars_flexible[nrow(vars_flexible),] <- df_temp
      last_i = i
    }
    else {
      last_vola = vars_flexible$Volatility[nrow(vars_flexible)]
      if (abs((sd(FF_daily$Mkt[(last_i+1):i]) - last_vola)/last_vola) >= deviation || i - last_i >= max_days_for_var) {
        ret_temp = 0
        rf_temp = 0
        for (j in c((last_i+1):i)) {
          ret_temp <- ((1 + FF_daily$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
          rf_temp <- ((1 + FF_daily$RF[j]/100)*(1 + rf_temp/100) - 1)*100
        }
        df_temp <- data.frame(i, var(FF_daily$Mkt[(last_i+1):i]), sd(FF_daily$Mkt[(last_i+1):i]), ret_temp, rf_temp)
        vars_flexible[nrow(vars_flexible) + 1,] <- df_temp
        last_i = i
      }
    }
    i = i + 1
  }
}
ret_temp = 0
for (j in c((last_i+1):(i-1))) {
  ret_temp <- ((1 + FF_daily$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
  rf_temp <- ((1 + FF_daily$RF[j]/100)*(1 + rf_temp/100) - 1)*100
}
df_temp <- data.frame(i, var(FF_daily$Mkt[(last_i+1):i]), sd(FF_daily$Mkt[(last_i+1):i]), ret_temp, rf_temp)
vars_flexible[nrow(vars_flexible) + 1,] <- df_temp

diff_test <- diff(vars_flexible$Index)
summary(diff_test)
sum(diff_test==7)

vars_flexible$`Mkt-RF` <- vars_flexible$Mkt - vars_flexible$RF

n_periods <- nrow(vars_flexible)
c_flex <- sqrt(var(vars_flexible$`Mkt-RF`[-1]) / var(vars_flexible$`Mkt-RF`[-1] / vars_flexible$Variance[-n_periods]))

weights <- c(1:(n_periods-1))
vola_managed_returns <- c(1:(n_periods-1))

for (period in 2:n_periods) {
  weights[period-1] <- c_flex/vars_flexible$Variance[period-1]
  vola_managed_returns[period-1] <- weights[period-1]*
    (vars_flexible$Mkt[period]-vars_flexible$RF[period])+vars_flexible$RF[period]
}
returns <- data.frame(vars_flexible$Mkt[2:n_periods], vola_managed_returns)
colnames(returns) <- c("Market Returns", "Vola Managed Returns")

plot(vars_flexible$Variance, type = "l")
plot(weights)
print(var(vola_managed_returns))
print(var(vars_flexible$Mkt[2:n_periods]))
print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39

tot_ret = c(1:n_periods)
tot_ret_VM = c(1:n_periods)

for (period in 2:n_periods) {
  tot_ret[period] = tot_ret[period - 1]*(1 + vars_flexible$Mkt[period]/100)
  tot_ret_VM[period] = tot_ret_VM[period - 1]*(1 + vola_managed_returns[period - 1]/100)
}
tot_ret[n_periods]
tot_ret_VM[n_periods]


# ***** Plot market and VM returns on log scale *****
performance <- data.frame(vars_flexible$Index, tot_ret, tot_ret_VM)
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
ggplot(performance, aes(vars_flexible$Index)) +
  geom_line(aes(y=tot_ret)) +
  geom_line(aes(y=tot_ret_VM)) +
  scale_y_continuous(trans = "log10",
                     breaks = trans_breaks('log10', function(x) 10^x),
                     minor_breaks = scale,
                     labels = trans_format('log10', math_format(10^.x)),
                     limits = c(0.1,100000),
                     expand = c(0,0)) +
  ggtitle("Cumulative Performance") + xlab("") + ylab("")

reg_flex <- lm(vola_managed_returns ~ vars_flexible$`Mkt-RF`[-1])
alpha <- reg_flex$coefficients[1]*12*(nrow(vars_flexible)-1)/1065
rmse <- sigma(reg_flex)
SR <- 12 *(nrow(vars_flexible)-1)/1065 * (mean(vola_managed_returns - vars_flexible$RF[-1])) /
  (sqrt(12 *(nrow(vars_flexible)-1)/1065) * sd(vola_managed_returns))
appr <- sqrt(12 *(nrow(vars_flexible)-1)/1065) * alpha / (rmse*22)
alpha
rmse
SR
appr

#######################################################################################################################
# Stefan 2. Strategy
# Calculate EWMA/GARCH Vola daily
# At a certain deviation, reallocate

# Parameters estimated with ML for 19260701 - 20190830
lambda = 0.937456814037826
omega = 0.00000125947711345891
alpha = 0.0986757938205847
beta = 0.890202868683846

output_names <- c("alpha_mkt", "beta_mkt", "R^2_mkt", "RMSE", "SR", "Appr_Ratio", "Performance")
combined_output_flex <- data.frame(row.names = output_names)

intervals <- c(5, 11, 22, 44, 66, 126, 252, 504)

for (i in intervals)
{
  interval = i

  daily_vars <- data.frame(0)
  daily_vars <- data.frame(FF_daily$Date[-c(1:interval)], FF_daily$Mkt[-c(1:interval)], FF_daily$RF[-c(1:interval)])
  colnames(daily_vars) <- c("Date", "Mkt", "RF")
  
  for (i in 1:nrow(daily_vars)) {
    if (interval >= 252) {  # Calculates EWMA and GARCH Vars only for period of interval, otherwise full period
      daily_ewma_var <- c(1:interval)
      daily_garch_var <- c(1:interval)
    
      daily_ewma_var[1] <- FF_daily$u_sq[i]
      daily_garch_var[1] <- FF_daily$u_sq[i]
    
      for (j in 2:interval) {
        daily_ewma_var[j] <- lambda*daily_ewma_var[j-1]+(1-lambda)*FF_daily$u_sq[i+j-1]
        daily_garch_var[j] <- omega+beta*daily_garch_var[j-1]+alpha*FF_daily$u_sq[i+j-1]
      }
      daily_vars$EWMAvars[i] <- daily_ewma_var[interval]
      daily_vars$GARCHvars[i] <- daily_garch_var[interval]
    }
    else {
      daily_vars$EWMAvars[i] <- FF_daily$EWMA_vars[i+interval]
      daily_vars$GARCHvars[i] <- FF_daily$GARCH_vars[i+interval]
    }
    daily_vars$Var[i] <- var(FF_daily$Mkt[i:(i+interval)])
  }
  
  daily_vars$EWMAvars <- daily_vars$EWMAvars*10000
  daily_vars$GARCHvars <- daily_vars$GARCHvars*10000
  ggplot(daily_vars, aes(Date)) + 
    geom_line(aes(y=EWMAvars, colour="red")) + 
    geom_line(aes(y=GARCHvars, colour="green")) + 
    geom_line(aes(y=Var, colour="blue"))
  
  daily_vars$Var_perc_dev[1] <- 0
  daily_vars$EWMA_perc_dev[1] <- 0
  daily_vars$GARCH_perc_dev[1] <- 0
  for (i in 2:nrow(daily_vars)) {
    daily_vars$Var_perc_dev[i] <- daily_vars$Var[i]/daily_vars$Var[i-1] - 1
    daily_vars$EWMA_perc_dev[i] <- daily_vars$EWMAvars[i]/daily_vars$EWMAvars[i-1] - 1
    daily_vars$GARCH_perc_dev[i] <- daily_vars$GARCHvars[i]/daily_vars$GARCHvars[i-1] - 1
  }
  summary(daily_vars$Var_perc_dev)
  summary(daily_vars$EWMA_perc_dev)
  summary(daily_vars$GARCH_perc_dev)
  
  # try to have as many reallocations as months
  quantiles <- apply(daily_vars[,7:9], 2, quantile, probs = c(0.2, 0.8)) #c(1/44, 1-1/44))
  quantiles
  
  # if only top or bottom 10%
  # quantiles[1,] <- c(1000, 1000, 1000) # sufficiently large number as max
  # quantiles[2,] <- apply(daily_vars[,7:9], 2, quantile, probs = 0.95)
  # quantiles
  
  # use the new boundaries
  vars_flexible_v2_Var <- data.frame(ymd("1900/01/01"), 0, 0, 0, 0)
  colnames(vars_flexible_v2_Var) <- c("Date", "Variance", "Volatility", "Mkt", "RF")
  vars_flexible_v2_EWMA <- data.frame(ymd("1900/01/01"), 0, 0, 0, 0)
  colnames(vars_flexible_v2_EWMA) <- c("Date", "Variance", "Volatility", "Mkt", "RF")
  vars_flexible_v2_GARCH <- data.frame(ymd("1900/01/01"), 0, 0, 0, 0)
  colnames(vars_flexible_v2_GARCH) <- c("Date", "Variance", "Volatility", "Mkt", "RF")
  i_Var = 1
  i_EWMA = 1
  i_GARCH = 1
  last_i_Var = 0
  last_i_EWMA = 0
  last_i_GARCH = 0
  while (i_Var < nrow(daily_vars)) {
    # Var
    if (daily_vars$Var_perc_dev[i_Var] > quantiles[2,1] || daily_vars$Var_perc_dev[i_Var] < quantiles[1,1]) { # ggf || and & austauschen
      ret_temp = 0
      rf_temp = 0
      for (j in c((last_i_Var+1):i_Var)) {
        ret_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
        rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
      }
      df_temp <- data.frame(FF_daily$Date[i_Var+interval], daily_vars$Var[i_Var], sqrt(daily_vars$Var[i_Var]), ret_temp, rf_temp)
      vars_flexible_v2_Var[nrow(vars_flexible_v2_Var) + 1,] <- df_temp
      last_i_Var = i_Var
    }
    i_Var = i_Var + 1
  }
  while (i_EWMA < nrow(daily_vars)) {
    # EWMA
    if (daily_vars$EWMA_perc_dev[i_EWMA] > quantiles[2,2] || daily_vars$EWMA_perc_dev[i_EWMA] < quantiles[1,2]) {
      ret_temp = 0
      rf_temp = 0
      for (j in c((last_i_EWMA+1):i_EWMA)) {
        ret_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
        rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
      }
      df_temp <- data.frame(FF_daily$Date[i_EWMA+interval], daily_vars$EWMAvars[i_EWMA], sqrt(daily_vars$EWMAvars[i_EWMA]), ret_temp, rf_temp)
      vars_flexible_v2_EWMA[nrow(vars_flexible_v2_EWMA) + 1,] <- df_temp
      last_i_EWMA = i_EWMA
    }
    i_EWMA = i_EWMA + 1
  }
  while (i_GARCH < nrow(daily_vars)) {
    # GARCH
    if (daily_vars$GARCH_perc_dev[i_GARCH] > quantiles[2,3] || daily_vars$GARCH_perc_dev[i_GARCH] < quantiles[1,3]) {
      ret_temp = 0
      rf_temp = 0
      for (j in c((last_i_GARCH+1):i_GARCH)) {
        ret_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
        rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
      }
      df_temp <- data.frame(FF_daily$Date[i_GARCH+interval], daily_vars$GARCHvars[i_GARCH], sqrt(daily_vars$GARCHvars[i_GARCH]), ret_temp, rf_temp)
      vars_flexible_v2_GARCH[nrow(vars_flexible_v2_GARCH) + 1,] <- df_temp
      last_i_GARCH = i_GARCH
    }
    i_GARCH = i_GARCH + 1
  }
  ret_temp = 0
  for (j in c((last_i_Var+1):i_Var)) {
    ret_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
    rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
  }
  df_temp <- data.frame(FF_daily$Date[i_Var+interval], daily_vars$Var[i_Var], sqrt(daily_vars$Var[i_Var]), ret_temp, rf_temp)
  vars_flexible_v2_Var[nrow(vars_flexible_v2_Var) + 1,] <- df_temp
  ret_temp = 0
  for (j in c((last_i_EWMA+1):i_EWMA)) {
    ret_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
    rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
  }
  df_temp <- data.frame(FF_daily$Date[i_EWMA+interval], daily_vars$EWMAvars[i_EWMA], sqrt(daily_vars$EWMAvars[i_EWMA]), ret_temp, rf_temp)
  vars_flexible_v2_EWMA[nrow(vars_flexible_v2_EWMA) + 1,] <- df_temp
  ret_temp = 0
  for (j in c((last_i_GARCH+1):i_GARCH)) {
    ret_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + ret_temp/100) - 1)*100
    rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
  }
  df_temp <- data.frame(FF_daily$Date[i_GARCH+interval], daily_vars$GARCHvars[i_GARCH], sqrt(daily_vars$GARCHvars[i_GARCH]), ret_temp, rf_temp)
  vars_flexible_v2_GARCH[nrow(vars_flexible_v2_GARCH) + 1,] <- df_temp
  
  vars_flexible_v2_Var$`Mkt-RF` <- vars_flexible_v2_Var$Mkt - vars_flexible_v2_Var$RF
  vars_flexible_v2_EWMA$`Mkt-RF` <- vars_flexible_v2_EWMA$Mkt - vars_flexible_v2_EWMA$RF
  vars_flexible_v2_GARCH$`Mkt-RF` <- vars_flexible_v2_GARCH$Mkt - vars_flexible_v2_GARCH$RF
  
  vars_flexible_v2_Var <- vars_flexible_v2_Var[-1,]
  vars_flexible_v2_EWMA <- vars_flexible_v2_EWMA[-1,]
  vars_flexible_v2_GARCH <- vars_flexible_v2_GARCH[-1,]
  
  n_Var <- nrow(vars_flexible_v2_Var)
  n_EWMA <- nrow(vars_flexible_v2_EWMA)
  n_GARCH <- nrow(vars_flexible_v2_GARCH)
  names_flex <- c("Var", "EWMA", "GARCH")
  
  # Set Scale Denominator
  denom_Var <- vars_flexible_v2_Var$Variance[-n_Var]
  denom_EWMA <- vars_flexible_v2_EWMA$Variance[-n_EWMA]
  denom_GARCH <- vars_flexible_v2_GARCH$Variance[-n_GARCH]
  
  mktrf_Var <- vars_flexible_v2_Var$`Mkt-RF`[-1]
  mktrf_EWMA <- vars_flexible_v2_EWMA$`Mkt-RF`[-1]
  mktrf_GARCH <- vars_flexible_v2_GARCH$`Mkt-RF`[-1]
  
  rf_Var <- vars_flexible_v2_Var$RF[-1]
  rf_EWMA <- vars_flexible_v2_EWMA$RF[-1]
  rf_GARCH <- vars_flexible_v2_GARCH$RF[-1]
  
  # Calculate c with Midnight Formula
  c_var <- data.frame(matrix(ncol = length(names_flex)))
  colnames(c_var) <- names_flex
  
  a_qe <- c(1:length(names_flex))
  b_qe <- c(1:length(names_flex))
  c_qe <- c(1:length(names_flex))
  
  for (i in 1:length(names_flex)) {
    a_qe[i] <- var(get(paste("mktrf_", names_flex[i], sep = ""))/get(paste("denom_", names_flex[i], sep = "")))
    b_qe[i] <- 2*cov(get(paste("mktrf_", names_flex[i], sep = ""))/get(paste("denom_", names_flex[i], sep = "")), 
                     get(paste("rf_", names_flex[i], sep = "")))
    c_qe[i] <- var(get(paste("rf_", names_flex[i], sep = "")))-var(get(paste("mktrf_", names_flex[i], sep = "")) + 
                                                                get(paste("rf_", names_flex[i], sep = "")))
    
    c_var[i] <- 1/(2*a_qe[i])*(-b_qe[i]+sqrt((b_qe[i])^2-4*a_qe[i]*c_qe[i]))
  }
  
  weights_Var <- data.frame(matrix(ncol = 1, nrow = n_Var - 1))
  weights_EWMA <- data.frame(matrix(ncol = 1, nrow = n_EWMA - 1))
  weights_GARCH <- data.frame(matrix(ncol = 1, nrow = n_GARCH - 1))
  
  for (period in 1:(n_Var-1)) {
    weights_Var[period,] <- c_var$Var / denom_Var[period]
  }
  for (period in 1:(n_EWMA-1)) {
    weights_EWMA[period,] <- c_var$EWMA / denom_EWMA[period]
  }
  for (period in 1:(n_GARCH-1)) {
    weights_GARCH[period,] <- c_var$GARCH / denom_GARCH[period]
  }
  quantile(weights_Var[,1], probs = c(0.5, 0.75, 0.9, 0.99))
  quantile(weights_EWMA[,1], probs = c(0.5, 0.75, 0.9, 0.99))
  quantile(weights_GARCH[,1], probs = c(0.5, 0.75, 0.9, 0.99))
  
  # Calculate volatility managed returns (incl. various transaction costs)
  vars_flexible_v2_Var$VMR <- 0
  vars_flexible_v2_Var$VMR_1bps <- 0
  vars_flexible_v2_Var$VMR_10bps <- 0
  vars_flexible_v2_Var$VMR_14bps <- 0
  vars_flexible_v2_EWMA$VMR <- 0
  vars_flexible_v2_EWMA$VMR_1bps <- 0
  vars_flexible_v2_EWMA$VMR_10bps <- 0
  vars_flexible_v2_EWMA$VMR_14bps <- 0
  vars_flexible_v2_GARCH$VMR <- 0
  vars_flexible_v2_GARCH$VMR_1bps <- 0
  vars_flexible_v2_GARCH$VMR_10bps <- 0
  vars_flexible_v2_GARCH$VMR_14bps <- 0
  for (i in 2:n_Var) {
    vars_flexible_v2_Var$VMR[i] <- weights_Var[i-1,1]*vars_flexible_v2_Var$`Mkt-RF`[i]+vars_flexible_v2_Var$RF[i]
    if (i == 2) {
      vars_flexible_v2_Var$VMR_1bps[i] <- weights_Var[i-1,1]*(vars_flexible_v2_Var$`Mkt-RF`[i] - 0.01)+vars_flexible_v2_Var$RF[i]
      vars_flexible_v2_Var$VMR_10bps[i] <- weights_Var[i-1,1]*(vars_flexible_v2_Var$`Mkt-RF`[i] - 0.10)+vars_flexible_v2_Var$RF[i]
      vars_flexible_v2_Var$VMR_14bps[i] <- weights_Var[i-1,1]*(vars_flexible_v2_Var$`Mkt-RF`[i] - 0.14)+vars_flexible_v2_Var$RF[i]
    }
    else {
      vars_flexible_v2_Var$VMR_1bps[i] <- weights_Var[i-1,1]*vars_flexible_v2_Var$`Mkt-RF`[i]+vars_flexible_v2_Var$RF[i] -
        abs(weights_Var[i-1,1]-weights_Var[i-2,1]) * 0.01
      vars_flexible_v2_Var$VMR_10bps[i] <- weights_Var[i-1,1]*vars_flexible_v2_Var$`Mkt-RF`[i]+vars_flexible_v2_Var$RF[i] -
        abs(weights_Var[i-1,1]-weights_Var[i-2,1]) * 0.10
      vars_flexible_v2_Var$VMR_14bps[i] <- weights_Var[i-1,1]*vars_flexible_v2_Var$`Mkt-RF`[i]+vars_flexible_v2_Var$RF[i] -
        abs(weights_Var[i-1,1]-weights_Var[i-2,1]) * 0.14
    }
  }
  for (i in 2:n_EWMA) {
    vars_flexible_v2_EWMA$VMR[i] <- weights_EWMA[i-1,1]*vars_flexible_v2_EWMA$`Mkt-RF`[i]+vars_flexible_v2_EWMA$RF[i]
    if (i == 2) {
      vars_flexible_v2_EWMA$VMR_1bps[i] <- weights_EWMA[i-1,1]*(vars_flexible_v2_EWMA$`Mkt-RF`[i] - 0.01)+vars_flexible_v2_EWMA$RF[i]
      vars_flexible_v2_EWMA$VMR_10bps[i] <- weights_EWMA[i-1,1]*(vars_flexible_v2_EWMA$`Mkt-RF`[i] - 0.10)+vars_flexible_v2_EWMA$RF[i]
      vars_flexible_v2_EWMA$VMR_14bps[i] <- weights_EWMA[i-1,1]*(vars_flexible_v2_EWMA$`Mkt-RF`[i] - 0.14)+vars_flexible_v2_EWMA$RF[i]
    }
    else {
      vars_flexible_v2_EWMA$VMR_1bps[i] <- weights_EWMA[i-1,1]*vars_flexible_v2_EWMA$`Mkt-RF`[i]+vars_flexible_v2_EWMA$RF[i] -
        abs(weights_EWMA[i-1,1]-weights_EWMA[i-2,1]) * 0.01
      vars_flexible_v2_EWMA$VMR_10bps[i] <- weights_EWMA[i-1,1]*vars_flexible_v2_EWMA$`Mkt-RF`[i]+vars_flexible_v2_EWMA$RF[i] -
        abs(weights_EWMA[i-1,1]-weights_EWMA[i-2,1]) * 0.10
      vars_flexible_v2_EWMA$VMR_14bps[i] <- weights_EWMA[i-1,1]*vars_flexible_v2_EWMA$`Mkt-RF`[i]+vars_flexible_v2_EWMA$RF[i] -
        abs(weights_EWMA[i-1,1]-weights_EWMA[i-2,1]) * 0.14
    }
  }
  for (i in 2:n_GARCH) {
    vars_flexible_v2_GARCH$VMR[i] <- weights_GARCH[i-1,1]*vars_flexible_v2_GARCH$`Mkt-RF`[i]+vars_flexible_v2_GARCH$RF[i]
    if (i == 2) {
      vars_flexible_v2_GARCH$VMR_1bps[i] <- weights_GARCH[i-1,1]*(vars_flexible_v2_GARCH$`Mkt-RF`[i] - 0.01)+vars_flexible_v2_GARCH$RF[i]
      vars_flexible_v2_GARCH$VMR_10bps[i] <- weights_GARCH[i-1,1]*(vars_flexible_v2_GARCH$`Mkt-RF`[i] - 0.10)+vars_flexible_v2_GARCH$RF[i]
      vars_flexible_v2_GARCH$VMR_14bps[i] <- weights_GARCH[i-1,1]*(vars_flexible_v2_GARCH$`Mkt-RF`[i] - 0.14)+vars_flexible_v2_GARCH$RF[i]
    }
    else {
      vars_flexible_v2_GARCH$VMR_1bps[i] <- weights_GARCH[i-1,1]*vars_flexible_v2_GARCH$`Mkt-RF`[i]+vars_flexible_v2_GARCH$RF[i] -
        abs(weights_GARCH[i-1,1]-weights_GARCH[i-2,1]) * 0.01
      vars_flexible_v2_GARCH$VMR_10bps[i] <- weights_GARCH[i-1,1]*vars_flexible_v2_GARCH$`Mkt-RF`[i]+vars_flexible_v2_GARCH$RF[i] -
        abs(weights_GARCH[i-1,1]-weights_GARCH[i-2,1]) * 0.10
      vars_flexible_v2_GARCH$VMR_14bps[i] <- weights_GARCH[i-1,1]*vars_flexible_v2_GARCH$`Mkt-RF`[i]+vars_flexible_v2_GARCH$RF[i] -
        abs(weights_GARCH[i-1,1]-weights_GARCH[i-2,1]) * 0.14
    }
  }
  
  ggplot() + 
    geom_line(aes(y = weights_Var[,1], x = c(1:nrow(weights_Var)), colour = "red")) + 
    geom_line(aes(y = weights_EWMA[,1], x = c(1:nrow(weights_EWMA)), colour = "green")) + 
    geom_line(aes(y = weights_GARCH[,1], x = c(1:nrow(weights_GARCH)), colour = "blue")) +
    ggtitle("Weights") + xlab("") + ylab("")
  
  apply(vars_flexible_v2_Var[-1,c(4,7)], 2, var)
  apply(vars_flexible_v2_EWMA[-1,c(4,7)], 2, var)
  apply(vars_flexible_v2_GARCH[-1,c(4,7)], 2, var)
  
  
  # Calculate total returns (incl. various transaction costs)
  tot_ret_VM_Var = data.frame(c(as.Date("1926/07/01"),vars_flexible_v2_Var$Date), 
                              c(1:(n_Var+1)), c(1:(n_Var+1)), c(1:(n_Var+1)), c(1:(n_Var+1)))
  colnames(tot_ret_VM_Var) <- c("Date", "Performance Var", "Performance Var 1bps", "Performance Var 10bps", "Performance Var 14bps")
  tot_ret_VM_EWMA = data.frame(c(as.Date("1926/07/01"),vars_flexible_v2_EWMA$Date), 
                               c(1:(n_EWMA+1)), c(1:(n_EWMA+1)), c(1:(n_EWMA+1)), c(1:(n_EWMA+1)))
  colnames(tot_ret_VM_EWMA) <- c("Date", "Performance EWMA", "Performance EWMA 1bps", "Performance EWMA 10bps", "Performance EWMA 14bps")
  tot_ret_VM_GARCH = data.frame(c(as.Date("1926/07/01"),vars_flexible_v2_GARCH$Date), 
                                c(1:(n_GARCH+1)), c(1:(n_GARCH+1)), c(1:(n_GARCH+1)), c(1:(n_GARCH+1)))
  colnames(tot_ret_VM_GARCH) <- c("Date", "Performance GARCH", "Performance GARCH 1bps", "Performance GARCH 10bps", "Performance GARCH 14bps")
  
  for (period in 2:(n_Var+1)) {
    tot_ret_VM_Var[period, 2] = tot_ret_VM_Var[period - 1, 2]*(1 + vars_flexible_v2_Var$VMR[period - 1]/100)
    tot_ret_VM_Var[period, 3] = tot_ret_VM_Var[period - 1, 3]*(1 + vars_flexible_v2_Var$VMR_1bps[period - 1]/100)
    tot_ret_VM_Var[period, 4] = tot_ret_VM_Var[period - 1, 4]*(1 + vars_flexible_v2_Var$VMR_10bps[period - 1]/100)
    tot_ret_VM_Var[period, 5] = tot_ret_VM_Var[period - 1, 5]*(1 + vars_flexible_v2_Var$VMR_14bps[period - 1]/100)
  }
  for (period in 2:(n_EWMA+1)) {
    tot_ret_VM_EWMA[period, 2] = tot_ret_VM_EWMA[period - 1, 2]*(1 + vars_flexible_v2_EWMA$VMR[period - 1]/100)
    tot_ret_VM_EWMA[period, 3] = tot_ret_VM_EWMA[period - 1, 3]*(1 + vars_flexible_v2_EWMA$VMR_1bps[period - 1]/100)
    tot_ret_VM_EWMA[period, 4] = tot_ret_VM_EWMA[period - 1, 4]*(1 + vars_flexible_v2_EWMA$VMR_10bps[period - 1]/100)
    tot_ret_VM_EWMA[period, 5] = tot_ret_VM_EWMA[period - 1, 5]*(1 + vars_flexible_v2_EWMA$VMR_14bps[period - 1]/100)
  }
  for (period in 2:(n_GARCH+1)) {
    tot_ret_VM_GARCH[period, 2] = tot_ret_VM_GARCH[period - 1, 2]*(1 + vars_flexible_v2_GARCH$VMR[period - 1]/100)
    tot_ret_VM_GARCH[period, 3] = tot_ret_VM_GARCH[period - 1, 3]*(1 + vars_flexible_v2_GARCH$VMR_1bps[period - 1]/100)
    tot_ret_VM_GARCH[period, 4] = tot_ret_VM_GARCH[period - 1, 4]*(1 + vars_flexible_v2_GARCH$VMR_10bps[period - 1]/100)
    tot_ret_VM_GARCH[period, 5] = tot_ret_VM_GARCH[period - 1, 5]*(1 + vars_flexible_v2_GARCH$VMR_14bps[period - 1]/100)
  }
  
  # Market return
  market_returns <- data.frame(FF_monthly$Date, c(1:nrow(FF_monthly)))
  colnames(market_returns) <- c("Date", "Performance Market")
  for (i in 2:nrow(FF_monthly)) {
    market_returns[i,2] <- market_returns[i-1,2]*(1 + FF_monthly$Mkt[i]/100)
  }
  
  
  # ***** Plot market and VM returns on log scale *****
  scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
             200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
             20000,30000,40000,50000,60000,70000,80000,90000,100000)
  dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
  ggplot() +
    geom_line(aes(y=market_returns$`Performance Market`, x=market_returns$Date, colour="Buy and Hold")) +
    geom_line(aes(y=tot_ret_VM_Var$`Performance Var`, x=tot_ret_VM_Var$Date, colour="Realized Variance")) +
    geom_line(aes(y=tot_ret_VM_EWMA$`Performance EWMA`, x=tot_ret_VM_EWMA$Date, colour="EWMA")) +
    geom_line(aes(y=tot_ret_VM_GARCH$`Performance GARCH`, x=tot_ret_VM_GARCH$Date, colour="GARCH")) +
    scale_x_date(limits = as.Date(c("1926-7-1", "2019-8-1")),
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
    theme(legend.position="bottom") +
    ggtitle(paste("Cumulative Performance",interval)) + xlab("") + ylab("") +
    scale_color_manual(name = "Strategies", values = c("Buy and Hold" = "black", "Realized Variance" = "red",
                                                       "EWMA" = "green", "GARCH" = "blue"))
  
  reg_flex_Var <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_Var)
  reg_flex_Var_1bps <- lm(VMR_1bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_Var)
  reg_flex_Var_10bps <- lm(VMR_10bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_Var)
  reg_flex_Var_14bps <- lm(VMR_14bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_Var)
  reg_flex_EWMA <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_EWMA)
  reg_flex_EWMA_1bps <- lm(VMR_1bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_EWMA)
  reg_flex_EWMA_10bps <- lm(VMR_10bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_EWMA)
  reg_flex_EWMA_14bps <- lm(VMR_14bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_EWMA)
  reg_flex_GARCH <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_GARCH)
  reg_flex_GARCH_1bps <- lm(VMR_1bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_GARCH)
  reg_flex_GARCH_10bps <- lm(VMR_10bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_GARCH)
  reg_flex_GARCH_14bps <- lm(VMR_14bps[-1] - RF[-1] ~ `Mkt-RF`[-1], vars_flexible_v2_GARCH)
  
  names_flex_cost <- c("Var", "Var_1bps", "Var_10bps", "Var_14bps", "EWMA", "EWMA_1bps", "EWMA_10bps", "EWMA_14bps",
                       "GARCH", "GARCH_1bps", "GARCH_10bps", "GARCH_14bps")
  reg_output_flex <- data.frame(matrix(ncol = length(names_flex_cost), nrow = length(output_names)))
  colnames(reg_output_flex) <- names_flex_cost
  rownames(reg_output_flex) <- output_names
  
  trading_days <- 252
  factors <- c(1:12)
  factors[1:4] <- nrow(vars_flexible_v2_Var) / (nrow(FF_daily)-interval) * trading_days
  factors[5:8] <- nrow(vars_flexible_v2_EWMA) / (nrow(FF_daily)-interval) * trading_days
  factors[9:12] <- nrow(vars_flexible_v2_GARCH) / (nrow(FF_daily)-interval) * trading_days
  
  for (i in 1:length(names_flex_cost)) {
    reg_output_flex["alpha_mkt", i] <- get(paste("reg_flex_", names_flex_cost[i], sep = ""))$coefficients[1] * factors[i]
    reg_output_flex["beta_mkt", i] <- get(paste("reg_flex_", names_flex_cost[i], sep = ""))$coefficients[2]
    reg_output_flex["R^2_mkt", i] <- summary(get(paste("reg_flex_", names_flex_cost[i], sep = "")))$r.squared
    reg_output_flex["RMSE", i] <- sigma(get(paste("reg_flex_", names_flex_cost[i], sep = ""))) * factors[i]
    reg_output_flex["SR", i] <- 1 #factors[i] * 
      #mean(get(paste("vars_flexible_v2_", names_flex_cost[i], sep = ""))$VMR - 
      #        get(paste("vars_flexible_v2_", names_flex_cost[i], sep = ""))$RF) / 
      #(sqrt(factors[i]) * sd(get(paste("vars_flexible_v2_", names_flex_cost[i], sep = ""))$VMR))
    reg_output_flex["Appr_Ratio", i] <- sqrt(factors[i]) * reg_output_flex["alpha_mkt", i] / reg_output_flex["RMSE", i]
    reg_output_flex["Performance", i] <- tail(get(paste("tot_ret_VM_", names_flex[ceiling(i/4)], sep = ""))[,c(5,2,3,4)[i%%4+1]], n = 1)
  }

  # round(reg_output_flex, 2)
  
  for (i in 1:ncol(reg_output_flex)) {
    combined_output_flex[,ncol(combined_output_flex)+1] <- round(reg_output_flex[,i],2)
    col = ncol(combined_output_flex)
    colnames(combined_output_flex)[col] <- paste(colnames(reg_output_flex)[i], interval, sep = "_")
  }
}
combined_output_flex

# alpha and apprasial ratio table (costs <-> interval)
output_flex_alpha <- data.frame(matrix(ncol = length(intervals), nrow = length(names_flex_cost)))
rownames(output_flex_alpha) <- names_flex_cost
colnames(output_flex_alpha) <- c("1 week", "1/2 months", "1 months", "2 months", "3 months", "6 months", "1 year", "2 years")
output_flex_appr_ratio <- data.frame(matrix(ncol = length(intervals), nrow = length(names_flex_cost)))
rownames(output_flex_appr_ratio) <- names_flex_cost
colnames(output_flex_appr_ratio) <- c("1 week", "1/2 months", "1 months", "2 months", "3 months", "6 months", "1 year", "2 years")
for (i in 1:length(intervals)) {
  for (j in 1:length(names_flex_cost)) {
    output_flex_alpha[j,i] <- combined_output_flex[1,j+(i-1)*12]
    output_flex_appr_ratio[j,i] <- combined_output_flex[6,j+(i-1)*12]
  }
}
# Alphas
output_flex_alpha
# Apprasial Ratios
output_flex_appr_ratio


# improve table (without costs)
output_flex_var <- data.frame(matrix(ncol = length(intervals), nrow = length(output_names)))
rownames(output_flex_var) <- output_names
colnames(output_flex_var) <- c("1 week", "1/2 months", "1 months", "2 months", "3 months", "6 months", "1 year", "2 years")
for (i in 1:length(intervals)) {
  output_flex_var[,i] <- combined_output_flex[,i*3-2]
}

output_flex_ewma <- data.frame(matrix(ncol = 3, nrow = length(output_names)))
rownames(output_flex_ewma) <- output_names
colnames(output_flex_ewma) <- c("From start", "1 year", "2 years")
output_flex_ewma[,1] <- combined_output_flex[,8]
output_flex_ewma[,2] <- combined_output_flex[,20]
output_flex_ewma[,3] <- combined_output_flex[,23]

output_flex_garch <- data.frame(matrix(ncol = 3, nrow = length(output_names)))
rownames(output_flex_garch) <- output_names
colnames(output_flex_garch) <- c("From start", "1 year", "2 years")
output_flex_garch[,1] <- combined_output_flex[,9]
output_flex_garch[,2] <- combined_output_flex[,21]
output_flex_garch[,3] <- combined_output_flex[,24]

output_flex_var
output_flex_ewma
output_flex_garch


filter(FF_monthly, Mkt < -10)
filter(vars_flexible_v2_Var, VMR < -10)
filter(vars_flexible_v2_EWMA, VMR < -10)
filter(vars_flexible_v2_GARCH, VMR < -10)

#********************************************************************************************************




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
var_diff <- diff(monthly_vars$variance)
var_factor2[i] <- 1
var_periods <- 12
var_weight <- 1.2
for (i in 1:n_months-1) {
  if (i <= var_periods) var_factor2[i+1] <- 1
  else {
    if (abs(var_diff[i]) <= 
        mean(abs(var_diff[(i - var_periods):(i - 1)]))) {
      var_factor2[i+1] <- monthly_vars$variance[i+1]
    }
    else {var_factor2[i+1] <- monthly_vars$recent_var[i+1]}
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
denom[,15] <- var_factor2[-n_months]

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

ggplot() +
  geom_line(aes(y=diff(monthly_vars$volatility[-(1:499)]), x=months[-(1:500)]), color = "red") +
  geom_line(aes(y=diff(FF_monthly$`Mkt-RF`[-(1:499)]), x=months[-(1:500)])) +
  geom_line(aes(y=monthly_vars$vol_of_vol_d[-(1:500)], x=months[-(1:500)]), color = "blue")
summary(lm(FF_monthly$`Mkt-RF`[-1] ~ monthly_vars$volatility[-n_months] ))
cov(diff(monthly_vars$variance[-n_months]), FF_monthly$`Mkt-RF`[-(1:2)])
summary(lm(FF_monthly$`Mkt-RF`[-(1:2)] ~ diff(monthly_vars$variance[-n_months])))

