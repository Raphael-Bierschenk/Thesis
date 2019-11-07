# Set Up: Initiate Packages
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

# Import Data
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)
VIX_daily <- read_csv("VIXCLS.csv", col_names = TRUE, 
                      col_types = cols(Date = "D", VIXCLS = "n"))

# Rename and Crop Data
FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)
VIX_daily <- VIX_daily %>% rename(Date = DATE, VIX = VIXCLS)

first_day <- 19260701
first_month <- 192607
last_day <- 20150430
last_month <- 201504

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)
# VIX_daily <- VIX_daily %>% subset(subset = Date <= ymd(last_day) & Date >= ymd(first_day))

# Mutate Data
n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))

FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")

# VIX_daily <- VIX_daily[VIX_daily$Date %in% FF_daily$Date,]
# 
# VIX_monthly <- data.frame(matrix(ncol = ncol(VIX_daily), nrow = n_months))
# colnames(VIX_monthly) <- colnames(VIX_daily)
# j <- 1
# VIX_monthly <- VIX_monthly %>% mutate(Date = FF_monthly$Date)
# 
# for (i in 1:as.integer(count(VIX_daily))) {
#   if (i == 1) {
#     VIX_monthly$VIX[j] <- VIX_daily$VIX[j]
#     j <- j + 1
#   }
#   else if (month(VIX_daily$Date[i]) != month(VIX_daily$Date[i-1])) {
#     VIX_monthly$VIX[j] <- ifelse(!is.na(VIX_daily$VIX[i-1]), 
#                                  VIX_daily$VIX[i-1], VIX_daily$VIX[i])
#     j <- j + 1
#   }
# }

################################################################################
#******************************* Monthly Level *******************************#

# Calculate Monthly Variances
trading_days <- 22
trading_months <- 12

var_m <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(`Mkt-RF`) * 21)

var_m <- var_m %>% mutate(volatility = sqrt(variance), 
                          log_var = log(variance), log_vol = log(volatility))
plot(var_m$volatility * sqrt(trading_months), type = "l")

# Scenario 2: ARIMA Model
variance_ts_m <- xts(var_m$variance, order.by = FF_monthly$Date)
min_obs <- 12
max_obs <- 200
ARMA_model_m <- vector(mode = "list", length = n_months - min_obs)

var_m$ARMA_var <- c(1:n_months)
for (i in 1:n_months) {
  if (i <= min_obs) var_m$ARMA_var[i] <- variance_ts_m[i]
  else {
    ARMA_model_m[[i-min_obs]] <- arima(variance_ts_m[max(1,(i-max_obs+1)):i], 
                                       order = c(1,1,0), method = "ML")
    var_m$ARMA_var[i] <- as.numeric(forecast(ARMA_model_m[[i-min_obs]], h = 1)$mean)
  }
}

# p <- vector(length = n_months - min_obs)
# d <- vector(length = n_months - min_obs)
# q <- vector(length = n_months - min_obs)
# for (i in 1:(n_months - min_obs)) {
#   p[i] <- arimaorder(model[[i]])[1]
#   d[i] <- arimaorder(model[[i]])[2]
#   q[i] <- arimaorder(model[[i]])[3]
# }
# plot(p, type = "l")
# plot(d, type = "l")
# plot(q, type = "l")

# Scenario 3: EWMA Model
n <- 0.94
var_m$EWMA_var <- c(1:n_months)
for (i in 1:n_months) {
  if (i == 1) var_m$EWMA_var[i] <- var_m$variance[i]
  else {
    var_m$EWMA_var[i] <- 
      n * var_m$EWMA_var[i - 1] + (1 - n) * var_m$variance[i]
  }
}

# Scenario 4: GARCH
omega = 0.000001
alpha = 0.074715
beta = 0.9175
var_m$Garch_var <- c(1:n_months)
var_m$Garch_var[1] <- var_m$variance[1]
for (i in 2:n_months) {
  var_m$Garch_var[i] <- omega + alpha*var_m$variance[i]
  + beta*var_m$Garch_var[i-1]
}

# Strategy Names
names <- c("var_managed", "vol_managed", "ARMA_var_managed", "EWMA_var_managed",
           "GARCH_var_managed")

# Set Scale Denominator
denom_m <- data.frame(matrix(ncol = length(names), nrow = n_months - 1))
colnames(denom_m) <- names

denom_m[,1] <- var_m$variance[-n_months]
denom_m[,2] <- var_m$volatility[-n_months]
denom_m[,3] <- var_m$ARMA_var[-n_months]
denom_m[,4] <- var_m$EWMA_var[-n_months]
denom_m[,5] <- var_m$Garch_var[-n_months]

# Calculate c with Midnight Formula
c_m <- data.frame(matrix(ncol = length(names)))
colnames(c_m) <- names

a_qe <- c(1:length(names))
b_qe <- c(1:length(names))
c_qe <- c(1:length(names))

for (i in 1:length(names)) {
  a_qe[i] <- var(FF_monthly$`Mkt-RF`[-1]/denom_m[,i])
  b_qe[i] <- 2*cov(FF_monthly$`Mkt-RF`[-1]/denom_m[,i], FF_monthly$RF[-1])
  c_qe[i] <- var(FF_monthly$RF[-1])-var(FF_monthly$Mkt[-1])
  
  c_m[i] <- 1/(2*a_qe[i])*(-b_qe[i]+sqrt((b_qe[i])^2-4*a_qe[i]*c_qe[i]))
}

# Calculate Weights and Returns
weights_m <- data.frame(matrix(ncol = length(names), nrow = n_months - 1))
colnames(weights_m) <- names

returns_m <- data.frame(matrix(ncol = length(names) + 3, nrow = n_months - 1))
colnames(returns_m) <- c("Date", "Mkt", "rf", names)
returns_m <- returns_m %>% mutate(Date = FF_monthly$Date[-1], 
                                  Mkt = FF_monthly$Mkt[-1], rf = FF_monthly$RF[-1])

for (i in 1:(n_months-1)) {
  for (j in 1:length(names)) {
    weights_m[i,j] <- c_m[j] / denom_m[i,j]
    returns_m[i,names[j]] <- weights_m[i, j] *
      (returns_m$Mkt[i] - returns_m$rf[i]) + returns_m$rf[i]
  }
}

# Check Variance and Display Weight Quantiles
print(apply(returns_m[,-c(1,3)], 2, var))
print(apply(weights_m, 2, quantile, probs = c(0.5, 0.75, 0.9, 0.99)))

# Calculate Total Return
tot_ret_m <- data.frame(matrix(ncol = length(names) + 2, nrow = n_months))
colnames(tot_ret_m) <- c("Date", "Mkt", names)
tot_ret_m <- tot_ret_m %>% mutate(Date = FF_monthly$Date)

tot_ret_m[1,-1] <- 1

for (i in 2:n_months) {
  tot_ret_m$Mkt[i] <- tot_ret_m$Mkt[i-1] * 
    (1 + (returns_m$Mkt[i-1]/100))
  for (j in 1:length(names)) {
    tot_ret_m[i,names[j]] <- tot_ret_m[i-1, names[j]] * 
      (1 + (returns_m[i-1, names[j]]/100))
  }
}

tot_ret_m[n_months,]

# ***** Plot market and VM returns on log scale *****
months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = "month")
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(tot_ret_m, aes(months)) +
  geom_line(aes(y=Mkt)) +
  geom_line(aes(y=var_managed)) +
  geom_line(aes(y=vol_managed)) +
  geom_line(aes(y=ARMA_var_managed)) +
  geom_line(aes(y=EWMA_var_managed)) +
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

# Compute Alpha and Ratios
reg_mkt_m <- vector(mode = "list", length = length(names))
reg_FF3_m <- vector(mode = "list", length = length(names))
b <- 12 * (returns_m$Mkt - returns_m$rf)
b1 <- 12 * (FF_monthly$SMB[-1])
b2 <- 12 * (FF_monthly$HML[-1])

for (i in 1:length(names)) {
  a <- 12 * (returns_m[, names[i]] - returns_m$rf)
  reg_mkt_m[[i]] <- lm(a ~ b)
}

for (i in 1:length(names)) {
  a <- 12 * (returns_m[, names[i]] - returns_m$rf)
  reg_FF3_m[[i]] <- lm(a ~ b + b1 + b2)
}

output_names <- c("alpha_mkt", "R^2_mkt", "RMSE", "SR", "Appr_Ratio", "alpha_FF3")
reg_output <- data.frame(matrix(ncol = length(names), nrow = length(output_names)))
colnames(reg_output) <- names
rownames(reg_output) <- output_names

for (i in 1:length(names)) {
  reg_output["alpha_mkt", i] <- reg_mkt_m[[i]]$coefficients[1]
  reg_output["R^2_mkt", i] <- summary(reg_mkt_m[[i]])$r.squared
  reg_output["RMSE", i] <- sigma(reg_mkt_m[[i]])
  reg_output["SR", i] <- 12 * (mean(returns_m[,names[i]] - returns_m$rf)) / 
    (sqrt(trading_months) * sd(returns_m[,names[i]]))
  reg_output["Appr_Ratio", i] <- sqrt(trading_months) * reg_output["alpha_mkt", i] /
    reg_output["RMSE", i]
  reg_output["alpha_FF3", i] <- reg_FF3_m[[i]]$coefficients[1]
}

round(reg_output, 2)

################################################################################
#******************************** Custom Level ********************************#

# Calculate Custom Variances
frequ <- 5
n_custom <- as.integer(n_days / frequ)
var_c <- data.frame(matrix(ncol = 2, nrow = n_custom))
colnames(var_c) <- c("Date", "variance")

var_c <- var_c %>% mutate(Date = FF_daily$Date[1:n_custom])

for (i in 1:n_days) {
  if(i %% frequ == 0) {
    j <- i / frequ
    var_c$Date[j] <- FF_daily$Date[i]
    var_c$variance[j] <- var(FF_daily$`Mkt-RF`[(i - frequ + 1):i])
  }
}

var_c <- var_c %>% mutate(volatility = sqrt(variance))

# Set Up Return Vector and Calculate Market Return
returns_c <- data.frame(matrix(ncol = 3 + length(names), nrow = n_custom - 1))
colnames(returns_c) <- c("Date", "Mkt", "RF", names)

returns_c <- returns_c %>% mutate(Date = FF_daily$Date[1:(n_custom-1)])

for (i in 1:(n_days - frequ)) {
  if(i %% frequ == 0) {
    j <- i / frequ
    returns_c$Date[j] <- FF_daily$Date[i]
    ret_temp <- 0
    rf_temp <- 0
    for (a in 1:frequ) {
      ret_temp <- ((1 + FF_daily$Mkt[i + a]/100)*(1 + ret_temp/100) - 1)*100
      rf_temp <- ((1 + FF_daily$RF[i + a]/100)*(1 + rf_temp/100) - 1)*100
    }
    returns_c$Mkt[j] <- ret_temp
    returns_c$RF[j] <- rf_temp
  }
}

returns_c <- returns_c %>% mutate("Mkt-RF" = Mkt - RF)

# Scenario 2: ARIMA Model
variance_ts_c <- xts(var_c$variance, order.by = var_c$Date)
min_obs <- 12
max_obs <- 200
ARMA_model_c <- vector(mode = "list", length = n_custom - min_obs)

var_c$ARMA_var <- c(1:n_custom)
for (i in 1:n_custom) {
  if (i <= min_obs) var_c$ARMA_var[i] <- variance_ts_c[i]
  else {
    ARMA_model_c[[i-min_obs]] <- arima(variance_ts_c[max(1,(i-max_obs+1)):i],
                                       order = c(1,0,1), method = "ML")
    var_c$ARMA_var[i] <- as.numeric(forecast(ARMA_model_c[[i-min_obs]], h = 1)$mean)
  }
}

# Scenario 3: EWMA Model
var_c$EWMA_var <- c(1:n_custom)
for (i in 1:n_custom) {
  if (i == 1) var_c$EWMA_var[i] <- var_c$variance[i]
  else {
    var_c$EWMA_var[i] <- 
      n * var_c$EWMA_var[i - 1] + (1 - n) * var_c$variance[i]
  }
}

# Scenario 4: GARCH
var_c$Garch_var <- c(1:n_custom)
var_c$Garch_var[1] <- var_c$variance[1]
for (i in 2:n_custom) {
  var_c$Garch_var[i] <- omega + alpha*var_c$variance[i]
  + beta*var_c$Garch_var[i-1]
}

# Set Scale Denominator
denom_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
colnames(denom_c) <- names

denom_c[,1] <- var_c$variance[-n_custom]
denom_c[,2] <- var_c$volatility[-n_custom]
denom_c[,3] <- var_c$ARMA_var[-n_custom]
denom_c[,4] <- var_c$EWMA_var[-n_custom]
denom_c[,5] <- var_c$Garch_var[-n_custom]

# Calculate c with Midnight Formula
c_c <- data.frame(matrix(ncol = length(names)))
colnames(c_c) <- names

a_qe <- c(1:length(names))
b_qe <- c(1:length(names))
c_qe <- c(1:length(names))

for (i in 1:length(names)) {
  a_qe[i] <- var(returns_c$`Mkt-RF`/denom_c[,i])
  b_qe[i] <- 2*cov(returns_c$`Mkt-RF`/denom_c[,i], returns_c$RF)
  c_qe[i] <- var(returns_c$RF)-var(returns_c$Mkt)
  
  c_c[i] <- 1/(2*a_qe[i])*(-b_qe[i]+sqrt((b_qe[i])^2-4*a_qe[i]*c_qe[i]))
}

# Calculate Weights and Returns
weights_c_th <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
colnames(weights_c_th) <- names

for (i in 1:(n_custom-1)) {
  for (j in 1:length(names)) {
    weights_c_th[i,j] <- (c_c[j] / denom_c[i,j])
  }
}

# Check whether necessary to adjust weight
diff_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 2))
colnames(diff_c) <- names

weights_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
colnames(weights_c) <- names

adjust_n <- 24

for (i in 1:length(names)) {
  diff_c[,names[i]] <- diff(weights_c_th[,names[i]])
}

for (i in 1:(n_custom-1)) {
  for (j in 1:length(names)) {
    if (i <= (adjust_n+1)) {
      weights_c[i,j] <- weights_c_th[i,j]
    }
    else {
      if (diff_c[(i-1),j] > mean(diff_c[(i-1-adjust_n):(i-2),j]) + 
        sd(diff_c[(i-1-adjust_n):(i-2),j]) || diff_c[(i-1),j] < 
        mean(diff_c[(i-1-adjust_n):(i-2),j]) - sd(diff_c[(i-1-adjust_n):(i-2),j])) {
        weights_c[i,j] <- weights_c_th[i,j] # theoretical (correct) weight
      }
      else {
        weights_c[i,j] <- weights_c[(i-1),j] # previous month's weight
      }
    }
  }
}

# Calculate Return --> adjust to weights_c_th here to take out diff
for (i in 1:(n_custom-1)) {
  for (j in 1:length(names)) {
    returns_c[i,names[j]] <- weights_c[i, j] *
      (returns_c$Mkt[i] - returns_c$RF[i]) + returns_c$RF[i]
  }
}

# Check Variance and Display Weight Quantiles
print(apply(returns_c[,-c(1,3,(3+length(names)+1))], 2, var))
print(apply(weights_c, 2, quantile, probs = c(0.5, 0.75, 0.9, 0.99)))

# Calculate Total Return
tot_ret_c <- data.frame(matrix(ncol = length(names) + 2, nrow = n_custom))
colnames(tot_ret_c) <- c("Date", "Mkt", names)
tot_ret_c <- tot_ret_c %>% mutate(Date = var_c$Date)

tot_ret_c[1,-1] <- 1

for (i in 2:n_custom) {
  tot_ret_c$Mkt[i] <- tot_ret_c$Mkt[i-1] * 
    (1 + (returns_c$Mkt[i-1]/100))
  for (j in 1:length(names)) {
    tot_ret_c[i,names[j]] <- tot_ret_c[i-1, names[j]] * 
      (1 + (returns_c[i-1, names[j]]/100))
  }
}

tot_ret_c[n_custom,]

# Compute Alpha and Ratios
reg_mkt_c <- vector(mode = "list", length = length(names))
b <- 252 / frequ * (returns_c$`Mkt-RF`)

for (i in 1:length(names)) {
  a <- 252 / frequ * (returns_c[, names[i]] - returns_c$RF)
  reg_mkt_c[[i]] <- lm(a ~ b)
}

output_names_c <- c("alpha_mkt", "R^2_mkt", "RMSE", "Appr_Ratio")
reg_output_c <- data.frame(matrix(ncol = length(names), nrow = length(output_names_c)))
colnames(reg_output_c) <- names
rownames(reg_output_c) <- output_names_c

for (i in 1:length(names)) {
  reg_output_c["alpha_mkt", i] <- reg_mkt_c[[i]]$coefficients[1]
  reg_output_c["R^2_mkt", i] <- summary(reg_mkt_c[[i]])$r.squared
  reg_output_c["RMSE", i] <- sigma(reg_mkt_c[[i]])
  reg_output_c["Appr_Ratio", i] <- sqrt(252 / frequ) * reg_output_c["alpha_mkt", i] /
    reg_output_c["RMSE", i]
}

round(reg_output_c, 2)

# plot log scale
months <- seq(as.Date("1926/7/1"), as.Date("2015/4/1"), by = frequ)
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(tot_ret_c, aes(tot_ret_c$Date)) +
  geom_line(aes(y=Mkt)) +
  geom_line(aes(y=var_managed)) +
  geom_line(aes(y=vol_managed)) +
  geom_line(aes(y=ARMA_var_managed)) +
  geom_line(aes(y=EWMA_var_managed)) +
  geom_line(aes(y=GARCH_var_managed)) +
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

################################################################################
#******************************** Daily Level ********************************#

# Calculate Daily Variances
var_d <- data.frame(matrix(ncol = 2, nrow = n_days))
colnames(var_d) <- c("Date", "variance")
var_d$variance <- log(1+FF_daily$Mkt/100)^2

var_d <- var_d %>% mutate(volatility = sqrt(variance), Date = FF_daily$Date)

# Scenario 2: ARIMA Model
variance_ts_d <- xts(var_d$variance, order.by = var_d$Date)
min_obs <- 12
max_obs <- 200
ARMA_model_d <- vector(mode = "list", length = n_days - min_obs)

var_d$ARMA_var <- c(1:n_days)
for (i in 1:n_days) {
  if (i <= min_obs) var_d$ARMA_var[i] <- variance_ts_d[i]
  else {
    ARMA_model_d[[i-min_obs]] <- arima(variance_ts_d[max(1,(i-max_obs+1)):i],
                                       order = c(1,0,1), method = "ML")
    var_d$ARMA_var[i] <- as.numeric(forecast(ARMA_model_d[[i-min_obs]], h = 1)$mean)
  }
}

# Scenario 3: EWMA Model
var_d$EWMA_var <- c(1:n_days)
for (i in 1:n_days) {
  if (i == 1) var_d$EWMA_var[i] <- var_d$variance[i]
  else {
    var_d$EWMA_var[i] <- 
      n * var_d$EWMA_var[i - 1] + (1 - n) * var_d$variance[i]
  }
}

# Scenario 4: GARCH
var_d$Garch_var <- c(1:n_days)
var_d$Garch_var[1] <- var_d$variance[1]
for (i in 2:n_days) {
  var_d$Garch_var[i] <- omega + alpha*var_d$variance[i]
  + beta*var_d$Garch_var[i-1]
}

# Set Scale Denominator
names_d <- c("ARMA_var_managed", "EWMA_var_managed", "GARCH_var_managed")

denom_d <- data.frame(matrix(ncol = length(names_d), nrow = n_days - 1))
colnames(denom_d) <- names_d

denom_d[,1] <- var_d$ARMA_var[-n_days]
denom_d[,2] <- var_d$EWMA_var[-n_days]
denom_d[,3] <- var_d$Garch_var[-n_days]

# Calculate c with Midnight Formula
c_d <- data.frame(matrix(ncol = length(names_d)))
colnames(c_d) <- names_d

a_qe <- c(1:length(names_d))
b_qe <- c(1:length(names_d))
c_qe <- c(1:length(names_d))

for (i in 1:length(names_d)) {
  a_qe[i] <- var(FF_daily$`Mkt-RF`[-1]/denom_d[,i])
  b_qe[i] <- 2*cov(FF_daily$`Mkt-RF`[-1]/denom_d[,i], FF_daily$RF[-1])
  c_qe[i] <- var(FF_daily$RF[-1])-var(FF_daily$Mkt[-1])
  
  c_d[i] <- 1/(2*a_qe[i])*(-b_qe[i]+sqrt((b_qe[i])^2-4*a_qe[i]*c_qe[i]))
}

# Calculate Weights and Returns
weights_d <- data.frame(matrix(ncol = length(names_d), nrow = n_days - 1))
colnames(weights_d) <- names_d

returns_d <- data.frame(matrix(ncol = 4 + length(names_d), nrow = n_days - 1))
colnames(returns_d) <- c("Date", "Mkt", "RF", "Mkt-RF", names_d)

returns_d <- returns_d %>% mutate(Date = FF_daily$Date[-1], Mkt = FF_daily$Mkt[-1],
                                  RF = FF_daily$RF[-1], "Mkt-RF" = Mkt - RF)

for (i in 1:(n_days-1)) {
  for (j in 1:length(names_d)) {
    weights_d[i,j] <- c_d[j] / denom_d[i,j]
    returns_d[i,names_d[j]] <- weights_d[i, j] *
      (returns_d$Mkt[i] - returns_d$RF[i]) + returns_d$RF[i]
  }
}

# Check Variance and Display Weight Quantiles
print(apply(returns_d[,-c(1,3)], 2, var))
print(apply(weights_d, 2, quantile, probs = c(0.5, 0.75, 0.9, 0.99)))

# Calculate Total Return
tot_ret_d <- data.frame(matrix(ncol = length(names_d) + 2, nrow = n_days))
colnames(tot_ret_d) <- c("Date", "Mkt", names_d)
tot_ret_d <- tot_ret_d %>% mutate(Date = FF_daily$Date)

tot_ret_d[1,-1] <- 1

for (i in 2:n_days) {
  tot_ret_d$Mkt[i] <- tot_ret_d$Mkt[i-1] * 
    (1 + (returns_d$Mkt[i-1]/100))
  for (j in 1:length(names_d)) {
    tot_ret_d[i,names_d[j]] <- tot_ret_d[i-1, names_d[j]] * 
      (1 + (returns_d[i-1, names_d[j]]/100))
  }
}

tot_ret_d[n_days,]

# Compute Alpha and Ratios
reg_mkt_d <- vector(mode = "list", length = length(names_d))
reg_FF3_d <- vector(mode = "list", length = length(names_d))
b <- 252 * (returns_d$`Mkt-RF`)
b1 <- 252 * (FF_daily$SMB[-1])
b2 <- 252 * (FF_daily$HML[-1])

for (i in 1:length(names_d)) {
  a <- 252 * (returns_d[, names_d[i]] - returns_d$RF)
  reg_mkt_d[[i]] <- lm(a ~ b)
}

for (i in 1:length(names_d)) {
  a <- 252 * (returns_d[, names_d[i]] - returns_d$RF)
  reg_FF3_d[[i]] <- lm(a ~ b + b1 + b2)
}

reg_output_d <- data.frame(matrix(ncol = length(names_d), nrow = length(output_names)))
colnames(reg_output_d) <- names_d
rownames(reg_output_d) <- output_names

for (i in 1:length(names_d)) {
  reg_output_d["alpha_mkt", i] <- reg_mkt_d[[i]]$coefficients[1]
  reg_output_d["R^2_mkt", i] <- summary(reg_mkt_d[[i]])$r.squared
  reg_output_d["RMSE", i] <- sigma(reg_mkt_d[[i]])
  reg_output_d["SR", i] <- 252 * (mean(returns_d[,names_d[i]] - returns_d$RF)) / 
    (sqrt(252) * sd(returns_d[,names_d[i]]))
  reg_output_d["Appr_Ratio", i] <- sqrt(252) * reg_output_d["alpha_mkt", i] /
    reg_output_d["RMSE", i]
  reg_output_d["alpha_FF3", i] <- reg_FF3_d[[i]]$coefficients[1]
}

round(reg_output_d, 2)

# plot log scale
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(tot_ret_d, aes(tot_ret_d$Date)) +
  geom_line(aes(y=Mkt)) +
  geom_line(aes(y=ARMA_var_managed)) +
  geom_line(aes(y=EWMA_var_managed)) +
  geom_line(aes(y=GARCH_var_managed)) +
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