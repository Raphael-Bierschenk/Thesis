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
library(broom)
library(ggthemes)
library(moments)
library(optimx)
library(stargazer)
library(rugarch)

# Import Data
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)
recession <- read_csv("USREC.csv", col_names = TRUE)

# Rename and Crop Data
FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)
recession <- recession %>% rename(Date = DATE, indicator = USREC)

first_day <- 19260701
first_month <- 192607
last_day <- 20190830
last_month <- 201908

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)
recession <- recession %>% subset(subset = Date <= ymd(last_day) & Date >= ymd(first_day))

# Mutate Data
n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))

FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")
FF_monthly$Date <- ymd(FF_monthly$Date)

FF_daily$u_sq <- log(1+FF_daily$`Mkt-RF`/100)^2
FF_monthly$u_sq <- log(1+FF_monthly$`Mkt-RF`/100)^2

################################################################################
#******************************* Monthly Level *******************************#

# Calculate Monthly Variances
trading_days <- 22
trading_months <- 12

var_m <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(`Mkt-RF`) * trading_days)

var_m <- as.data.frame(var_m)
var_m$Date <- FF_monthly$Date
var_m <- var_m %>% mutate(volatility = sqrt(variance))
plot(var_m$volatility * sqrt(trading_months), type = "l")

# Scenario 2: ARIMA Model
variance_ts_m <- xts(var_m$variance, order.by = FF_monthly$Date)

kpss.test(var_m$variance)
kpss.test(diff(var_m$variance))
adf.test(diff(var_m$variance))
pp.test(diff(var_m$variance))

Acf(diff(variance_ts_m), lag = 24, main = "ACF Differenced Monthly Variance")
Pacf(diff(variance_ts_m), lag = 24, main = "PACF Differenced Monthly Variance")

ARMA_model_m <- auto.arima(variance_ts_m,stepwise = FALSE, approximation = FALSE)
jarque.test(as.vector(ARMA_model_m$residuals))

var_m$ARMA_var <- c(fitted(ARMA_model_m)[-1], forecast(ARMA_model_m, h = 1)$mean)
var_m$ARMA_vol <- sqrt(var_m$ARMA_var)

# Scenario 3: EWMA Model
# Optimize Lambda
EWMA_function_m <- function(lambda)
{
  EWMA_var <- c(1:n_months)
  EWMA_var[1] <- FF_monthly$u_sq[1]
  for(i in 2:n_months) {
    EWMA_var[i] <- lambda * EWMA_var[i-1] + (1 - lambda) * FF_monthly$u_sq[i]
  }
  EWMA_likelihood <- c(1:(n_months-1))
  for(i in 1:(n_months-1)) {
    EWMA_likelihood[i] <- -log(EWMA_var[i])-FF_monthly$u_sq[i+1]/EWMA_var[i]
  }
  return (sum(EWMA_likelihood))
}
ewma_max_m <- optimize(EWMA_function_m, interval = c(0, 1), maximum = TRUE, 
                       tol = 0.000000000000001)
lambda_m <- ewma_max_m$maximum

# Calculate EWMA Variance
var_m$EWMA_var <- c(1:n_months)
var_m$EWMA_var[1] <- var_m$variance[1]
for (i in 2:n_months) {
  var_m$EWMA_var[i] <- lambda_m * var_m$EWMA_var[i - 1] + (1 - lambda_m) * 
    var_m$variance[i]
  }
var_m$EWMA_vol <- sqrt(var_m$EWMA_var)

# Scenario 4: GARCH
# Optimize Parameter
GARCH_model <- ugarchfit(ugarchspec(
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)), data = sqrt(FF_monthly$u_sq),
  solver = "nloptr", solver.control = list(solver = 4))

GARCH_function_m <- function(alpha, beta) {
  omega <- max(0,mean(FF_monthly$u_sq)*(1-alpha-beta))
  GARCH_var <- c(1:n_months)
  GARCH_var[1] <- FF_monthly$u_sq[1]
  for(i in 2:n_months) {
    GARCH_var[i] <- omega + beta*GARCH_var[i-1] + alpha*FF_monthly$u_sq[i]
  }
  GARCH_likelihood <- c(1:(n_months - 1))
  for(i in 1:(n_months-1)) {
    GARCH_likelihood[i] <- -log(GARCH_var[i])-FF_monthly$u_sq[i+1]/GARCH_var[i]
  }
  return (sum(GARCH_likelihood))
}

GARCH_max_m <- optimx(c(0.1, 0.9), function(x) GARCH_function_m(x[1], x[2]), 
                      method = "Nelder-Mead", control = list(maximize = TRUE))
alpha_m <- GARCH_max_m$p1
beta_m <- GARCH_max_m$p2
omega_m <- max(0,mean(FF_monthly$u_sq)*(1-alpha_m-beta_m))

# Calculate GARCH Variance
var_m$GARCH_var <- c(1:n_months)
var_m$GARCH_var[1] <- FF_monthly$u_sq[1] * 10000
for (i in 2:n_months) {
  var_m$GARCH_var[i] <- (omega_m +
    alpha_m*FF_monthly$u_sq[i] + beta_m * var_m$GARCH_var[i-1] / 10000) * 10000
}
var_m$GARCH_vol <- sqrt(var_m$GARCH_var)

# Analyze Variance Estimates
ggplot(var_m, aes(x = Date)) +
  geom_line(aes(y = volatility * sqrt(12), color = "Realized Variance")) +
  geom_line(aes(y = ARMA_vol * sqrt(12), color = "ARMA")) +
  geom_line(aes(y = EWMA_vol * sqrt(12), color = "EWMA")) +
  geom_line(aes(y = GARCH_vol * sqrt(12), color = "GARCH")) +
  theme_minimal() +
  ggtitle("Annualized Volatility") +
  ylab("") +
  scale_color_manual(name = "", values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

var_names <- c("variance", "ARMA_var", "EWMA_var", "GARCH_var")
var_var_reg_m <- vector(mode = "list", length = length(names))
ret_var_reg_m <- vector(mode = "list", length = length(names))
retvar_var_reg_m <- vector(mode = "list", length = length(names))

for (i in 1:length(var_names)) {
  var_var_reg_m[[i]] <- lm(var_m$variance[-1] ~ var_m[-n_months, var_names[i]])
  ret_var_reg_m[[i]] <- lm(FF_monthly$`Mkt-RF`[-1] ~ var_m[-n_months, var_names[i]])
  retvar_var_reg_m[[i]] <- lm((FF_monthly$`Mkt-RF`[-1] / var_m$variance[-1]) ~
                                var_m[-n_months, var_names[i]])
  names(var_var_reg_m[[i]]$coefficients) <- 
    c(names(var_var_reg_m[[i]]$coefficients[1]), var_names[i])
  names(ret_var_reg_m[[i]]$coefficients) <- 
    c(names(var_var_reg_m[[i]]$coefficients[1]), var_names[i])
  names(retvar_var_reg_m[[i]]$coefficients) <- 
    c(names(var_var_reg_m[[i]]$coefficients[1]), var_names[i])
}

stargazer(var_var_reg_m[[1]], var_var_reg_m[[2]], var_var_reg_m[[3]], var_var_reg_m[[4]],
          ret_var_reg_m[[1]], ret_var_reg_m[[2]], ret_var_reg_m[[3]], ret_var_reg_m[[4]],
          retvar_var_reg_m[[1]], retvar_var_reg_m[[2]], retvar_var_reg_m[[3]], retvar_var_reg_m[[4]],
          type = "text", dep.var.labels = c("Variance Next Month","Return Next Month",
                                            "Return per Unit of Variance"), 
          covariate.labels = c("Realized Variance", "ARMA", "EWMA", "GARCH"),
          out = "table.htm")
stargazer(ret_var_reg_m[[1]], ret_var_reg_m[[2]], ret_var_reg_m[[3]], ret_var_reg_m[[4]])
stargazer(retvar_var_reg_m[[1]], retvar_var_reg_m[[2]], retvar_var_reg_m[[3]], retvar_var_reg_m[[4]])

# Strategy Names
names <- c("var_managed","ARMA_var_managed", "EWMA_var_managed",
           "GARCH_var_managed")

# Set Scale Denominator
denom_m <- data.frame(matrix(ncol = length(names), nrow = n_months - 1))
colnames(denom_m) <- names

denom_m[,1] <- var_m$variance[-n_months]
denom_m[,2] <- var_m$ARMA_var[-n_months]
denom_m[,3] <- var_m$EWMA_var[-n_months]
denom_m[,4] <- var_m$GARCH_var[-n_months]

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
trading_cost <- 0.1
w_abs_m <- data.frame(matrix(ncol = length(names), nrow = n_months - 1))
colnames(w_abs_m) <- names
w_abs_m[1,] <- 1

weights_m <- data.frame(matrix(ncol = length(names) + 1, nrow = n_months - 1))
colnames(weights_m) <- c("Date", names)
weights_m <- weights_m %>% mutate(Date = FF_monthly$Date[-1])

returns_m <- data.frame(matrix(ncol = length(names) + 3, nrow = n_months - 1))
colnames(returns_m) <- c("Date", "Mkt", "rf", names)
returns_m <- returns_m %>% mutate(Date = FF_monthly$Date[-1], 
                                  Mkt = FF_monthly$Mkt[-1], rf = FF_monthly$RF[-1])

returns_cost_m <- data.frame(matrix(ncol = length(names) + 1, nrow = n_months - 1))
colnames(returns_cost_m) <- c("Date", names)
returns_cost_m <- returns_cost_m %>% mutate(Date = FF_monthly$Date[-1])

for (i in 1:(n_months-1)) {
  for (j in 1:length(names)) {
    weights_m[i,names[j]] <- c_m[j] / denom_m[i,j]
    if (i != 1) {
      w_abs_m[i,names[j]] <- abs(weights_m[i,names[j]] - weights_m[i-1,names[j]])
    }
    returns_m[i,names[j]] <- weights_m[i, names[j]] *
      (returns_m$Mkt[i] - returns_m$rf[i]) + returns_m$rf[i]
    returns_cost_m[i,names[j]] <- returns_m[i,names[j]] - 
      w_abs_m[i,names[j]] * trading_cost
  }
}

# Check Variance
print(apply(returns_m[,-c(1,3)], 2, var))
quantile_m <- as.data.frame(apply(weights_m[,-1], 2, quantile, 
                                  probs = c(0.5, 0.75, 0.9, 0.99)))

# Calculate Total Return
tot_ret_m <- data.frame(matrix(ncol = length(names) + 2, nrow = n_months))
colnames(tot_ret_m) <- c("Date", "Mkt", names)
tot_ret_m <- tot_ret_m %>% mutate(Date = FF_monthly$Date)

tot_ret_m[1,-1] <- 1

tot_ret_cost_m <- data.frame(matrix(ncol = length(names) + 2, nrow = n_months))
colnames(tot_ret_cost_m) <- c("Date", "Mkt", names)
tot_ret_cost_m <- tot_ret_cost_m %>% mutate(Date = FF_monthly$Date)

tot_ret_cost_m[1,-1] <- 1

for (i in 2:n_months) {
  tot_ret_m$Mkt[i] <- tot_ret_m$Mkt[i-1] * 
    (1 + (returns_m$Mkt[i-1]/100))
  tot_ret_cost_m$Mkt[i] <- tot_ret_m$Mkt[i]
  for (j in 1:length(names)) {
    tot_ret_m[i,names[j]] <- tot_ret_m[i-1, names[j]] * 
      (1 + (returns_m[i-1, names[j]]/100))
    tot_ret_cost_m[i,names[j]] <- tot_ret_cost_m[i-1, names[j]] * 
      (1 + (returns_cost_m[i-1, names[j]]/100))
    }
}

tot_ret_m[n_months,]
tot_ret_cost_m[n_months,]

# ***** Plot market and VM returns on log scale *****
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(tot_ret_m, aes(x = Date)) +
  geom_line(aes(y=Mkt, color = "Buy and Hold")) +
  geom_line(aes(y=var_managed, color = "Realized Variance")) +
  geom_line(aes(y=ARMA_var_managed, color = "ARMA")) +
  geom_line(aes(y=EWMA_var_managed, color = "EWMA")) +
  geom_line(aes(y=GARCH_var_managed, color = "GARCH")) +
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
  ggtitle("Cumulative Performance") + 
  xlab("Date") +
  ylab("Performance") +
  theme_stata() + 
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

# Compute Alpha and Ratios
output_names <- c("alpha_mkt", "beta_mkt", "R^2_mkt", "RMSE", "SR", "Appr_Ratio", "alpha_FF3")
reg_mkt_m <- vector(mode = "list", length = length(names))
reg_FF3_m <- vector(mode = "list", length = length(names))
reg_mkt_cost_m <- vector(mode = "list", length = length(names))
reg_FF3_cost_m <- vector(mode = "list", length = length(names))
b <- trading_months * (returns_m$Mkt - returns_m$rf)
b_rec <- trading_months * (returns_m$Mkt - returns_m$rf) * recession$indicator[-1]
rec <- recession$indicator[-1]
b1 <- trading_months * (FF_monthly$SMB[-1])
b2 <- trading_months * (FF_monthly$HML[-1])

for (i in 1:length(names)) {
  a <- trading_months * (returns_m[, names[i]] - returns_m$rf)
  reg_mkt_m[[i]] <- lm(a ~ b)
  a <- trading_months * (returns_cost_m[, names[i]] - returns_m$rf)
  reg_mkt_cost_m[[i]] <- lm(a ~ b)
}

for (i in 1:length(names)) {
  a <- trading_months * (returns_m[, names[i]] - returns_m$rf)
  reg_FF3_m[[i]] <- lm(a ~ b + b1 + b2)
  a <- trading_months * (returns_cost_m[, names[i]] - returns_m$rf)
  reg_FF3_cost_m[[i]] <- lm(a ~ b + b1 + b2)
}

reg_output_m <- data.frame(matrix(ncol = length(names), nrow = length(output_names)))
colnames(reg_output_m) <- names
rownames(reg_output_m) <- output_names

reg_output_cost_m <- data.frame(matrix(ncol = length(names), nrow = length(output_names)))
colnames(reg_output_cost_m) <- names
rownames(reg_output_cost_m) <- output_names

for (i in 1:length(names)) {
  reg_output_m["alpha_mkt", i] <- reg_mkt_m[[i]]$coefficients[1]
  reg_output_m["beta_mkt", i] <- reg_mkt_m[[i]]$coefficients[2]
  reg_output_m["R^2_mkt", i] <- summary(reg_mkt_m[[i]])$r.squared
  reg_output_m["RMSE", i] <- stats::sigma(reg_mkt_m[[i]])
  reg_output_m["SR", i] <- trading_months * 
    (mean(returns_m[,names[i]] - returns_m$rf)) / 
    (sqrt(trading_months) * sd(returns_m[,names[i]] - returns_m$rf))
  reg_output_m["Appr_Ratio", i] <- sqrt(trading_months) * 
    reg_output_m["alpha_mkt", i] / reg_output_m["RMSE", i]
  reg_output_m["alpha_FF3", i] <- reg_FF3_m[[i]]$coefficients[1]
}

for (i in 1:length(names)) {
  reg_output_cost_m["alpha_mkt", i] <- reg_mkt_cost_m[[i]]$coefficients[1]
  reg_output_cost_m["beta_mkt", i] <- reg_mkt_cost_m[[i]]$coefficients[2]
  reg_output_cost_m["R^2_mkt", i] <- summary(reg_mkt_cost_m[[i]])$r.squared
  reg_output_cost_m["RMSE", i] <- stats::sigma(reg_mkt_cost_m[[i]])
  reg_output_cost_m["SR", i] <- trading_months * 
    (mean(returns_cost_m[,names[i]] - returns_m$rf)) / 
    (sqrt(trading_months) * sd(returns_cost_m[,names[i]] - returns_m$rf))
  reg_output_cost_m["Appr_Ratio", i] <- sqrt(trading_months) * 
    reg_output_cost_m["alpha_mkt", i] / reg_output_cost_m["RMSE", i]
  reg_output_cost_m["alpha_FF3", i] <- reg_FF3_cost_m[[i]]$coefficients[1]
}

stargazer(reg_output_m, type = "text", summary = FALSE, out = "test.htm")
stargazer(reg_output_cost_m, type = "text", summary = FALSE)

# Compute Breakeven Cost
breakeven_function <- function(cost, strategy) {
  returns_be_m <- c(1:(n_months-1))
  for (i in 1:(n_months-1)) {
      returns_be_m[i] <- returns_m[i,names[strategy]] - 
        w_abs_m[i,names[strategy]] * cost
  }
  a <- trading_months * (returns_be_m - returns_m$rf)
  return(lm(a~b)$coefficients[1])
}

cost_m <- data.frame(matrix(nrow = length(names), ncol = 6))
colnames(cost_m) <- c("alpha", "|Delta w|", "1bps", "10bps", "14bps", "Break Even")
rownames(cost_m) <- names

for (i in 1:length(names)) {
  cost_m[i,1] <- breakeven_function(0, i)
  cost_m[i,2] <- mean(w_abs_m[,names[i]])
  cost_m[i,3] <- breakeven_function(0.01, i)
  cost_m[i,4] <- breakeven_function(0.1, i)
  cost_m[i,5] <- breakeven_function(0.14, i)
  cost_m[i,6] <- uniroot(breakeven_function, strategy = i, lower = 0, upper = 100)$root
}

stargazer(cost_m, type = "text", summary = FALSE, out = "test.htm")

# Leverage Constraint
leverage_function <- function(leverage, strategy) {
  returns_le_m <- c(1:(n_months-1))
  for (i in 1:(n_months-1)) {
    returns_le_m[i] <- min(weights_m[i, names[strategy]], 1 + leverage) *
      (returns_m$Mkt[i] - returns_m$rf[i]) + returns_m$rf[i]
  }
  a <- trading_months * (returns_le_m - returns_m$rf)
  return(lm(a~b)$coefficients[1])
}

leverage_m <- data.frame(matrix(nrow = length(names), ncol = 5))
colnames(leverage_m) <- c("alpha", "Median Weight", "100%", "50%", "0%")
rownames(leverage_m) <- names

for (i in 1:length(names)) {
  leverage_m[i,1] <- breakeven_function(0, i)
  leverage_m[i,2] <- median(weights_m[,names[i]])
  leverage_m[i,3] <- leverage_function(1, i)
  leverage_m[i,4] <- leverage_function(0.5, i)
  leverage_m[i,5] <- leverage_function(0, i)
}

stargazer(leverage_m, type = "text", summary = FALSE, out = "test.htm")

# Recession
reg_rec_m <- vector(mode = "list", length = length(names))

for (i in 1:length(names)) {
  a <- trading_months * (returns_m[, names[i]] - returns_m$rf)
  reg_rec_m[[i]] <- lm(a ~ rec + b + b_rec)
}

stargazer(reg_rec_m[[1]], reg_rec_m[[2]], reg_rec_m[[3]], reg_rec_m[[4]], 
          dep.var.labels = names,covariate.labels = c("Recession Indicator", "Market",
                                                     "Market:Recession"),
            type = "text", out = "test.htm", keep.stat = c("n", "rsq"))

# Analysis of Results
stat_names <- c("Mean", "SD", "Skewness", "Kurtosis", "Min", "Q1", "Median",
                "Q3", "Max")
stat_ret_m <- data.frame(matrix(ncol = length(names), nrow = length(stat_names),
                            dimnames = list(stat_names, names)))
stat_ret_m[1,] <- apply(returns_m[,names], 2, mean)
stat_ret_m[2,] <- apply(returns_m[,names], 2, sd)
stat_ret_m[3,] <- apply(returns_m[,names], 2, skewness)
stat_ret_m[4,] <- apply(returns_m[,names], 2, kurtosis)
stat_ret_m[5:9,] <- as.data.frame(apply(returns_m[,names], 2, quantile))

stat_ret_cost_m <- data.frame(matrix(ncol = length(names), nrow = length(stat_names),
                                dimnames = list(stat_names, names)))
stat_ret_cost_m[1,] <- apply(returns_cost_m[,names], 2, mean)
stat_ret_cost_m[2,] <- apply(returns_cost_m[,names], 2, sd)
stat_ret_cost_m[3,] <- apply(returns_cost_m[,names], 2, skewness)
stat_ret_cost_m[4,] <- apply(returns_cost_m[,names], 2, kurtosis)
stat_ret_cost_m[5:9,] <- as.data.frame(apply(returns_cost_m[,names], 2, quantile))

stargazer(stat_ret_m, type = "text",summary = FALSE, out = "table.htm")
round(stat_ret_cost_m, 4)

stat_weight_m <- data.frame(matrix(ncol = length(names), nrow = length(stat_names),
                                dimnames = list(stat_names, names)))
stat_weight_m[1,] <- apply(weights_m[,names], 2, mean)
stat_weight_m[2,] <- apply(weights_m[,names], 2, sd)
stat_weight_m[3,] <- apply(weights_m[,names], 2, skewness)
stat_weight_m[4,] <- apply(weights_m[,names], 2, kurtosis)
stat_weight_m[5:9,] <- as.data.frame(apply(weights_m[,names], 2, quantile))

round(stat_weight_m, 4)

ggplot(returns_m) +
  geom_density(aes(x=Mkt, color = "Buy and Hold"), adjust = 1) +
  geom_density(aes(x=var_managed, color = "Realized Variance"), adjust = 1) +
  geom_density(aes(x=ARMA_var_managed, color = "ARMA"), adjust = 1) +
  geom_density(aes(x=EWMA_var_managed, color = "EWMA"), adjust = 1) +
  geom_density(aes(x=GARCH_var_managed, color = "GARCH"), adjust = 1) +
  xlim(-40,40) +
  theme_stata() +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

returns_rolling_m <- data.frame(matrix(nrow = n_months - 1, ncol = length(names) + 2))
colnames(returns_rolling_m) <- c("Date", "Mkt", names)
returns_rolling_m$Date <- returns_m$Date

for (i in 1:(n_months -1 )) {
  returns_rolling_m[i,"Mkt"] <- mean(returns_m[(i-min(i,11)):i,"Mkt"])
  for (j in 1:length(names)) {
    returns_rolling_m[i,names[j]] <- mean(returns_m[(i-min(i,11)):i,names[j]])
  }
}

ggplot(returns_rolling_m, aes(x=Date)) +
  geom_line(aes(y = Mkt, color = "Buy and Hold")) +
  geom_line(aes(y = var_managed, color = "Realized Variance")) +
  geom_line(aes(y = ARMA_var_managed, color = "ARMA")) +
  geom_line(aes(y = EWMA_var_managed, color = "EWMA")) +
  geom_line(aes(y = GARCH_var_managed, color = "GARCH")) +
  theme_stata() +
  ggtitle("Rolling One Year Return") +
  ylab("") +
  scale_color_manual(name = "", values = c("Buy and Hold" = "black",
                                           "Realized Variance" = "red", "ARMA" = "blue",
                                           "EWMA" = "green", "GARCH" = "yellow"))
# Calculate Drawdowns
for (i in 1:nrow(tot_ret_m)) {
  if (i==1) {
    tot_ret_m$DD_Mkt[i] <- 0
    tot_ret_m$DD_var_managed[i] <- 0
    tot_ret_m$DD_ARMA_var_managed[i] <- 0
    tot_ret_m$DD_EWMA_var_managed[i] <- 0
    tot_ret_m$DD_GARCH_var_managed[i] <- 0
  }
  else {
    tot_ret_m$DD_Mkt[i] <- min(0,(tot_ret_m$Mkt[i]-max(tot_ret_m$Mkt[1:i]))/
                                 max(tot_ret_m$Mkt[1:i]))
    tot_ret_m$DD_var_managed[i] <- 
      min(0,(tot_ret_m$var_managed[i]-max(tot_ret_m$var_managed[1:i]))/
                                         max(tot_ret_m$var_managed[1:i]))
    tot_ret_m$DD_ARMA_var_managed[i] <- 
      min(0,(tot_ret_m$ARMA_var_managed[i]-max(tot_ret_m$ARMA_var_managed[1:i]))/
                                              max(tot_ret_m$ARMA_var_managed[1:i]))
    tot_ret_m$DD_EWMA_var_managed[i] <- 
      min(0,(tot_ret_m$EWMA_var_managed[i]-max(tot_ret_m$EWMA_var_managed[1:i]))/
                                              max(tot_ret_m$EWMA_var_managed[1:i]))
    tot_ret_m$DD_GARCH_var_managed[i] <- 
      min(0,(tot_ret_m$GARCH_var_managed[i]-max(tot_ret_m$GARCH_var_managed[1:i]))/
                                               max(tot_ret_m$GARCH_var_managed[1:i]))
  }
}

ggplot(tot_ret_m, aes(x = Date)) +
  geom_line(aes(y=DD_Mkt, color = "Buy and Hold")) +
  geom_line(aes(y=DD_var_managed, color = "Realized Variance")) +
  geom_line(aes(y=DD_ARMA_var_managed, color = "ARMA")) +
  geom_line(aes(y=DD_EWMA_var_managed, color = "EWMA")) +
  geom_line(aes(y=DD_GARCH_var_managed, color = "GARCH")) +
  scale_x_date(limits = as.Date(c("1926-7-1", "2019-8-1")),
               expand = c(0,0),
               breaks = dates,
               labels = format(dates, "%Y"),
               minor_breaks = NULL) +
  theme(legend.position="bottom") +
  theme_stata() +
  ggtitle("Drawdown") + 
  xlab("Date") +
  ylab("Drawdown") +
  scale_color_manual(name = "", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

# Source of Alpha
strat_reg_m <- vector(mode = "list", length = 3)
m <- trading_months * (returns_m[, names[1]] - returns_m$rf)
strat_reg_m[[1]] <- lm((trading_months * (returns_m[, names[2]] - returns_m$rf)) ~ m + b)
strat_reg_m[[2]] <- lm((trading_months * (returns_m[, names[3]] - returns_m$rf)) ~ m + b)
strat_reg_m[[3]] <- lm((trading_months * (returns_m[, names[4]] - returns_m$rf)) ~ m + b)

stargazer(strat_reg_m[[1]], strat_reg_m[[2]], strat_reg_m[[3]],
          type = "text", dep.var.labels = c("ARMA","EWMA", "GARCH"), 
          covariate.labels = c("Market", "Variance Managed"), out = "table.htm",
          keep.stat = c("n", "rsq"))

# Considering costs
regression_function <- function(cost, strategy) {
  returns_var_m <- c(1:(n_months-1))
  for (i in 1:(n_months-1)) {
    returns_var_m[i] <- returns_m[i,names[1]] - 
      w_abs_m[i,names[1]] * cost
  }
  
  returns_c_m <- c(1:(n_months-1))
  for (i in 1:(n_months-1)) {
    returns_c_m[i] <- returns_m[i,names[strategy]] - 
      w_abs_m[i,names[strategy]] * cost
  }
  a <- trading_months * (returns_c_m - returns_m$rf)
  m <- trading_months * (returns_var_m - returns_m$rf)
  return(lm(a ~ m + b))
}
cost_vector <- c(0.01, 0.1, 0.14, 0.3)
cost_m <- data.frame(matrix(nrow = length(cost_vector), ncol = length(names[-1])))
colnames(cost_m) <- names[-1]
rownames(cost_m) <- paste(cost_vector * 100, "bps")
for (i in 1:length(names[-1])) {
  for (j in 1: length(cost_vector)) {
    cost_m[j,i] <- paste(round(coefficients(regression_function(cost_vector[j],i+1))[1],2),
                         " (", 
                         round(summary(regression_function(cost_vector[j],i+1))$coefficient[1,4], 3),
                         ")", sep = "")
  }
}

stargazer(cost_m, type = "text", summary = FALSE, out = "test.htm")

# Analysis of times when var managed does work well and whether improved strategies perform better
bd_times <- returns_m
bd_times$weight <- weights_m$var_managed
bd_times$weight_ARMA <- weights_m$ARMA_var_managed
bd_times$weight_EWMA <- weights_m$EWMA_var_managed
bd_times$weight_GARCH <- weights_m$GARCH_var_managed

# Function RMSD: Root Mean Square Deviation
RMSD <- function (x) {
  n <- ncol(x)
  matrix <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      matrix[i,j] <- sqrt(mean((x[,i]-x[,j])^2))
    }
  }
  row.names(matrix) <- colnames(x)
  colnames(matrix) <- colnames(x)
  return(round(matrix,3))
}

cor(returns_m[c(2,4:7)]-returns_m[,3])

# Analysis of Correlations and RMSD
analysis_matrix_rownames_cor <- c("No. of obs.", "Corr Mkt-Var ", "Corr Mkt-ARMA ",
                              "Corr Mkt-EWMA ", "Corr Mkt-GARCH ", "Corr Var-ARMA ",
                              "Corr Var-EWMA ", "Corr Var-GARCH ")
analysis_matrix_rownames_rmsd <- c("No. of obs.", "RMSD Mkt-Var ", "RMSD Mkt-ARMA ",
                                   "RMSD Mkt-EWMA ", "RMSD Mkt-GARCH ", "RMSD Var-ARMA ",
                                   "RMSD Var-EWMA ", "RMSD Var-GARCH ")
analysis_matrix_colnames <- c("Whole Sample", "| Mkt < 0%, w > 1", "| Mkt < -5%, w > 1",
                              "| Var < -10%, w > 1", "| Var < -15%, w > 1")
analysis_matrix_cor <- data.frame(matrix(nrow = length(analysis_matrix_rownames_cor), 
                                     ncol = length(analysis_matrix_colnames)))
analysis_matrix_rmsd <- data.frame(matrix(nrow = length(analysis_matrix_rownames_rmsd), 
                                         ncol = length(analysis_matrix_colnames)))
rownames(analysis_matrix_cor) <- analysis_matrix_rownames_cor
rownames(analysis_matrix_rmsd) <- analysis_matrix_rownames_rmsd
colnames(analysis_matrix_cor) <- analysis_matrix_colnames
colnames(analysis_matrix_rmsd) <- analysis_matrix_colnames

# Whole Sample
analysis_matrix_cor[1,1] <- nrow(returns_m)
analysis_matrix_rmsd[1,1] <- nrow(returns_m)

for (i in 1:7) {
  if (i <= 4) {
    col = 1
    j = i
  } else {
    col = 2
    j = i %% 4
  }
  analysis_matrix_cor[i + 1, 1] <- cor(returns_m[c(2,4:7)])[j + col, col]
  analysis_matrix_rmsd[i + 1, 1] <- RMSD(returns_m[c(2,4:7)])[j + col, col]
}

# Extreme Periods
extreme_periods <- list()
extreme_periods[[1]] <- filter(bd_times, Mkt < 0 & weight > 0)
extreme_periods[[2]] <- filter(bd_times, Mkt < -5 & weight > 1)
extreme_periods[[3]] <- filter(bd_times, var_managed < -10 & weight > 1)
extreme_periods[[4]] <- filter(bd_times, var_managed < -15 & weight > 1)

for (period in 1:length(extreme_periods)) {
  analysis_matrix_cor[1, period + 1] <- nrow(extreme_periods[[period]])
  analysis_matrix_rmsd[1, period + 1] <- nrow(extreme_periods[[period]])
  for (i in 1:7) {
    if (i <= 4) {
      col = 1
      j = i
    } else {
      col = 2
      j = i %% 4
    }
    analysis_matrix_cor[i + 1, period + 1] <- 
      cor(extreme_periods[[period]][c(2,4:7)])[j + col, col]
    analysis_matrix_rmsd[i + 1, period + 1] <- 
      RMSD(extreme_periods[[period]][c(2,4:7)])[j + col, col]
  }
}

stargazer(analysis_matrix_cor, summary = FALSE, type = "text", out = "test.htm")
stargazer(analysis_matrix_rmsd, summary = FALSE, type = "text", out = "test.htm")

# Analysis of Delta w
analysis_matrix_rownames_2 <- c("Var", "ARMA", "EWMA", "GARCH")
analysis_matrix_colnames_2 <- c("w > 1", "| Mkt < 0%, w > 1", "| Mkt < -5%, w > 1", 
                                "| Var < -10%, w > 1", "| Var < -15%, w > 1")
analysis_matrix_w <- data.frame(matrix(nrow = length(analysis_matrix_rownames_2), 
                                       ncol = length(analysis_matrix_colnames_2)))
rownames(analysis_matrix_w) <- analysis_matrix_rownames_2
colnames(analysis_matrix_w) <- analysis_matrix_colnames_2

analysis_matrix_rel_w <- data.frame(matrix(nrow = length(analysis_matrix_rownames_2), 
                                       ncol = length(analysis_matrix_colnames_2)))
rownames(analysis_matrix_rel_w) <- analysis_matrix_rownames_2
colnames(analysis_matrix_rel_w) <- analysis_matrix_colnames_2

analysis_matrix_w[1, 1] <- mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_w[2, 1] <- mean(filter(bd_times, weight > 1)$weight_ARMA)
analysis_matrix_w[3, 1] <- mean(filter(bd_times, weight > 1)$weight_EWMA)
analysis_matrix_w[4, 1] <- mean(filter(bd_times, weight > 1)$weight_GARCH)

for (period in 1:length(extreme_periods)) {
  analysis_matrix_w[1, period + 1] <- mean(extreme_periods[[period]]$weight)
  analysis_matrix_w[2, period + 1] <- mean(extreme_periods[[period]]$weight_ARMA)
  analysis_matrix_w[3, period + 1] <- mean(extreme_periods[[period]]$weight_EWMA)
  analysis_matrix_w[4, period + 1] <- mean(extreme_periods[[period]]$weight_GARCH)
}

analysis_matrix_rel_w[1, 1] <- 
  mean(filter(bd_times, weight > 1)$weight - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_rel_w[2, 1] <- 
  mean(filter(bd_times, weight > 1)$weight_ARMA - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_rel_w[3, 1] <- 
  mean(filter(bd_times, weight > 1)$weight_EWMA - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_rel_w[4, 1] <- 
  mean(filter(bd_times, weight > 1)$weight_GARCH - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)

for (period in 1:length(extreme_periods)) {
  analysis_matrix_rel_w[1, period + 1] <- 
    mean(extreme_periods[[period]]$weight - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
  analysis_matrix_rel_w[2, period + 1] <- 
    mean(extreme_periods[[period]]$weight_ARMA - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
  analysis_matrix_rel_w[3, period + 1] <- 
    mean(extreme_periods[[period]]$weight_EWMA - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
  analysis_matrix_rel_w[4, period + 1] <- 
    mean(extreme_periods[[period]]$weight_GARCH - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
}

stargazer(analysis_matrix_w, summary = FALSE, type = "text", out = "test.htm")
stargazer(analysis_matrix_rel_w, summary = FALSE, type = "text", out = "test.htm")

################################################################################
#******************************** Custom Level ********************************#

# Set Up List to Store Outputs from Different Frequencies
min_frequ <- 5
max_frequ <- 50
reg_output_c <- vector(mode = "list", length = max_frequ - min_frequ)
reg_output_cost_c <- vector(mode = "list", length = max_frequ - min_frequ)
quantile_c <- vector(mode = "list", length = max_frequ - min_frequ)

# Loop to Test Different Frequencies

for (frequ in min_frequ:max_frequ) {
  # Calculate Custom Variances
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
  returns_c <- data.frame(matrix(ncol = 3 + length(names), nrow = (n_custom-1)))
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
  u_sq_c <- c(1:n_custom)
  ret_temp <- 0
  rf_temp <- 0
  for (a in 1:frequ) {
    ret_temp <- ((1 + FF_daily$Mkt[a]/100)*(1 + ret_temp/100) - 1)*100
    rf_temp <- ((1 + FF_daily$RF[a]/100)*(1 + rf_temp/100) - 1)*100
  }
  u_sq_c[1] <- log(1 + (ret_temp - rf_temp) / 100)^2
  u_sq_c[2:n_custom] <- log( 1 + returns_c$`Mkt-RF` / 100)^2
  
  # Scenario 2: ARIMA Model
  variance_ts_c <- xts(var_c$variance, order.by = var_c$Date)
  
  Acf(variance_ts_c, lag = 22)
  Pacf(variance_ts_c, lag = 22)
  
  ARMA_model_c <- auto.arima(variance_ts_c, stepwise = FALSE, approximation = FALSE)
  jarque.test(as.vector(ARMA_model_c$residuals))
  print(ARMA_model_c)
  
  var_c$ARMA_var <- c(fitted(ARMA_model_c)[-1], forecast(ARMA_model_c, h = 1)$mean)
  
  # Scenario 3: EWMA Model
  EWMA_function_c <- function(lambda)
  {
    EWMA_var <- c(1:n_custom)
    EWMA_var[1] <- u_sq_c[1]
    for(i in 2:n_custom) {
      EWMA_var[i] <- lambda * EWMA_var[i-1] + (1 - lambda) * u_sq_c[i]
    }
    EWMA_likelihood <- c(1:(n_custom-1))
    for(i in 1:(n_custom-1)) {
      EWMA_likelihood[i] <- -log(EWMA_var[i])-u_sq_c[i+1]/EWMA_var[i]
    }
    return (sum(EWMA_likelihood))
  }
  ewma_max_c <- optimize(EWMA_function_c, interval = c(0, 1), maximum = TRUE, 
                         tol = 0.000000000000001)
  lambda_c <- ewma_max_c$maximum
  
  # Calculate EWMA Variance
  var_c$EWMA_var <- c(1:n_custom)
  var_c$EWMA_var[1] <- var_c$variance[1]
  for (i in 2:n_custom) {
    var_c$EWMA_var[i] <- lambda_c * var_c$EWMA_var[i - 1] + (1 - lambda_c) * 
      var_c$variance[i]
  }

  # Scenario 4: GARCH
  # Optimize Parameters
  GARCH_function_c <- function(alpha, beta)
  {
    omega <- max(0,mean(u_sq_c)*(1-alpha-beta))
    GARCH_var <- c(1:n_custom)
    GARCH_var[1] <- u_sq_c[1]
    for(i in 2:n_custom) {
      GARCH_var[i] <- omega + beta*GARCH_var[i-1] + alpha*u_sq_c[i]
    }
    GARCH_likelihood <- c(1:(n_custom - 1))
    for(i in 1:(n_custom-1)) {
      GARCH_likelihood[i] <- -log(GARCH_var[i])-u_sq_c[i+1]/GARCH_var[i]
    }
    return (sum(GARCH_likelihood))
  }
  GARCH_max_c <- optimx(c(0.1, 0.9), function(x) GARCH_function_c(x[1], x[2]), 
                      method = "Nelder-Mead", control = list(maximize = TRUE))
  alpha_c <- GARCH_max_c$p1
  beta_c <- GARCH_max_c$p2
  omega_c <- max(0,mean(u_sq_c)*(1-alpha_c-beta_c))
  
  # Calculate GARCH Variance
  var_c$GARCH_var <- c(1:n_custom)
  var_c$GARCH_var[1] <- u_sq_c[1]
  for (i in 2:n_custom) {
    var_c$GARCH_var[i] <- omega_c +
      alpha_c*u_sq_c[i] + beta_c*var_c$GARCH_var[i-1]
  }
  
  # Set Scale Denominator
  denom_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
  colnames(denom_c) <- names
  
  denom_c[,1] <- var_c$variance[-n_custom]
  denom_c[,2] <- var_c$ARMA_var[-n_custom]
  denom_c[,3] <- var_c$EWMA_var[-n_custom]
  denom_c[,4] <- var_c$GARCH_var[-n_custom]
  
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
  trading_cost <- 0.1
  w_abs_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
  colnames(w_abs_c) <- names
  w_abs_c[1,] <- 1
  
  returns_cost_c <- data.frame(matrix(ncol = length(names) + 1, nrow = (n_custom-1)))
  colnames(returns_cost_c) <- c("Date", names)
  returns_cost_c <- returns_cost_c %>% mutate(Date = FF_daily$Date[1:(n_custom-1)])
  
  weights_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
  colnames(weights_c) <- names
  
  for (i in 1:(n_custom-1)) {
    for (j in 1:length(names)) {
      weights_c[i,j] <- (c_c[j] / denom_c[i,j])
      if (i != 1) {
        w_abs_c[i,names[j]] <- abs(weights_c[i,j] - weights_c[i-1,j])
      }
      returns_c[i,names[j]] <- weights_c[i, j] *
        (returns_c$Mkt[i] - returns_c$RF[i]) + returns_c$RF[i]
      returns_cost_c[i,names[j]] <- returns_c[i,names[j]] - 
        w_abs_c[i,names[j]] * trading_cost
    }
  }
  
  # Check Variance and Calculate Weight Quantiles
  print(apply(returns_c[,-c(1,3,(3+length(names)+1))], 2, var))
  mod_num <- frequ + 1 - min_frequ
  quantile_c[[mod_num]] <- as.data.frame(apply(weights_c, 2, quantile, 
                                                probs = c(0.5, 0.75, 0.9, 0.99)))
  
  # Calculate Total Return
  tot_ret_c <- data.frame(matrix(ncol = length(names) + 2, nrow = n_custom))
  colnames(tot_ret_c) <- c("Date", "Mkt", names)
  tot_ret_c <- tot_ret_c %>% mutate(Date = var_c$Date)
  
  tot_ret_cost_c <- data.frame(matrix(ncol = length(names) + 2, nrow = n_custom))
  colnames(tot_ret_cost_c) <- c("Date", "Mkt", names)
  tot_ret_cost_c <- tot_ret_cost_c %>% mutate(Date = var_c$Date)
  
  tot_ret_c[1,-1] <- 1
  tot_ret_cost_c[1,-1] <- 1
  
  for (i in 2:n_custom) {
    tot_ret_c$Mkt[i] <- tot_ret_c$Mkt[i-1] * 
      (1 + (returns_c$Mkt[i-1]/100))
    tot_ret_cost_c$Mkt[i] <- tot_ret_c$Mkt[i]
    for (j in 1:length(names)) {
      tot_ret_c[i,names[j]] <- tot_ret_c[i-1, names[j]] * 
        (1 + (returns_c[i-1, names[j]]/100))
      tot_ret_cost_c[i,names[j]] <- tot_ret_cost_c[i-1, names[j]] * 
        (1 + (returns_cost_c[i-1, names[j]]/100))
    }
  }
  
  tot_ret_c[n_custom,]
  tot_ret_cost_c[n_custom,]
  
  # Compute Alpha and Ratios
  reg_mkt_c <- vector(mode = "list", length = length(names))
  reg_mkt_cost_c <- vector(mode = "list", length = length(names))
  b <- 252 / frequ * (returns_c$`Mkt-RF`)
  
  for (i in 1:length(names)) {
    a <- 252 / frequ * (returns_c[, names[i]] - returns_c$RF)
    reg_mkt_c[[i]] <- lm(a ~ b)
    a <- 252 / frequ * (returns_cost_c[, names[i]] - returns_c$RF)
    reg_mkt_cost_c[[i]] <- lm(a ~ b)
  }
  
  output_names_c <- c("alpha_mkt", "R^2_mkt", "RMSE", "Appr_Ratio")
  reg_output_c[[mod_num]] <- data.frame(matrix(ncol = length(names), 
                                               nrow = length(output_names_c)))
  colnames(reg_output_c[[mod_num]]) <- names
  rownames(reg_output_c[[mod_num]]) <- output_names_c
  
  reg_output_cost_c[[mod_num]] <- data.frame(matrix(ncol = length(names), 
                                               nrow = length(output_names_c)))
  colnames(reg_output_cost_c[[mod_num]]) <- names
  rownames(reg_output_cost_c[[mod_num]]) <- output_names_c
  
  for (i in 1:length(names)) {
    reg_output_c[[mod_num]]["alpha_mkt", i] <- reg_mkt_c[[i]]$coefficients[1]
    reg_output_c[[mod_num]]["R^2_mkt", i] <- summary(reg_mkt_c[[i]])$r.squared
    reg_output_c[[mod_num]]["RMSE", i] <- stats::sigma(reg_mkt_c[[i]])
    reg_output_c[[mod_num]]["Appr_Ratio", i] <- sqrt(252 / frequ) * 
      reg_output_c[[mod_num]]["alpha_mkt", i] / reg_output_c[[mod_num]]["RMSE", i]
  }
  
  for (i in 1:length(names)) {
    reg_output_cost_c[[mod_num]]["alpha_mkt", i] <- reg_mkt_cost_c[[i]]$coefficients[1]
    reg_output_cost_c[[mod_num]]["R^2_mkt", i] <- summary(reg_mkt_cost_c[[i]])$r.squared
    reg_output_cost_c[[mod_num]]["RMSE", i] <- stats::sigma(reg_mkt_cost_c[[i]])
    reg_output_cost_c[[mod_num]]["Appr_Ratio", i] <- sqrt(252 / frequ) * 
      reg_output_cost_c[[mod_num]]["alpha_mkt", i] / reg_output_cost_c[[mod_num]]["RMSE", i]
  }
  
  round(reg_output_c[[mod_num]], 2)
  round(reg_output_cost_c[[mod_num]], 2)
  
  # plot log scale
  ggplot(tot_ret_c, aes(x = Date)) +
    geom_line(aes(y=Mkt, color = "Buy and Hold")) +
    geom_line(aes(y=var_managed, color = "Realized Variance")) +
    geom_line(aes(y=ARMA_var_managed, color = "ARMA")) +
    geom_line(aes(y=EWMA_var_managed, color = "EWMA")) +
    geom_line(aes(y=GARCH_var_managed, color = "GARCH")) +
    scale_y_log10() +
    theme_stata() + 
    ggtitle("Cumulative Performance") + 
    xlab("Date") +
    ylab("Performance") +
    scale_color_manual(name = "Strategies", 
                       values = c("Buy and Hold" = "black", 
                                  "Realized Variance" = "red", "ARMA" = "blue", 
                                  "EWMA" = "green", "GARCH" = "yellow"))
}


################################################################################
#******************************** Daily Level ********************************#

# Calculate Daily Variances
var_d <- data.frame(matrix(ncol = 2, nrow = n_days))
colnames(var_d) <- c("Date", "variance")
var_d$variance <- FF_daily$u_sq

var_d <- var_d %>% mutate(volatility = sqrt(variance), Date = FF_daily$Date)

# simple average
var_d$simple_6 <- c(1:n_days)
var_d$simple_6[1] <- var_d$variance[1]
for (i in 2:n_days) {
  var_d$simple_6[i] <- mean(var_d$variance[(i-min(6,i)+1):i])
}

# Scenario 2: ARIMA Model
variance_ts_d <- xts(var_d$variance, order.by = FF_daily$Date)

Acf(variance_ts_d, lag = 22)
Pacf(variance_ts_d, lag = 22)

ARMA_model_d <- auto.arima(variance_ts_d, stepwise = FALSE, approximation = FALSE)
jarque.test(as.vector(ARMA_model_d$residuals))

var_d$ARMA_var <- c(fitted(ARMA_model_d)[-1], forecast(ARMA_model_d, h = 1)$mean)

# Scenario 3: EWMA Model
EWMA_function_d <- function(lambda)
{
  EWMA_var <- c(1:n_days)
  EWMA_var[1] <- FF_daily$u_sq[1]
  for(i in 2:n_days) {
    EWMA_var[i] <- lambda * EWMA_var[i-1] + (1 - lambda) * FF_daily$u_sq[i]
  }
  EWMA_likelihood <- c(1:(n_days-1))
  for(i in 1:(n_days-1)) {
    EWMA_likelihood[i] <- -log(EWMA_var[i])-FF_daily$u_sq[i+1]/EWMA_var[i]
  }
  return (sum(EWMA_likelihood))
}
ewma_max_d <- optimize(EWMA_function_d, interval = c(0, 1), maximum = TRUE, 
                       tol = 0.000000000000001)
lambda_d <- ewma_max_d$maximum

# Calculate EWMA Variance
var_d$EWMA_var <- c(1:n_days)
var_d$EWMA_var[1] <- var_d$variance[1]
for (i in 2:n_days) {
  var_d$EWMA_var[i] <- lambda_d * var_d$EWMA_var[i - 1] + (1 - lambda_d) * 
    var_d$variance[i]
}

# Scenario 4: GARCH
# Optimize Parameter
GARCH_function_d <- function(alpha, beta)
{
  omega <- max(0,mean(FF_daily$u_sq)*(1-alpha-beta))
  GARCH_var <- c(1:n_days)
  GARCH_var[1] <- FF_daily$u_sq[1]
  for(i in 2:n_days) {
    GARCH_var[i] <- omega + beta*GARCH_var[i-1] + alpha*FF_daily$u_sq[i]
  }
  GARCH_likelihood <- c(1:(n_days - 1))
  for(i in 1:(n_days-1)) {
    GARCH_likelihood[i] <- -log(GARCH_var[i])-FF_daily$u_sq[i+1]/GARCH_var[i]
  }
  return (sum(GARCH_likelihood))
}
GARCH_max_d <- optimx(c(0.1, 0.9), function(x) GARCH_function_d(x[1], x[2]), 
                      method = "Nelder-Mead", control = list(maximize = TRUE))
alpha_d <- GARCH_max_d$p1
beta_d <- GARCH_max_d$p2
omega_d <- max(0,mean(FF_daily$u_sq)*(1-alpha_d-beta_d))

# Calculate GARCH Variance
var_d$GARCH_var <- c(1:n_days)
var_d$GARCH_var[1] <- var_d$variance[1]
for (i in 2:n_days) {
  var_d$GARCH_var[i] <- omega_d + alpha_d*var_d$variance[i] + beta_d*var_d$GARCH_var[i-1]
}

# Set Scale Denominator
names_d <- c("ARMA_var_managed", "EWMA_var_managed", "GARCH_var_managed")

denom_d <- data.frame(matrix(ncol = length(names_d), nrow = n_days - 1))
colnames(denom_d) <- names_d

denom_d[,1] <- var_d$ARMA_var[-n_days]
denom_d[,2] <- var_d$EWMA_var[-n_days]
denom_d[,3] <- var_d$GARCH_var[-n_days]

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

# Check Variance and Calculate Weight Quantiles
print(apply(returns_d[,-c(1,3)], 2, var))
quantiles_d <- apply(weights_d, 2, quantile, probs = c(0.5, 0.75, 0.9, 0.99))

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
  reg_output_d["beta_mkt", i] <- reg_mkt_d[[i]]$coefficients[2]
  reg_output_d["R^2_mkt", i] <- summary(reg_mkt_d[[i]])$r.squared
  reg_output_d["RMSE", i] <- stats::sigma(reg_mkt_d[[i]])
  reg_output_d["SR", i] <- 252 * (mean(returns_d[,names_d[i]] - returns_d$RF)) / 
    (sqrt(252) * sd(returns_d[,names_d[i]]))
  reg_output_d["Appr_Ratio", i] <- sqrt(252) * reg_output_d["alpha_mkt", i] /
    reg_output_d["RMSE", i]
  reg_output_d["alpha_FF3", i] <- reg_FF3_d[[i]]$coefficients[1]
}

round(reg_output_d, 2)

# plot log scale
ggplot(tot_ret_d, aes(x = Date)) +
  geom_line(aes(y=Mkt, col = "Buy and Hold")) +
  geom_line(aes(y=ARMA_var_managed, col = "ARMA")) +
  geom_line(aes(y=EWMA_var_managed, col = "EWMA")) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH")) +
  scale_y_log10() +
  theme_stata() + 
  ggtitle("Cumulative Performance") + 
  xlab("Date") +
  ylab("Performance") +
  scale_color_manual(name = "Strategies", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

################################################################################
#************************** Export Tables and Graphs **************************#

# Check whether can use tidy() for data.frames, otherwise formattable

alpha_c <- data.frame(matrix(ncol = length(names), nrow = max_frequ - min_frequ))
colnames(alpha_c) <- names

alpha_cost_c <- data.frame(matrix(ncol = length(names), nrow = max_frequ - min_frequ))
colnames(alpha_cost_c) <- names

appr_c <- data.frame(matrix(ncol = length(names), nrow = max_frequ - min_frequ))
colnames(appr_c) <- names

appr_cost_c <- data.frame(matrix(ncol = length(names), nrow = max_frequ - min_frequ))
colnames(appr_cost_c) <- names

for (i in 1:(max_frequ - min_frequ + 1)) {
  for (j in 1:length(names)) {
    alpha_c[i,j] <- reg_output_c[[i]][1,j]
    appr_c[i,j] <- reg_output_c[[i]][4,j]
    alpha_cost_c[i,j] <- reg_output_cost_c[[i]][1,j]
    appr_cost_c[i,j] <- reg_output_cost_c[[i]][4,j]
  }
}

ggplot(alpha_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance")) +
  geom_line(aes(y=ARMA_var_managed, col = "ARMA")) +
  geom_line(aes(y=EWMA_var_managed, col = "EWMA")) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH")) +
  theme_stata() + 
  ggtitle("Alpha (MKT) Depending on Frequency") +
  xlab("Frequency in Days") + 
  ylab("Alpha (in %)") +
  ylim(0, 6) +
  scale_color_manual(name = "Strategies", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

ggplot(alpha_cost_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance")) +
  geom_line(aes(y=ARMA_var_managed, col = "ARMA")) +
  geom_line(aes(y=EWMA_var_managed, col = "EWMA")) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH")) +
  theme_stata() + 
  ggtitle("Alpha (MKT) Depending on Frequency") +
  xlab("Frequency in Days") + 
  ylab("Alpha (in %)") +
  ylim(0, 6) +
  scale_color_manual(name = "Strategies", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

ggplot(appr_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance")) +
  geom_line(aes(y=ARMA_var_managed, col = "ARMA")) +
  geom_line(aes(y=EWMA_var_managed, col = "EWMA")) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH")) +
  theme_stata() + 
  ggtitle("Appraisal Ratio (MKT) Depending on Frequency") +
  xlab("Frequency in Days") + 
  ylab("Appraisal Ratio") +
  ylim(0, 0.5) +
  scale_color_manual(name = "Strategies", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))

ggplot(appr_cost_c, aes(x = c(min_frequ:max_frequ))) +
  geom_line(aes(y=var_managed, col = "Realized Variance")) +
  geom_line(aes(y=ARMA_var_managed, col = "ARMA")) +
  geom_line(aes(y=EWMA_var_managed, col = "EWMA")) +
  geom_line(aes(y=GARCH_var_managed, col = "GARCH")) +
  theme_stata() + 
  ggtitle("Appraisal Ratio (MKT) Depending on Frequency") +
  xlab("Frequency in Days") + 
  ylab("Appraisal Ratio") +
  ylim(0, 0.5) +
  scale_color_manual(name = "Strategies", 
                     values = c("Buy and Hold" = "black", 
                                "Realized Variance" = "red", "ARMA" = "blue", 
                                "EWMA" = "green", "GARCH" = "yellow"))
