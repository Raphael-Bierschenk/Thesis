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
library(optimx)

# Import Data
FF_daily <- read_csv("F-F_Research_Data_Factors_daily.CSV", col_names = TRUE, skip = 3)
FF_monthly <- read_csv("F-F_Research_Data_Factors.CSV", col_names = TRUE, skip = 3)

# Rename and Crop Data
FF_daily <- FF_daily %>% rename(Date = X1)
FF_monthly <- FF_monthly %>% rename(Date = X1)

first_day <- 19260701
first_month <- 192607
last_day <- 20190831
last_month <- 201908

FF_daily <- FF_daily %>% subset(subset = Date <= last_day & Date >= first_day)
FF_monthly <- FF_monthly %>% subset(subset = Date <= last_month & Date >= first_month)

# Mutate Data
n_days <- as.integer(count(FF_daily))
n_months <- as.integer(count(FF_monthly))

FF_monthly <- FF_monthly %>% mutate(Mkt = `Mkt-RF` + RF)
FF_daily <- FF_daily %>% mutate(Mkt = `Mkt-RF` + RF)

FF_daily$Date <- ymd(FF_daily$Date)
FF_monthly$Date <- as.character(FF_monthly$Date)
FF_monthly$Date <- parse_date_time(FF_monthly$Date, "ym")
FF_monthly$Date <- as.Date(FF_monthly$Date)

FF_daily$u <- log(1+FF_daily$Mkt/100)
FF_daily$u_sq <- FF_daily$u^2
FF_monthly$u <- log(1+FF_monthly$Mkt/100)
FF_monthly$u_sq <- FF_monthly$u^2

#############################################
# Optimization of EWMA and GARCH parameters #
#############################################

# EWMA daily
ewma_function_daily <- function(lambda)
{
  ewma_variances <- c(1:nrow(FF_daily))
  ewma_variances[1] <- FF_daily$u_sq[1]
  for(i in 2:nrow(FF_daily)) {
    ewma_variances[i] <- lambda*ewma_variances[i-1] + (1-lambda)*FF_daily$u_sq[i]
  }
  ewma_likelihood <- c(1:(nrow(FF_daily)-1))
  for(i in 1:(nrow(FF_daily)-1)) {
    ewma_likelihood[i] <- -log(ewma_variances[i])-FF_daily$u_sq[i+1]/ewma_variances[i]
  }
  return (sum(ewma_likelihood))
}
print(ewma_max_daily <- optimize(ewma_function_daily, interval = c(0, 1), maximum = TRUE, tol = 0.000000000000001))
lambda_daily <- ewma_max_daily$maximum

# EWMA monthly
ewma_function_monthly <- function(lambda)
{
  ewma_variances <- c(1:nrow(FF_monthly))
  ewma_variances[1] <- FF_monthly$u_sq[1]
  for(i in 2:nrow(FF_monthly)) {
    ewma_variances[i] <- lambda*ewma_variances[i-1] + (1-lambda)*FF_monthly$u_sq[i]
  }
  ewma_likelihood <- c(1:(nrow(FF_monthly)-1))
  for(i in 1:(nrow(FF_monthly)-1)) {
    ewma_likelihood[i] <- -log(ewma_variances[i])-FF_monthly$u_sq[i+1]/ewma_variances[i]
  }
  return (sum(ewma_likelihood))
}
print(ewma_max_monthly <- optimize(ewma_function_monthly, interval = c(0, 1), maximum = TRUE, tol = 0.000000000000001))
lambda_monthly <- ewma_max_monthly$maximum

################## Alternative
ewma_function_monthly1 <- function(lambda)
{
  ewma_variances <- c(1:n_months)
  ewma_variances[1] <- var_m$variance[1]
  for(i in 2:n_months) {
    ewma_variances[i] <- lambda*ewma_variances[i-1] + (1-lambda)*var_m$variance[i]
  }
  ewma_likelihood <- c(1:(n_months-1))
  for(i in 1:(n_months-1)) {
    ewma_likelihood[i] <- -log(ewma_variances[i])-var_m$variance[i+1]/ewma_variances[i]
  }
  return (sum(ewma_likelihood))
}
print(ewma_max_monthly1 <- optimize(ewma_function_monthly1, interval = c(0, 1), maximum = TRUE, tol = 0.000000000000001))
lambda_monthly1 <- ewma_max_monthly1$maximum




# GARCH daily
long_term_variance <- mean(FF_daily$u_sq)
garch_function_daily <- function(alpha, beta)
{
  omega <- max(0,long_term_variance*(1-alpha-beta))
  garch_variances <- c(1:nrow(FF_daily))
  garch_variances[1] <- FF_daily$u_sq[1]
  for(i in 2:nrow(FF_daily)) {
    garch_variances[i] <- omega + beta*garch_variances[i-1] + alpha*FF_daily$u_sq[i]
  }
  garch_likelihood <- c(1:(nrow(FF_daily)-1))
  for(i in 1:(nrow(FF_daily)-1)) {
    garch_likelihood[i] <- -log(garch_variances[i])-FF_daily$u_sq[i+1]/garch_variances[i]
  }
  return (sum(garch_likelihood))
}
print(garch_max_daily <- optimx(c(0.1, 0.9), function(x) garch_function_daily(x[1], x[2]), 
                                method = "Nelder-Mead", control = list(maximize = TRUE)))
alpha_daily <- garch_max_daily$p1
beta_daily <- garch_max_daily$p2
omega_daily <- max(0,long_term_variance*(1-alpha_daily-beta_daily))

# GARCH monthly
long_term_variance1 <- mean(var_m$variance)
garch_function_monthly1 <- function(alpha, beta)
{
  omega <- max(0,long_term_variance1*(1-alpha-beta))
  garch_variances <- c(1:n_months)
  garch_variances[1] <- var_m$variance[1]
  for(i in 2:nrow(FF_monthly)) {
    garch_variances[i] <- omega + beta*garch_variances[i-1] + alpha*var_m$variance[i]
  }
  garch_likelihood <- c(1:(n_months - 1))
  for(i in 1:(n_months-1)) {
    garch_likelihood[i] <- -log(garch_variances[i])-var_m$variance[i+1]/garch_variances[i]
  }
  return (sum(garch_likelihood))
}
print(garch_max_monthly1 <- optimx(c(0.1, 0.9), function(x) garch_function_monthly1(x[1], x[2]), 
                                method = "Nelder-Mead", control = list(maximize = TRUE)))
alpha_monthly <- garch_max_monthly1$p1
beta_monthly <- garch_max_monthly1$p2
omega_monthly <- max(0,long_term_variance*(1-alpha_monthly-beta_monthly))
