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
library(stargazer)
library(lmtest)
library(sandwich)

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

FF_daily$u <- log(1+FF_daily$`Mkt-RF`/100)
FF_daily$u_sq <- FF_daily$u^2

# Calculate EWMA Variances
# Parameters estimated with ML for 19260701 - 20190830
lambda = 0.937456814037826
FF_daily$EWMA_vars <- c(1:nrow(FF_daily))
FF_daily$EWMA_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$EWMA_vars[i] <- lambda*FF_daily$EWMA_vars[i-1]+(1-lambda)*FF_daily$u_sq[i]
}

# Calculate GARCH Variances
# Parameters estimated with ML for 19260701 - 20190830
omega = 0.00000125947711345891
alpha = 0.0986757938205847
beta = 0.890202868683846
FF_daily$GARCH_vars <- c(1:nrow(FF_daily))
FF_daily$GARCH_vars[1] <- FF_daily$u_sq[1]
for (i in 2:nrow(FF_daily)) {
  FF_daily$GARCH_vars[i] <- omega+alpha*FF_daily$u_sq[i]+beta*FF_daily$GARCH_vars[i-1]
}

################################################################################
#                    Strategies with flexible time periods                     #
################################################################################

output_names <- c("alpha_mkt", "alpha_mkt_se", "alpha_mkt_pval", "alpha_ff3",
                  "alpha_ff3_se", "alpha_ff3_pval", "beta_mkt", "beta_mkt_se",
                  "beta_mkt_pval", "N", "R^2_mkt", "RMSE_mkt", "SR_new",
                  "Appr_Ratio", "Performance")
combined_output_flex <- data.frame(row.names = output_names)

# Determine intervals for variances
min_days_for_var = 5
intervals <- c(min_days_for_var, 11, 22, 44, 66, 126, 252, 504)

# Create data frama that contains daily variances for all strategies
daily_vars <- data.frame(0)
daily_vars <- data.frame(FF_daily$Date[-c(1:(min_days_for_var-1))],
                         FF_daily$Mkt[-c(1:(min_days_for_var-1))],
                         FF_daily$RF[-c(1:(min_days_for_var-1))],
                         FF_daily$SMB[-c(1:(min_days_for_var-1))],
                         FF_daily$HML[-c(1:(min_days_for_var-1))])
colnames(daily_vars) <- c("Date", "Mkt", "RF", "SMB", "HML")

# For example, on the fifth day, vars are calculated for days 1 to 5
for (i in 1:nrow(daily_vars))
{
  for (interval in intervals)
  {
    col_name <- paste("Var_", interval, sep = "")
    # Calculate vars for period smaller than interval by using returns from
    # day one to day i until length of interval is reached
    if (i <= interval - min_days_for_var + 1) {
      daily_vars[[col_name]][i] <- var(FF_daily$Mkt[1:(i+min_days_for_var-1)])
    }
    else {
      daily_vars[[col_name]][i] <- 
        var(FF_daily$Mkt[(i-interval+min_days_for_var):(i+min_days_for_var-1)])
    }
  }
  # Set EWMA and GARCH variance
  daily_vars$EWMA[i] <- FF_daily$EWMA_vars[i+min_days_for_var-1]*10000
  daily_vars$GARCH[i] <- FF_daily$GARCH_vars[i+min_days_for_var-1]*10000
}

ggplot(daily_vars, aes(Date)) + 
  geom_line(aes(y=sqrt(252*EWMA), colour="EWMA")) + 
  geom_line(aes(y=sqrt(252*GARCH), colour="GARCH")) + 
  geom_line(aes(y=sqrt(252*Var_5), colour="Var_5")) + 
  geom_line(aes(y=sqrt(252*Var_11), colour="Var_11")) + 
  geom_line(aes(y=sqrt(252*Var_22), colour="Var_22")) + 
  geom_line(aes(y=sqrt(252*Var_44), colour="Var_44")) + 
  geom_line(aes(y=sqrt(252*Var_66), colour="Var_66")) + 
  geom_line(aes(y=sqrt(252*Var_126), colour="Var_126")) +
  geom_line(aes(y=sqrt(252*Var_252), colour="Var_252")) +
  geom_line(aes(y=sqrt(252*Var_504), colour="Var_504"))

# Calculate percentage deviations of all variances
for (i in 6:15) {
  daily_vars[[paste(colnames(daily_vars[i]), "_perc_dev", sep = "")]][1] <- 0
}
for (i in 2:nrow(daily_vars)) {
  for (j in 16:25) {
    daily_vars[i,j] <- daily_vars[i,j-10]/daily_vars[i-1,j-10] - 1
  }
}

# Create a list containing quantiles of variances
quantiles <- list()
quantiles[[1]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(1/44, 1-1/44))
quantiles[[2]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.05, 0.95))
quantiles[[3]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.10, 0.90))
quantiles[[4]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.15, 0.85))
quantiles[[5]] <- apply(daily_vars[,16:25], 2, quantile, probs = c(0.20, 0.80))
quantiles

strategies <- colnames(daily_vars[16:25])

quantile = 1
strategy <- "Var_5_perc_dev"

for (quantile in 1:length(quantiles))
{
  for (strategy in strategies)
  {
    returns_flex <- data.frame(ymd("1900/01/01"), 0, 0, 0, 0, 0)
    colnames(returns_flex) <- c("Date", "Variance", "Mkt", "RF", "SMB", "HML")
    
    last_i = 0
    for (i in 1:nrow(daily_vars)) {
      if (daily_vars[[strategy]][i] > quantiles[[quantile]][2,strategy] || 
          daily_vars[[strategy]][i] < quantiles[[quantile]][1,strategy]) {
        mkt_temp = 0
        rf_temp = 0
        smb_temp = 0
        hml_temp = 0
        for (j in c((last_i+1):i)) {
          mkt_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + mkt_temp/100) - 1)*100
          rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
          smb_temp <- ((1 + daily_vars$SMB[j]/100)*(1 + smb_temp/100) - 1)*100
          hml_temp <- ((1 + daily_vars$HML[j]/100)*(1 + hml_temp/100) - 1)*100
        }
        df_temp <- data.frame(daily_vars$Date[i],
                              daily_vars[[substr(strategy,1,nchar(strategy)-9)]][i],
                              mkt_temp, rf_temp, smb_temp, hml_temp)
        returns_flex[nrow(returns_flex) + 1,] <- df_temp
        last_i = i
      }
    }
    # if last day's variance was not outside interval, add last period manually
    if (last_i != nrow(daily_vars)) {
      mkt_temp = 0
      rf_temp = 0
      smb_temp = 0
      hml_temp = 0
      for (j in c((last_i+1):i)) {
        mkt_temp <- ((1 + daily_vars$Mkt[j]/100)*(1 + mkt_temp/100) - 1)*100
        rf_temp <- ((1 + daily_vars$RF[j]/100)*(1 + rf_temp/100) - 1)*100
        smb_temp <- ((1 + daily_vars$SMB[j]/100)*(1 + smb_temp/100) - 1)*100
        hml_temp <- ((1 + daily_vars$HML[j]/100)*(1 + hml_temp/100) - 1)*100
      }
      df_temp <- data.frame(daily_vars$Date[i],
                            daily_vars[[substr(strategy,1,nchar(strategy)-9)]][i],
                            mkt_temp, rf_temp, smb_temp, hml_temp)
      returns_flex[nrow(returns_flex) + 1,] <- df_temp
    }
    
    returns_flex$`Mkt-RF` <- returns_flex$Mkt - returns_flex$RF
    
    returns_flex <- returns_flex[-1,]
    
    
    # Set shortcuts
    n <- nrow(returns_flex)
    denom <- returns_flex$Variance[-n]
    mktrf <- returns_flex$`Mkt-RF`[-1]
    rf <- returns_flex$RF[-1]
    
    # Calculate c
    a_qe_flex <- var(mktrf/denom)
    b_qe_flex <- 2*cov(mktrf/denom,rf)
    c_qe_flex <- var(rf)-var(mktrf+rf)
    
    c_flex <- 1/(2*a_qe_flex)*(-b_qe_flex+sqrt((b_qe_flex^2-4*a_qe_flex*c_qe_flex)))
    
    # Calculate weights
    weights_flex <- c_flex/denom
    quantile(weights_flex, probs = c(0.5, 0.75, 0.9, 0.99))
    
    returns_flex$Weight <- c(1,weights_flex)
    
    # Calculate volatility managed returns (incl. various transaction costs)
    returns_flex <- returns_flex %>%
      mutate(AbsDeltaW = abs(diff(c(0,Weight))),
             VMR = Weight*`Mkt-RF` + RF,
             VMR_1bps = Weight*(`Mkt-RF`) + RF - AbsDeltaW*0.01,
             VMR_10bps = Weight*(`Mkt-RF`) + RF - AbsDeltaW*0.10,
             VMR_14bps = Weight*(`Mkt-RF`) + RF - AbsDeltaW*0.14)
    
    # Check whether VMR and Mkt have the same variance
    apply(returns_flex[-1,c(3,10)], 2, var)
    
    # Calculate total returns (incl. various transaction costs)
    for (i in 1:n) {
      returns_flex$tot_ret_VM[i] <- c(1,returns_flex$tot_ret_VM)[i]*
        (1 + returns_flex$VMR[i]/100)
      returns_flex$tot_ret_VM_1bps[i] <- c(1,returns_flex$tot_ret_VM_1bps)[i]*
        (1 + returns_flex$VMR_1bps[i]/100)
      returns_flex$tot_ret_VM_10bps[i] <- c(1,returns_flex$tot_ret_VM_10bps)[i]*
        (1 + returns_flex$VMR_10bps[i]/100)
      returns_flex$tot_ret_VM_14bps[i] <- c(1,returns_flex$tot_ret_VM_14bps)[i]*
        (1 + returns_flex$VMR_14bps[i]/100)
    }
    
    # Annualize returns
    trading_days <- 252
    factor <- (nrow(returns_flex) - 1) / nrow(daily_vars) * trading_days
    returns_flex <- returns_flex %>%
      mutate(Mkt = factor*Mkt,
             RF = factor*RF,
             SMB = factor*SMB,
             HML = factor*HML,
             `Mkt-RF` = factor*`Mkt-RF`,
             VMR = factor*VMR,
             VMR_1bps = factor*VMR_1bps,
             VMR_10bps = factor*VMR_10bps,
             VMR_14bps = factor*VMR_14bps)
    
    # Regressions on Mkt factor
    reg_flex_mkt <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    reg_flex_mkt_1bps <- lm(VMR_1bps[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    reg_flex_mkt_10bps <- lm(VMR_10bps[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    reg_flex_mkt_14bps <- lm(VMR_14bps[-1] - RF[-1] ~ `Mkt-RF`[-1], returns_flex)
    
    robust_se_mkt <- sqrt(diag(vcovHC(model_mkt, type = "HC")))
    
    stargazer(reg_flex_EWMA, type = "text", out = "testt.htm",
              dep.var.labels = "Volatility-managed return",
              covariate.labels = "Market return",
              omit.stat = c("f", "adj.rsq"))
    
    # Regressions on FF3 factors
    reg_flex_ff3 <- lm(VMR[-1] - RF[-1] ~ `Mkt-RF`[-1] + SMB[-1] + HML[-1], returns_flex)
    reg_flex_ff3_1bps <- lm(VMR_1bps[-1] - RF[-1] ~ `Mkt-RF`[-1] + SMB[-1] + HML[-1], returns_flex)
    reg_flex_ff3_10bps <- lm(VMR_10bps[-1] - RF[-1] ~ `Mkt-RF`[-1] + SMB[-1] + HML[-1], returns_flex)
    reg_flex_ff3_14bps <- lm(VMR_14bps[-1] - RF[-1] ~ `Mkt-RF`[-1] + SMB[-1] + HML[-1], returns_flex)
  }
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
  #geom_line(aes(y=tot_ret$var_managed, x=tot_ret$Date, colour="Base Strategy")) +
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
                                                     "EWMA" = "green", "GARCH" = "blue", "Base Strategy" = "orange"))


names_flex_cost <- c("Var", "Var_1bps", "Var_10bps", "Var_14bps", "EWMA", "EWMA_1bps", "EWMA_10bps", "EWMA_14bps",
                     "GARCH", "GARCH_1bps", "GARCH_10bps", "GARCH_14bps")
reg_output_flex <- data.frame(matrix(ncol = length(names_flex_cost), nrow = length(output_names)))
colnames(reg_output_flex) <- names_flex_cost
rownames(reg_output_flex) <- output_names



output_names <- c("alpha_mkt", "alpha_mkt_se", "alpha_mkt_pval", "alpha_ff3", "alpha_ff3_se", "alpha_ff3_pval", "beta_mkt",
                  "beta_mkt_se", "beta_mkt_pval", "N", "R^2_mkt", "RMSE_mkt", "SR_new", "Appr_Ratio", "Performance")
for (i in 1:length(names_flex_cost)) {
  model_mkt <- get(paste("reg_flex_", names_flex_cost[i], sep = ""))
  model_ff3 <- get(paste("reg_flex_ff3_", names_flex_cost[i], sep = ""))
  robust_se_mkt <- sqrt(diag(vcovHC(model_mkt, type = "HC")))
  robust_se_ff3 <- sqrt(diag(vcovHC(model_ff3, type = "HC")))
  reg_output_flex["alpha_mkt", i] <- model_mkt$coefficients[1] * factors[i]
  reg_output_flex["alpha_mkt_se", i] <- coeftest(model_mkt, vcovHC(model_mkt, type = "HC"))[1,2] * factors[i]
  reg_output_flex["alpha_mkt_pval", i] <- coeftest(model_mkt, vcovHC(model_mkt, type = "HC"))[1,4]
  reg_output_flex["alpha_ff3", i] <- model_ff3$coefficients[1] * factors[i]
  reg_output_flex["alpha_ff3_se", i] <- coeftest(model_ff3, vcovHC(model_ff3, type = "HC"))[1,2] * factors[i]
  reg_output_flex["alpha_ff3_pval", i] <- coeftest(model_ff3, vcovHC(model_ff3, type = "HC"))[1,4]
  reg_output_flex["beta_mkt", i] <- model_mkt$coefficients[2]
  reg_output_flex["beta_mkt_se", i] <- coeftest(model_mkt, vcovHC(model_mkt, type = "HC"))[2,2] * factors[i]
  reg_output_flex["beta_mkt_pval", i] <- coeftest(model_mkt, vcovHC(model_mkt, type = "HC"))[2,4]
  reg_output_flex["N", i] <- sum(summary(model_mkt)$df[1:2])
  reg_output_flex["R^2_mkt", i] <- summary(model_mkt)$r.squared
  reg_output_flex["RMSE_mkt", i] <- sigma(model_mkt) * factors[i]
  reg_output_flex["SR_new", i] <- factors[i] * 
    mean(get(paste("vars_flexible_v2_", names_flex[ceiling(i/4)], sep = ""))[[c("VMR_14bps", "VMR", "VMR_1bps", "VMR_10bps")[i%%4+1]]] - 
            get(paste("vars_flexible_v2_", names_flex[ceiling(i/4)], sep = ""))$RF) / (sqrt(factors[i]) *
            sd(get(paste("vars_flexible_v2_", names_flex[ceiling(i/4)], sep = ""))[[c("VMR_14bps", "VMR", "VMR_1bps", "VMR_10bps")[i%%4+1]]] -
                                     get(paste("vars_flexible_v2_", names_flex[ceiling(i/4)], sep = ""))$RF))
  reg_output_flex["Appr_Ratio", i] <- sqrt(factors[i]) * reg_output_flex["alpha_mkt", i] / reg_output_flex["RMSE", i]
  reg_output_flex["Performance", i] <- tail(get(paste("tot_ret_VM_", names_flex[ceiling(i/4)], sep = ""))[,c(5,2,3,4)[i%%4+1]], n = 1)
}

round(reg_output_flex, 4)

for (i in 1:ncol(reg_output_flex)) {
  combined_output_flex[,ncol(combined_output_flex)+1] <- round(reg_output_flex[,i],2)
  col = ncol(combined_output_flex)
  colnames(combined_output_flex)[col] <- paste(colnames(reg_output_flex)[i], interval, sep = "_")
}

#combined_output_flex


# Market return
market_returns <- data.frame(FF_monthly$Date, c(1:nrow(FF_monthly)))
colnames(market_returns) <- c("Date", "Performance Market")
for (i in 2:nrow(FF_monthly)) {
  market_returns[i,2] <- market_returns[i-1,2]*(1 + FF_monthly$Mkt[i]/100)
}

stargazer(reg_flex_EWMA, type = "text", out = "testt.htm",
          dep.var.labels = "Volatility-managed return",
          covariate.labels = "Market return",
          omit.stat = c("f", "adj.rsq"))


# alpha and apprasial ratio table (costs <-> interval)
names_cost <- c("Default", "1 bps", "10 bps", "14 bps")
output_flex_alpha <- data.frame(matrix(ncol = length(intervals)+2, nrow = length(names_cost)))
rownames(output_flex_alpha) <- names_cost
colnames(output_flex_alpha) <- c("1 week", "1/2 months", "1 months", "2 months", "3 months", "6 months", "1 year", "2 years", "EWMA", "GARCH")
output_flex_appr_ratio <- data.frame(matrix(ncol = length(intervals)+2, nrow = length(names_cost)))
rownames(output_flex_appr_ratio) <- names_cost
colnames(output_flex_appr_ratio) <- c("1 week", "1/2 months", "1 months", "2 months", "3 months", "6 months", "1 year", "2 years", "EWMA", "GARCH")
for (i in 1:length(intervals)) {
  for (j in 1:length(names_cost)) {
    output_flex_alpha[j,i] <- combined_output_flex[1,j+(i-1)*12]
    output_flex_appr_ratio[j,i] <- combined_output_flex[14,j+(i-1)*12]
    output_flex_alpha[j,9] <- combined_output_flex[1,4+j]
    output_flex_alpha[j,10] <- combined_output_flex[1,8+j]
    output_flex_appr_ratio[j,9] <- combined_output_flex[14,4+j]
    output_flex_appr_ratio[j,10] <- combined_output_flex[14,8+j]
  }
}
# Alphas
output_flex_alpha
# Apprasial Ratios
output_flex_appr_ratio


# improve table (without costs)
output_flex <- data.frame(matrix(ncol = length(intervals)+2, nrow = length(output_names)))
rownames(output_flex) <- output_names
colnames(output_flex) <- c("1 week", "1/2 months", "1 months", "2 months", "3 months", "6 months", "1 year", "2 years", "EWMA", "GARCH")
for (i in 1:(length(intervals)+2)) {
  if(i <= 8) { output_flex[,i] <- combined_output_flex[,i*12-11] }
  if(i == 9) { output_flex[,i] <- combined_output_flex[,5] }
  if(i == 10) { output_flex[,i] <- combined_output_flex[,9] }
}
output_flex


# Graphs

# Distribution of GARCH percentage deviations
daily_vars %>%
  mutate(GARCH_perc_dev_temp = ifelse(GARCH_perc_dev > 1, 1, GARCH_perc_dev)) %>%
  ggplot(aes(GARCH_perc_dev_temp)) +
  geom_histogram(aes(log(1+GARCH_perc_dev_temp)), binwidth = 0.005) +
  geom_vline(xintercept = quantiles[1,3], color = "red", linetype = "dashed") +
  geom_vline(xintercept = quantiles[2,3], color = "red", linetype = "dashed") +
  ggtitle("Distribution of percentage deviation of variances for GARCH") + xlab("") + ylab("")


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
    model[[i-min_obs]] <- 1 # auto.arima(variance_ts[max(1,(i-max_obs)):i], d = 0)
    monthly_vars$ARMA_var[i] <- 1 # as.numeric(forecast(model[[i-min_obs]], h = 1)$mean)
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
months <- seq(as.Date("1926/7/1"), as.Date("2019/8/1"), by = "month")
scale <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,
           200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,
           20000,30000,40000,50000,60000,70000,80000,90000,100000)
dates <- seq.Date(from = as.Date("1930-1-1"), to = as.Date("2010-1-1"), by = "10 years")
ggplot(tot_ret, aes(months)) +
  geom_line(aes(y=Mkt)) +
  geom_line(aes(y=var_managed)) +
  geom_line(aes(y=var_factor)) +
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
