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

last_day <- 20151231
last_month <- 201512

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

# ***** Calculate Monthly Variances *****
trading_days <- 242
monthly_vars <- FF_daily %>%
  mutate(month = month(Date), year = year(Date)) %>%
  group_by(year, month) %>%
  summarise(variance = var(Mkt) * trading_days * count(Mkt))

monthly_vars <- monthly_vars %>% mutate(volatility = sqrt(variance))
plot(monthly_vars$volatility, type = "l") # Compare to Figure 2: Looks identical but numbers wrong by factor ~16

# try stuff out
monthly_vars <- monthly_vars %>% mutate(var1 = 1 / variance)
print(quantile(monthly_vars$var1, probs = c(0.5, 0.75, 0.9, 0.99)))
iceal_c <- 6.39 / 0.063990642 

FF_monthly$Mkt_y <- FF_monthly$Mkt * 12
c = sqrt(var(FF_monthly$Mkt[2:1074])/
           var(1/monthly_vars$variance[1:1073]*FF_monthly$Mkt[2:1074]))

weights <- c(1:1073)
vola_managed_returns <- c(1:1073)
for (month in 2:1074) {
  weights[month-1] <- c/monthly_vars$variance[month-1]
  vola_managed_returns[month-1] <- weights[month-1]*FF_monthly$Mkt[month]
}
returns <- data.frame(FF_monthly$Mkt[2:i], vola_managed_returns)

print(var(vola_managed_returns))
print(var(FF_monthly$Mkt[2:i]))
print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99)))




for (i in 1062:1074) { # Simulate for various end months
  j = i - 1
  # ***** Calculate c *****
  c = sqrt(var(FF_monthly$Mkt[2:i])/var(1/monthly_vars$variance[1:j]*FF_monthly$Mkt[2:i]))
  # Tried to hardcode c to obtain weight quantiles from paper but distribution itself seems to be unequal
  # -> Even if c is wrong, the problem occurs earlier
  
  
  # ***** Calculate weights and volatility managed returns *****
  weights <- c(1:j)
  vola_managed_returns <- c(1:j)
  
  for (month in 2:i) {
    weights[month-1] <- c/monthly_vars$variance[month-1]
    vola_managed_returns[month-1] <- weights[month-1]*FF_monthly$Mkt[month]
  }
  returns <- data.frame(FF_monthly$Mkt[2:i], vola_managed_returns)
  colnames(returns) <- c("Market Returns", "Vola Managed Returns")
  
  
  # ***** Some descriptive statistics *****
  print(var(vola_managed_returns))
  print(var(FF_monthly$Mkt[2:i]))
  print(quantile(weights, probs = c(0.5, 0.75, 0.9, 0.99))) # paper: 0.93 1.59 2.64 6.39
  # When looping through all months of 2015 as end dates, the distribution is still very similar
  # The 99% quantile does not vary much (ranges from 6.76 to 6.77)
}
plot(weights)
summary(FF_monthly$Mkt)
summary(vola_managed_returns)
sort(weights, decreasing = TRUE)[1:20]
sum(weights>=1) # No. of weights >= 1 (where leverage is needed)


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

