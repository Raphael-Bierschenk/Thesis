alpha_sub1 <- data.frame(matrix(ncol = length(names) + 1,
                               nrow = (n_months - trading_months * 30)))
colnames(alpha_sub1) <- c("Volatility", names_clean)

last_day_sub <- ymd(first_day)
year(last_day_sub) <- year(last_day_sub) + rolling_years
first_day_sub <- ymd(first_day)

rolling_years <- 30
for (period in 1:(n_months - rolling_years * trading_months)){
  
  month(last_day_sub) <- month(last_day_sub) + 1
  month(first_day_sub) <- month(first_day_sub) + 1
  print(period/(n_months - rolling_years * trading_months))
  
  FF_monthly_sub <- FF_monthly %>% subset(subset = date < last_day_sub &
                                            date >= first_day_sub)
  returns_m_sub <- returns_m %>% 
    subset(subset = date < last_day_sub & date > first_day_sub)
  
  returns_cost_m_sub <- data.frame(matrix(nrow = n_months_sub - 1,
                                          ncol = length(names) + 1))
  colnames(returns_cost_m_sub) <- c("date", names)
  returns_cost_m_sub <- returns_cost_m_sub %>% mutate(date = returns_m_sub$date)
  
  weights_m_sub <- weights_m %>% 
    subset(subset = date < last_day_sub & date > first_day_sub)
  
  w_abs_m_sub <- w_abs_m %>%
    subset(subset = date < last_day_sub & date > first_day_sub)
  
  var_m_sub <- var_m %>%
    subset(subset = date < last_day_sub & date >= first_day_sub)
  
  n_months_sub <- var_m_sub %>% count() %>% as.numeric()
  
  c_m_sub <- data.frame(matrix(ncol = length(names)))
  colnames(c_m_sub) <- names
  
  a_qe_m_sub <- c(1:length(names))
  b_qe_m_sub <- c(1:length(names))
  c_qe_m_sub <- c(1:length(names))
  
  for (i in 1:length(names)) {
    a_qe_m_sub[i] <- var(FF_monthly_sub$`mkt-rf`[-1] / 
                           var_m_sub[-n_months_sub,var_names[i]])
    b_qe_m_sub[i] <- 2 * cov(FF_monthly_sub$`mkt-rf`[-1] / 
                               var_m_sub[-n_months_sub,var_names[i]],
                             FF_monthly_sub$rf[-1])
    c_qe_m_sub[i] <- var(FF_monthly_sub$rf[-1]) - var(FF_monthly_sub$mkt[-1])
    
    c_m_sub[i] <- 1 / (2 * a_qe_m_sub[i]) * 
      (-b_qe_m_sub[i] + sqrt((b_qe_m_sub[i])^2 - 4 * a_qe_m_sub[i] * c_qe_m_sub[i]))
  }

  for (i in 1:length(names)) {
    weights_m_sub[,names[i]] <- as.numeric(c_m_sub[i]) / var_m_sub[-n_months_sub,var_names[i]]
    returns_m_sub[,names[i]] <- weights_m_sub[,names[i]] * 
      returns_m_sub$`mkt-rf` + returns_m_sub$rf
    w_abs_m_sub[1,names[i]] <- 1
    w_abs_m_sub[-1,names[i]] <- weights_m_sub[,names[i]] %>% diff() %>% abs()
    returns_cost_m_sub[,names[i]] <- returns_m_sub[,names[i]] - w_abs_m_sub[,names[i]] * trading_cost
  }
  
  # Run Regressions to Determine Alpha, Beta, etc.
  reg_mkt_m_sub <- vector(mode = "list", length = length(names))
  
  b_m_sub <- trading_months * returns_m_sub$`mkt-rf`
  
  a_m_sub <- as.data.frame(matrix(nrow = n_months_sub - 1, ncol = length(names)))
  colnames(a_m_sub) <- names
  
  for (i in 1:length(names)) {
    a_m_sub[,names[i]] <- trading_months * (returns_cost_m_sub[, names[i]] - returns_m_sub$rf)
    reg_mkt_m_sub[[i]] <- lm(a_m_sub[,names[i]] ~ b_m_sub)
  }
  
  for (i in 1:length(names)) {
    alpha_sub1[period, 1] <- sd(returns_m_sub[,names[1]])
    alpha_sub1[period, i + 1] <- reg_mkt_m_sub[[i]]$coefficients[1]
  }

}
