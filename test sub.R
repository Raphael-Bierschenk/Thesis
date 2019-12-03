sub_dates <- c(ymd(first_day), ymd(19570730), ymd(19880730), ymd(last_day))



last_day_sub <- ymd(first_day)
year(last_day_sub) <- year(last_day_sub) + rolling_years
first_day_sub <- ymd(first_day)


for (period in 1:(length(sub_dates)-1)){
  first_day_sub <- sub_dates[i]
  last_day_sub <- sub_dates[i+1]
  
  FF_monthly_sub <- FF_monthly %>% 
    subset(subset = date < last_day_sub & date >= first_day_sub)
  
  var_m_sub <- var_m %>%
    subset(subset = date < last_day_sub & date >= first_day_sub)
  
  returns_m_sub <- returns_m %>% 
    subset(subset = date < last_day_sub & date > first_day_sub)
  
  weights_m_sub <- weights_m %>% 
    subset(subset = date < last_day_sub & date > first_day_sub)

  
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
    weights_m_sub[,names[i]] <- 
      as.numeric(c_m_sub[i]) / var_m_sub[-n_months_sub,var_names[i]]
    returns_m_sub[,names[i]] <- 
      weights_m_sub[,names[i]] * returns_m_sub$`mkt-rf` + returns_m_sub$rf
  }
  
  # Run Regressions to Determine Alpha, Beta, etc.
  reg_mkt_m_sub <- vector(mode = "list", length = length(names))
  reg_mkt_se_m_sub <- vector(mode = "list", length = length(names))
  reg_FF3_m_sub <- vector(mode = "list", length = length(names))
  
  b_m_sub <- trading_months * FF_monthly_sub$`mkt-rf`[-1]
  b1_m_sub <- trading_months * FF_monthly_sub$SMB[-1]
  b2_m_sub <- trading_months * FF_monthly_sub$HML[-1]
  
  a_m_sub <- as.data.frame(matrix(nrow = n_months_sub - 1, ncol = length(names)))
  colnames(a_m_sub) <- names
  
  for (i in 1:length(names)) {
    a_m_sub[,names[i]] <- 
      trading_months * (returns_cost_m_sub[, names[i]] - returns_m_sub$rf)
    reg_mkt_m_sub[[i]] <- 
      lm(a_m_sub[,names[i]] ~ b_m_sub)
    reg_FF3_m_sub[[i]] <- 
      lm(a_m_sub[,names[i]] ~ b_m_sub + b1_m_sub + b2_m_sub)
    reg_mkt_se_m_sub[[i]] <- 
      reg_mkt_m_sub[[i]] %>% vcovHC(type = "HC") %>% sqrt() %>% diag()
  }
  
  # Create Ouput Table
  output_names <- c("alpha_mkt", "beta_mkt", "RMSE", "SR", "AR", "alpha_FF3",
                    "alpha_FF3_se")
  reg_output_m <- data.frame(matrix(ncol = length(names), nrow = length(output_names)))
  colnames(reg_output_m) <- names
  rownames(reg_output_m) <- output_names
  
  for (i in 1:length(names)) {
    reg_output_m["alpha_mkt", i] <- reg_mkt_m[[i]]$coefficients[1]
    reg_output_m["beta_mkt", i] <- reg_mkt_m[[i]]$coefficients[2]
    reg_output_m["RMSE", i] <- stats::sigma(reg_mkt_m[[i]])
    reg_output_m["SR", i] <- trading_months * 
      (mean(returns_m[,names[i]] - returns_m$rf)) / 
      (sqrt(trading_months) * sd(returns_m[,names[i]] - returns_m$rf))
    reg_output_m["AR", i] <- sqrt(trading_months) * 
      reg_output_m["alpha_mkt", i] / reg_output_m["RMSE", i]
    reg_output_m["alpha_FF3", i] <- reg_FF3_m[[i]]$coefficients[1]
    reg_output_m["alpha_FF3_se", i] <- coeftest(reg_FF3_m[[1]], 
                                                vcovHC(reg_FF3_m[[1]], 
                                                       type = "HC"))[1,2]
  }
  
  alpha_stars_m <- vector(length = length(names))
  for (i in 1:length(names)) {
    alpha_stars_m[i] <- ff3_alpha_stars(reg_FF3_m[[i]])
  }
  
  
}
