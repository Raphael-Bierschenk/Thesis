alpha_sub <- data.frame(matrix(ncol = length(names) + 1,
                               nrow = (n_months - trading_months * 30)))
colnames(alpha_sub) <- c("Volatility", names_clean)

last_day_sub <- ymd(first_day)
year(last_day_sub) <- year(last_day_sub) + rolling_years
first_day_sub <- ymd(first_day)

rolling_years <- 30
for (period in 1:(n_months - rolling_years * trading_months + 1)){

  FF_monthly_sub <- FF_monthly %>% subset(subset = date < last_day_sub &
                                            date >= first_day_sub)
  FF_daily_sub <- FF_daily %>% subset(subset = date < last_day_sub & 
                                        date >= first_day_sub)
  var_m_sub <- var_m %>% 
    subset(subset = date < last_day_sub & date >= first_day_sub) %>%
    select(date, var)
    
  n_days_sub <- FF_daily_sub %>% count() %>% as.numeric()
  n_months_sub <- FF_monthly_sub %>% count() %>% as.numeric()

  
  # Calculate ARIMA Variance 
  variance_ts_m_sub <- xts(var_m_sub$var, order.by = var_m_sub$date)
  
  ARIMA_model_m_sub <- auto.arima(variance_ts_m_sub,
                                  stepwise = FALSE, 
                                  approximation = FALSE)
  
  var_m_sub$ARIMA_var <- c(fitted(ARIMA_model_m_sub)[-1], 
                           forecast(ARIMA_model_m_sub, h = 1)$mean)
  
  # Optimize Lambda for EWMA Model
  EWMA_function_m_sub <- function(lambda) {
    EWMA_var <- c(1:n_months_sub)
    EWMA_var[1] <- FF_monthly_sub$u_sq[1]
    if (EWMA_var[1] == 0) EWMA_var[1] <- FF_monthly_sub$u_sq[2]
    for(i in 2:n_months_sub) {
      EWMA_var[i] <- lambda * EWMA_var[i-1] + (1 - lambda) * FF_monthly_sub$u_sq[i]
    }
    EWMA_likelihood <- c(1:(n_months_sub-1))
    for(i in 1:(n_months_sub-1)) {
      EWMA_likelihood[i] <- -log(EWMA_var[i]) - FF_monthly_sub$u_sq[i+1] / EWMA_var[i]
    }
    return (sum(EWMA_likelihood))
  }
  
  GARCH_max_m_sub <- optimize(EWMA_function_m_sub, interval = c(0, 1), maximum = TRUE, 
                              tol = 0.000000000000001)
  lambda_m_sub <- GARCH_max_m_sub$maximum
  
  # Calculate EWMA Variance and Volatility
  var_m_sub$EWMA_var <- c(1:n_months_sub)
  var_m_sub$EWMA_var[1] <- var_m_sub$var[1]
  
  for (i in 2:n_months_sub) {
    var_m_sub$EWMA_var[i] <- lambda_m * var_m_sub$EWMA_var[i - 1] + 
      (1 - lambda_m) * var_m_sub$var[i]
  }
  
  # Optimize GARCH Parameter
  GARCH_function_m_sub <- function(alpha, beta) {
    omega <- max(0,mean(FF_monthly_sub$u_sq) * (1 - alpha - beta))
    GARCH_var <- c(1:n_months_sub)
    GARCH_var[1] <- FF_monthly_sub$u_sq[1]
    if (GARCH_var[1] == 0) GARCH_var[1] <- FF_monthly_sub$u_sq[2]
    for(i in 2:n_months_sub) {
      GARCH_var[i] <- omega + beta * GARCH_var[i-1] + alpha * FF_monthly_sub$u_sq[i]
    }
    GARCH_likelihood <- c(1:(n_months_sub - 1))
    for(i in 1:(n_months_sub - 1)) {
      GARCH_likelihood[i] <- -log(GARCH_var[i]) - FF_monthly_sub$u_sq[i+1] / GARCH_var[i]
    }
    return (sum(GARCH_likelihood))
  }
  
  GARCH_max_m_sub <- optimx(c(0.1, 0.9), function(x) GARCH_function_m(x[1], x[2]), 
                            method = "Nelder-Mead", control = list(maximize = TRUE))
  
  alpha_m_sub <- GARCH_max_m_sub$p1
  beta_m_sub <- GARCH_max_m_sub$p2
  omega_m_sub <- max(0, mean(FF_monthly_sub$u_sq) * (1 - alpha_m_sub - beta_m_sub))
  
  # Calculate GARCH Variance and Volatility
  var_m_sub$GARCH_var <- c(1:n_months_sub)
  var_m_sub$GARCH_var[1] <- FF_monthly_sub$u_sq[1] + omega_m_sub
  for (i in 2:n_months_sub) {
    var_m_sub$GARCH_var[i] <- (omega_m_sub + alpha_m_sub*FF_monthly_sub$u_sq[i] + 
                                 beta_m_sub * var_m_sub$GARCH_var[i-1])
  }
  var_m_sub <- var_m_sub %>% mutate(GARCH_var = GARCH_var * 10000)
  
  # Calculate c
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
  # Calculate Weights, Absolute Weight Deviation and Returns
  weights_m_sub <- data.frame(matrix(ncol = length(names) + 1, nrow = n_months_sub - 1))
  colnames(weights_m_sub) <- c("date", names)
  weights_m_sub <- weights_m_sub %>% mutate(date = FF_monthly_sub$date[-1])
  
  returns_m_sub <- data.frame(matrix(ncol = length(names) + 4, nrow = n_months_sub - 1))
  colnames(returns_m_sub) <- c("date", "mkt", "rf", "mkt-rf", names)
  returns_m_sub[,1:4] <- returns_m %>%
    subset(subset = date < last_day_sub & date > first_day_sub) %>%
    select("date", "mkt", "rf", "mkt-rf")
  
  for (i in 1:length(names)) {
    weights_m_sub[,names[i]] <- as.numeric(c_m[i]) / var_m_sub[-n_months_sub,var_names[i]]
    returns_m_sub[,names[i]] <- weights_m_sub[,names[i]] * 
      returns_m_sub$`mkt-rf` + returns_m_sub$rf
  }
  
  
  # Run Regressions to Determine Alpha, Beta, etc.
  reg_mkt_m_sub <- vector(mode = "list", length = length(names))
  
  b_m_sub <- trading_months * returns_m_sub$`mkt-rf`
  
  a_m_sub <- as.data.frame(matrix(nrow = n_months_sub - 1, ncol = length(names)))
  colnames(a_m_sub) <- names
  
  for (i in 1:length(names)) {
    a_m_sub[,names[i]] <- trading_months * (returns_m_sub[, names[i]] - returns_m_sub$rf)
    reg_mkt_m_sub[[i]] <- lm(a_m_sub[,names[i]] ~ b_m_sub)
  }
  
  for (i in 1:length(names)) {
    alpha_sub[period, 1] <- sd(returns_m_sub$mkt)
    alpha_sub[period, i + 1] <- reg_mkt_m_sub[[i]]$coefficients[1]
  }
  
  month(last_day_sub) <- month(last_day_sub) + 1
  month(first_day_sub) <- month(first_day_sub) + 1
  print(period)
}
