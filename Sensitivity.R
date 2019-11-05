abc <- 1
abd <- 1
hallo <- matrix(nrow = 7, ncol = 7)
hallo1 <- matrix(nrow = 7, ncol = 7)
for (frequ in 4:10) {
  abc <- 1
  for (adjust_n in (c(6,12,24,36,48,60,72))) {
    

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
    
    names <- ("var_managed")
    
    # Set Scale Denominator
    denom_c <- data.frame(matrix(ncol = length(names), nrow = n_custom - 1))
    colnames(denom_c) <- names
    
    denom_c[,1] <- var_c$variance[-n_custom]

    
    # Calculate c with Midnight Formula
    c_c <- data.frame(matrix(ncol = length(names)))
    colnames(c_c) <- names
    
    for (i in 1:length(names)) {
      c_c <- sqrt(var(returns_c$`Mkt-RF`) /
                    var(returns_c$`Mkt-RF` / denom_c[,i]))
    }
    
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
            weights_c[i,j] <- weights_c_th[i,j]
          }
          else {
            weights_c[i,j] <- weights_c[(i-1),j]
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
    
    hallo[abc, abd] <- reg_mkt_c[[1]]$coefficients[1]
    hallo1[abc, abd] <- sqrt(252 / frequ) * hallo[abc,abd] /sigma(reg_mkt_c[[1]])
    
    abc <- abc + 1
    


  }
  abd <- abd +1
}