# EWMA Backup
j = 1
for (i in 2:n_days) {
  if (month(FF_daily$Date[i]) != month(FF_daily$Date[i-1])) {
    var_m$EWMA_var[j] <- var_d$EWMA_var[i-1] * 220000
    j <- j + 1
  }
  else if (i == n_days) {
    var_m$EWMA_var[j] <- var_d$EWMA_var[i] * 220000
  }
}



var_m$EWMA_var <- c(1:n_months)
for (i in 1:n_months) {
  if (i == 1) var_m$EWMA_var[i] <- FF_monthly$u_sq[i] * 10000
  else {
    var_m$EWMA_var[i] <- lambda_monthly * var_m$EWMA_var[i - 1] + 
      (1 - lambda_monthly) * FF_monthly$u_sq[i] * 10000
  }
}


# GARCH Backup
for (i in 2:n_days) {
  var_d$GARCH_var[i] <- omega + alpha*FF_daily$u_sq[i] + beta*var_d$GARCH_var[i-1]
}

var_m$GARCH_var <- c(1:n_months)
j = 1
for (i in 2:n_days) {
  if (month(FF_daily$Date[i]) != month(FF_daily$Date[i-1])) {
    var_m$GARCH_var[j] <- var_d$GARCH_var[i-1] * 220000
    j <- j + 1
  }
  else if (i == n_days) {
    var_m$GARCH_var[j] <- var_d$GARCH_var[i] * 220000
  }
}