

# Analysis of times when var managed does work well and whether improved strategies perform better
bd_times <- returns_m
bd_times$weight <- weights_m$var_managed
bd_times$weight_ARMA <- weights_m$ARMA_var_managed
bd_times$weight_EWMA <- weights_m$EWMA_var_managed
bd_times$weight_GARCH <- weights_m$GARCH_var_managed

# Idea: Look at correlations between Mkt and other strategies and compare it to correlations in extreme periods
# The lower the correlation to Mkt, the better the strategy because it acts less like var_managed then
# Fact: correlation between Mkt and strategy is almost equal to regression beta
# This is due to Variance of Mkt-rf and strategy almost being equal (Var(Mkt)=Var(strategy))

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
analysis_matrix_rownames <- c("No. of obs.", "Corr Mkt-Var ", "Corr Mkt-ARMA ",
                              "Corr Mkt-EWMA ", "Corr Mkt-GARCH ", "Corr Var-ARMA ",
                              "Corr Var-EWMA ", "Corr Var-GARCH ", "RMSD Mkt-Var ",
                              "RMSD Mkt-ARMA ", "RMSD Mkt-EWMA ", "RMSD Mkt-GARCH ",
                              "RMSD Var-ARMA ", "RMSD Var-EWMA ",  "RMSD Var-GARCH ")
analysis_matrix_colnames <- c("Whole Sample", "| Mkt < 0%, w > 1", "| Mkt < -5%, w > 1",
                              "| Var < -10%, w > 1", "| Var < -15%, w > 1", 
                              "| Var < -10%, w > 3", "| Var < -15%, w > 3")
analysis_matrix <- data.frame(matrix(nrow = length(analysis_matrix_rownames), 
                                     ncol = length(analysis_matrix_colnames)))
rownames(analysis_matrix) <- analysis_matrix_rownames
colnames(analysis_matrix) <- analysis_matrix_colnames

# Whole Sample
analysis_matrix[1,1] <- nrow(returns_m)
for (i in 1:7) {
  if (i <= 4) {
    col = 1
    j = i
  } else {
    col = 2
    j = i %% 4
  }
  analysis_matrix[i + 1, 1] <- cor(returns_m[c(2,4:7)])[j + col, col]
  analysis_matrix[i + 8, 1] <- RMSD(returns_m[c(2,4:7)])[j + col, col]
}
nrow(filter(bd_times, Mkt<0))

# Extreme Periods
extreme_periods <- list()
extreme_periods[[1]] <- filter(bd_times, Mkt < 0 & weight > 0)
extreme_periods[[2]] <- filter(bd_times, Mkt < -5 & weight > 1)
extreme_periods[[3]] <- filter(bd_times, var_managed < -10 & weight > 1)
extreme_periods[[4]] <- filter(bd_times, var_managed < -15 & weight > 1)
extreme_periods

for (period in 1:length(extreme_periods)) {
  analysis_matrix[1, period + 1] <- nrow(extreme_periods[[period]])
  for (i in 1:7) {
    if (i <= 4) {
      col = 1
      j = i
    } else {
      col = 2
      j = i %% 4
    }
    analysis_matrix[i + 1, period + 1] <- cor(extreme_periods[[period]][c(2,4:7)])[j + col, col]
    analysis_matrix[i + 8, period + 1] <- RMSD(extreme_periods[[period]][c(2,4:7)])[j + col, col]
  }
}

round(analysis_matrix, 4)

# Analysis of Delta w
analysis_matrix_rownames_2 <- c("Avg. w - Var", "Avg. w - ARMA", "Avg. w - EWMA", 
                                "Avg. w - GARCH", "Avg. rel. delta w - ARMA",
                                "Avg. rel. delta w - EWMA", "Avg. rel. delta w - GARCH")
analysis_matrix_colnames_2 <- c("w > 1", "| Mkt < 0%, w > 1", "| Mkt < -5%, w > 1", 
                                "| Var < -10%, w > 1", "| Var < -15%, w > 1")
analysis_matrix_2 <- data.frame(matrix(nrow = length(analysis_matrix_rownames_2), 
                                       ncol = length(analysis_matrix_colnames_2)))
rownames(analysis_matrix_2) <- analysis_matrix_rownames_2
colnames(analysis_matrix_2) <- analysis_matrix_colnames_2

analysis_matrix_2[1, 1] <- mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_2[2, 1] <- mean(filter(bd_times, weight > 1)$weight_ARMA)
analysis_matrix_2[3, 1] <- mean(filter(bd_times, weight > 1)$weight_EWMA)
analysis_matrix_2[4, 1] <- mean(filter(bd_times, weight > 1)$weight_GARCH)
analysis_matrix_2[5, 1] <- 
  mean(filter(bd_times, weight > 1)$weight_ARMA - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_2[6, 1] <- 
  mean(filter(bd_times, weight > 1)$weight_EWMA - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)
analysis_matrix_2[7, 1] <- 
  mean(filter(bd_times, weight > 1)$weight_GARCH - filter(bd_times, weight > 1)$weight) /
  mean(filter(bd_times, weight > 1)$weight)

for (period in 1:length(extreme_periods)) {
  analysis_matrix_2[1, period + 1] <- mean(extreme_periods[[period]]$weight)
  analysis_matrix_2[2, period + 1] <- mean(extreme_periods[[period]]$weight_ARMA)
  analysis_matrix_2[3, period + 1] <- mean(extreme_periods[[period]]$weight_EWMA)
  analysis_matrix_2[4, period + 1] <- mean(extreme_periods[[period]]$weight_GARCH)
  analysis_matrix_2[5, period + 1] <- 
    mean(extreme_periods[[period]]$weight_ARMA - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
  analysis_matrix_2[6, period + 1] <- 
    mean(extreme_periods[[period]]$weight_EWMA - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
  analysis_matrix_2[7, period + 1] <- 
    mean(extreme_periods[[period]]$weight_GARCH - extreme_periods[[period]]$weight) /
    mean(extreme_periods[[period]]$weight)
}

round(analysis_matrix_2, 4)



# 1
# Results for periods of negative return and weight > 1
# No. of obs. = 185
# Var and ARMA much higher corr to Mkt, ARMA higher corr to Var
# EWMA and GARCH lower corr to Mkt and lower corr to Var
#  -> behave less like Var -> better during these periods

# 2
# Results for periods of Mkt returns < -5% and weight > 1
# No. of obs. = 24
# Var and ARMA slightly higher corr to Mkt, ARMA lower corr to Var
# EWMA and GARCH basically uncorrelated to Mkt and strikingly lower corr to Var
#  -> behave much less like Var -> much better during these periods

# 3
# Results for periods of var managed returns < -10% and weight > 1
# No. of obs. = 21
# Var and ARMA higher corr to Mkt, ARMA lower corr to Var
# EWMA and GARCH strikingly lower corr to Mkt and Var
#  -> behave much less like Var -> much better during these periods

# 4
# Results for periods of var managed returns < -15% and weight > 1
# No. of obs. = 12
# Var and ARMA higher corr to Mkt, ARMA lower corr to Var
# With only the most extreme datapoints left, EWMA and GARCH are even negatively correlated to Mkt and var now
#  -> In most extreme negative cases, EWMA and GARCH behave in the opposite fashion of Var, which is super nice

# 5
extreme_period_5 <- print(filter(bd_times, var_managed < -10 & weight > 3))
cor(returns_m[c(2,4:7)])
cor(extreme_period_5[c(2,4:7)])
RMSD(returns_m[c(2,4:7)])
RMSD(extreme_period_5[c(2,4:7)])
# Results for periods of var managed returns < -10% and weight > 3
# No. of obs. = 11
# High weights already implies higher correlation -> not that meaningful anymore
# We take a look at root mean square deviation now
# RMSD of Var increased by far the most
# The other three increased not that much
# Relatively, ARMA increased the most of the three and EWMA the least

# 6
extreme_period_6 <- print(filter(bd_times, var_managed < -15 & weight > 3))
cor(returns_m[c(2,4:7)])
cor(extreme_period_6[c(2,4:7)])
RMSD(returns_m[c(2,4:7)])
RMSD(extreme_period_6[c(2,4:7)])
# Results for periods of var managed returns < -15% and weight > 3
# No. of obs. = 6
# RMSD of Var increased by far the most
# ARMA almost doubled, whereas EWMA and GARCH increased by roughly the same percentage



# Excel Output
bd_times$Mkt_performance <- tot_ret_m$Mkt[-1]
bd_times$var_performance <- tot_ret_m$var_managed[-1]
bd_times$ARMA_performance <- tot_ret_m$ARMA_var_managed[-1]
bd_times$EWMA_performance <- tot_ret_m$EWMA_var_managed[-1]
bd_times$GARCH_performance <- tot_ret_m$GARCH_var_managed[-1]
write_excel_csv(bd_times, "test_output.csv")
