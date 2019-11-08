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

adjust_n <- 24

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