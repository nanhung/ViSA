convergeplot <- function(data, index = "T-ci"){
  
  x<-as.numeric(names(data$S_ci))
  
  if (index == "M-ci"){
    plot(x, data$S_ci[1,], type = "b", xlab = "n", 
         ylab = "Convergence index", main = "Main effect", ylim = c(0,1)) 
    abline(0.1, 0, lty = 2)
    abline(0.05, 0, lty = 3)
    for (i in 2:ncol(data$S_si)){
      lines(x, data$S_ci[i,], type = "b", col = i)
    }
  } else if (index == "T-ci"){
    plot(x, data$S_ci[1,], type = "b", xlab = "n", 
         ylab = "Convergence index", main = "Total effect", ylim = c(0,1)) 
    abline(0.1, 0, lty = 2)
    abline(0.05, 0, lty = 3)
    for (i in 2:ncol(data$T_si)){
      lines(x, data$T_ci[i,], type = "b", col = i)
    }
  } else if (index == "M-si"){
    plot(x, data$S_si[1,], type = "b", xlab = "n",
         ylab = "Sensitivity index", main = "Main effect", ylim = c(0,1)) 
    abline(0.05, 0, lty = 2)
    abline(0.01, 0, lty = 3)
    for (i in 2:ncol(data$S_si)){
      lines(x, data$S_si[i,], type = "b", col = i)
    }
  } else if (index == "T-si"){
    plot(x, data$T_si[1,], type = "b", xlab = "n",
         ylab = "Sensitivity index", main = "Total effect", ylim = c(0,1)) 
    abline(0.05, 0, lty = 2)
    abline(0.01, 0, lty = 3)
    for (i in 2:ncol(data$T_si)){
      lines(x, data$S_si[i,], type = "b", col = i)
    }
  }
  
  legend("topright", legend=rownames(data$S_si), lty = 1, col=c(1:nrow(data$S_si)))
}
