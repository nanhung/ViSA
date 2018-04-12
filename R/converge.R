converge <- function(data, factors = NULL){
  p <- nrow(data$S)
  
  if (length(factors)==nrow(data$S)){
    X.labels <- factors
  } else if (is.null(factors)){
    X.labels <- paste("X", 1 : p, sep = "")
  } else {
    stop("The number of factors is not equal!")
  }
  
  X <- as.data.frame(matrix(nrow = p, ncol = 2))
  row.names(X) <- X.labels
  colnames(X) <- c("first order","total order")
  X[,1] <- data$S$`max. c.i.` - data$S$`min. c.i.`
  X[,2] <- data$T$`max. c.i.` - data$T$`min. c.i.`
  
  cat("\nModel runs: ", nrow(data$X))
  cat("\nConvergence of the sensitivity indices:\n")
  cat("\nMaximum of first order:", row.names(X)[which(X[,1] == max(X[,1]))], max(X[,1]), "\n")
  cat("Maximum of total order:", row.names(X)[which(X[,2] == max(X[,2]))], max(X[,2]), "\n")
  cat("\n")
  return(print(X))
}

plot.converge <- function(x){
  p <- 2
  pch = c(21, 24)
  n <- nrow(x)
  at <- 1 : n
  ylim <- c(0, 1)
  xlim <- c(1, n)
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", type = "n")
  axis(side = 1, at = at, labels = rownames(x))
  axis(side = 2)
  points(x$`first order`, xlim = c(1, n + 1), ylim = ylim, pch = pch[1])
  points(x$`total order`, xlim = c(1, n + 1), ylim = ylim, pch = pch[2])
  abline(0.1, 0, lty = 2)
  abline(0.05, 0, lty = 3)
  box()
  legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
}
