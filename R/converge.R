onverge <- function(data, factors = NULL){
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
  return(print(X))
  cat("\nMaximum of first order:", row.names(X)[which(X[,1] == max(X[,1]))], max(X[,1]), "\n")
  cat("Maximum of total order:", row.names(X)[which(X[,2] == max(X[,2]))], max(X[,2]), "\n")
}
