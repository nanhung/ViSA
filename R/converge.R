converge <- function(data){
  
  diff <- as.data.frame(matrix(nrow = nrow(data$S), ncol = 2))
  row.names(diff) <- row.names(data$S)
  colnames(diff) <- c("first order","total order")
  diff[,1] <- data$S$`max. c.i.` - data$S$`min. c.i.`
  diff[,2] <- data$T$`max. c.i.` - data$T$`min. c.i.`
  
  data$diff <- diff
  data$S$diff <- data$S$`max. c.i.` - data$S$`min. c.i.`
  data$T$diff <- data$T$`max. c.i.` - data$T$`min. c.i.`
  
  p <- 2
  pch = c(21, 24)
  n <- nrow(data$S)
  at <- 1 : n
  ylim <- c(0, 1)
  xlim <- c(1, n)
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", type = "n")
  axis(side = 1, at = at, labels = rownames(data$S))
  axis(side = 2)
  points(data$S$diff, xlim = c(1, n + 1), ylim = ylim, pch = pch[1])
  points(data$T$diff, xlim = c(1, n + 1), ylim = ylim, pch = pch[2])
  abline(0.1, 0, lty = 2)
  abline(0.05, 0, lty = 3)
  box()
  
  legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  cat("\nConvergence of the sensitivity indices:\n")
  cat("\nMaximum of first order:", row.names(data$diff)[which(data$diff[,1] == max(data$diff[,1]))], max(data$diff[,1]), "\n")
  cat("Maximum of total order:", row.names(data$diff)[which(data$diff[,2] == max(data$diff[,2]))], max(data$diff[,2]), "\n")
  cat("\n")
  print(data$diff)
}
