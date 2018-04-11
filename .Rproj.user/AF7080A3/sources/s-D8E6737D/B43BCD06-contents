convrge.plot<-function(df){
  x<-as.numeric(names(df))
  plot(x, df[1,], type = "b", xlab = "n", ylab = "Convergence index", ylim = c(0,1))
  abline(0.1, 0, lty = 2)
  abline(0.05, 0, lty = 3)
  for (i in 2:ncol(df)){
    lines(x, df[i,], type = "b", col = i)
  }
  legend("topright", legend=rownames(df), lty = 1, col=c(1:nrow(df)))
}
