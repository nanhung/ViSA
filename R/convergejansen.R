convergejansen <- function(n, n.factors, X, fun, nboot = 100, ...){
  df <- as.data.frame(matrix(nrow = n.factors, ncol=0))
  for (n in seq(n/10, n, n/10)){
    X1 <- X(n)
    X2 <- X(n)
    x <- soboljansen(model = NULL, X1, X2, nboot) #
    y <- fun(x$X)
    t <- tell(x,y)
    c <- converge(t, factors = row.names(t$S))
    df <- cbind(df, c$`first order`)
  }
  names(df)<-seq(n/10, n, n/10)
  row.names(df)<-row.names(c)
  
  return(df)
}
