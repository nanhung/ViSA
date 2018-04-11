convergejansen <- function(n, n.factors, fun){
  df <- as.data.frame(matrix(nrow = n.factors, ncol=0))
  for (n in seq(n/10, n, n/10)){
    X1 <- sobol.X(n)
    X2 <- sobol.X(n)
    x <- soboljansen(model = NULL, X1, X2, nboot = 100) #
    y <- fun(x$X)
    t <- tell(x,y)
    c <- converge(t)
    df <- cbind(df, c$`first order`)
  }
  names(df)<-seq(n/10, n, n/10)
  row.names(df)<-row.names(c)

  return(df)
}
