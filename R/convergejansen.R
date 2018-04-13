convergejansen <- function(n, n.factors, X, fun, nboot = 100){
  
  S_si <- as.data.frame(matrix(nrow = n.factors, ncol=0))
  T_ci <- S_ci <- T_si <- S_si
  
  for (n in seq(n/10, n, n/10)){
    X1 <- X(n)
    X2 <- X(n)
    x <- soboljansen(model = NULL, X1, X2, nboot) #
    y <- fun(x$X)
    tell(x,y)
    
    x$S$diff <- x$S$`max. c.i.` - x$S$`min. c.i.`
    x$T$diff <- x$T$`max. c.i.` - x$T$`min. c.i.`
    
    S_si <- cbind(S_si, x$S$original)
    T_si <- cbind(T_si, x$T$original)
    S_ci <- cbind(S_ci, x$S$diff)
    T_ci <- cbind(T_ci, x$T$diff)
  }
  names(S_si)<-names(T_si)<-names(S_ci)<-names(T_ci)<-seq(n/10, n, n/10)
  row.names(S_si)<-row.names(T_si)<-row.names(S_ci)<-row.names(T_ci)<-row.names(x$S)
  
  return(list(S_si = S_si, T_si = T_si, S_ci = S_ci, T_ci = T_ci))
}
