rfast99 <- function(factors, n, M = 4, omega = NULL,
                   q = NULL, q.arg = NULL, ...) {
  
  # factors numbers and names
  
  if (is.character(factors)) {
    X.labels <- factors
    p <- length(X.labels)
  } else {
    p <- factors
    X.labels <- paste("X", 1 : p, sep = "")
  }
  
  # quantiles
  
  if (is.null(q)) {
    q <- rep("qunif", p)
  } else if (length(q) == 1) {
    q <- rep(q, p)
  }
  if (is.null(q.arg)) {
    q.arg <- rep(list(), p)
  } else if (FALSE %in% sapply(q.arg, is.list)) { # q.arg isn't a list of lists
    q.arg <- rep(list(q.arg), p)
  }
  
  # set of frequencies
  
  if (is.null(omega)) {
    omega <- numeric(p)
    omega[1] <- floor((n - 1) / (2 * M))
    m <- floor(omega[1] / (2 * M))
    if (m >= p - 1) {
      omega[-1] <- floor(seq(from = 1, to = m, length.out = p - 1))
    } else {
      omega[-1] <- (0 : (p - 2)) %% m + 1
    }
  }
  
  # discretization of the s-space
  
  s <- 2 * pi / n * (0 : (n - 1))
  
  # transformation to get points in the x-space
  
  X <- as.data.frame(matrix(nrow = n * p, ncol = p))
  colnames(X) <- X.labels
  omega2 <- numeric(p)
  for (i in 1 : p) {
    omega2[i] <- omega[1]
    omega2[-i] <- omega[-1]
    l <- seq((i - 1) * n + 1, i * n)
    for (j in 1 : p) {
      v <- runif(1, min = 0, max = 2 * pi) # add random phase shift
      g <- 0.5 + 1 / pi * asin(sin(omega2[j] * s + v))
      X[l, j] <- do.call(q[j], c(list(p = g), q.arg[[j]]))
    }
  }
  
  # object of class "fast18"
  
  x <- list(M = M, s = s, omega = omega, X = X,
            call = match.call())
  class(x) <- "rfast99"
  
  return(x)
}

print.rfast99 <- function(x) {
  cat("")
}

tell.rfast99 <- function(x, y = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  p <- ncol(x$X)
  n <- length(x$s)
  
  V <- numeric(p)
  D1 <- numeric(p)
  Dt <- numeric(p)
  
  for (i in 1 : p) {
    l <- seq((i - 1) * n + 1, i * n)
    f <- fft(x$y[l], inverse = FALSE)
    Sp <- ( Mod(f[2 : (n / 2)]) / n )^2
    V[i] <- 2 * sum(Sp)
    D1[i] <- 2 * sum(Sp[(1 : x$M) * x$omega[1]])
    Dt[i] <- 2 * sum(Sp[1 : (x$omega[1] / 2)])
  }
  
  x$V <- V
  x$D1 <- D1
  x$Dt <- Dt
  assign(id, x, parent.frame())
}

repfast<-function(fun, factors, n = n, M = 4, q, q.arg, rep=100, ci=.95){
  
  m <- as.data.frame(matrix(ncol = factors * 2, nrow = 0))
  
  for (i in rep(n, rep)){
    x <- rfast99(factors = factors, n = n, M = M,
                 q = q, q.arg = q.arg)
    y <- fun(x$X)
    tell.rfast99(x, y)
    m <- rbind(m, c(x$D1 / x$V, 1 - x$Dt / x$V))
  }
  
  s1<-ncol(m)/2
  st<-ncol(m)/2+1
  
  p <- s1
  X.labels <- paste("X", 1 : p, sep = "")
  
  S1.Mean<- apply(m[,1:s1], 2, mean)
  S1.Lower.limit<- apply(m[,1:s1], 2, quantile, probs= c((1-ci)/2))
  S1.Upper.limit<- apply(m[,1:s1], 2, quantile, probs= c(1-(1-ci)/2))
  
  X1<-data.frame(S1.Mean, S1.Lower.limit, S1.Upper.limit, row.names = X.labels)
  names(X1)<-c("original", "min. c.i.", "max. c.i.")  
  
  St.Mean<- apply(m[,st:ncol(m)], 2, mean)
  St.Lower.limit<- apply(m[,st:ncol(m)], 2, quantile, probs= c((1-ci)/2))
  St.Upper.limit<- apply(m[,st:ncol(m)], 2, quantile, probs= c(1-(1-ci)/2))
  
  Xt<-data.frame(St.Mean, St.Lower.limit, St.Upper.limit, row.names = X.labels)
  names(Xt)<-names(X1)
  
  class(x)<-"rfast99"
  
  x$S <- X1
  x$T <- Xt
  
  cat("\n", "First order indices:", "\n")
  print(x$S)
  cat("\n", "Total order indices:", "\n")
  print(x$T)
  
  return(x)
}

