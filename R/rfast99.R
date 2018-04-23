rfast99 <- function(factors, n, M = 4, omega = NULL,
                    q = NULL, q.arg = NULL, rep = 1, ...) {
  
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
  
    dim <- c(n * p, rep, p)
  
    a <- array((n * p):( n * p * rep * p), dim = dim)

    #X <- as.data.frame(matrix(nrow = n * p, ncol = p))
    #X <- as.data.frame(x[,1,])
    
    for (k in 1:rep){
      X <- a[,k,]
      
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
      a[,k,] <- X
    }

  # object of class "fast18"
  
  x <- list(M = M, s = s, omega = omega, a = a, rep = rep, factors = factors,
            call = match.call())
  class(x) <- "rfast99"
  
  return(x)
}

tell.rfast99 <- function(x, y = NULL, ...) {
  
  
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  p <- dim(x$a)[3]
  n <- length(x$s)
  
  V <- array(numeric(p), dim = c(p, x$rep))
  D1 <- array(numeric(p), dim = c(p, x$rep))
  Dt <- array(numeric(p), dim = c(p, x$rep))
  
  for (j in 1 : x$rep){
    for (i in 1 : p) {
      l <- seq((i - 1) * n + 1, i * n)
      f <- fft(x$y[l, j], inverse = FALSE)
      Sp <- ( Mod(f[2 : (n / 2)]) / n )^2
      V[i, j] <- 2 * sum(Sp)
      D1[i, j] <- 2 * sum(Sp[(1 : x$M) * x$omega[1]])
      Dt[i, j] <- 2 * sum(Sp[1 : (x$omega[1] / 2)])
    }
   }
  
  S_original <- apply(D1 / V, 1, mean)
  S_min_ci <- apply(D1 / V, 1, quantile, probs= c((1-0.95)/2))
  S_max_ci <- apply(D1 / V, 1, quantile, probs= c(1-(1-0.95)/2))
  S <- data.frame(S_original, S_min_ci, S_max_ci)
    
  T_original <- apply(1 - Dt / V, 1, mean)
  T_min_ci <- apply(1 - Dt / V, 1, quantile, probs= c((1-0.95)/2))
  T_max_ci <- apply(1 - Dt / V, 1, quantile, probs= c(1-(1-0.95)/2))
  T <- data.frame(T_original, T_min_ci, T_max_ci)
  
  names(S) <- names(T) <- c("original", "min. c.i.", "max. c.i.")
  
  x$V <- V
  x$D1 <- D1
  x$Dt <- Dt
  x$S <- S
  x$T <- T
  
  assign(id, x, parent.frame())
}

print.rfast99 <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", dim(x$y)[1], "\n")
    cat("\nFirst order indices::\n")
    print(x$S)
    cat("\nFirst order indices::\n")
    print(x$T)
  } else {
    cat("(empty)\n")
  }
}
