rfast99 <- function(factors, n, M = 4, omega = NULL,
                    q = NULL, q.arg = NULL, replicate = 1, conf = 0.95, ...) {
  
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
  
  dim <- c(n * p, replicate, p)
  a <- array((n * p):( n * p * replicate * p), dim = dim)
  
  for (k in 1:replicate){
    X <- a[,k,]
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
  dimnames(a)[[3]] <- X.labels
  
  # object of class "fast18"
  
  x <- list(M = M, s = s, omega = omega, a = a, factors = factors, 
            replicate = replicate, conf = conf, call = match.call())
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
  S_se <- apply(D1 / V, 1, function(x) sqrt(var(x)/length(x)))
  S_min_ci <- apply(D1 / V, 1, quantile, probs= c((1-x$conf)/2))
  S_max_ci <- apply(D1 / V, 1, quantile, probs= c(1-(1-x$conf)/2))
  S <- data.frame(S_original, S_se, S_min_ci, S_max_ci)
  
  T_original <- apply(1 - Dt / V, 1, mean)
  T_se <- apply(1 - Dt / V, 1, function(x) sqrt(var(x)/length(x)))
  T_min_ci <- apply(1 - Dt / V, 1, quantile, probs= c((1-x$conf)/2))
  T_max_ci <- apply(1 - Dt / V, 1, quantile, probs= c(1-(1-x$conf)/2))
  T <- data.frame(T_original, T_se, T_min_ci, T_max_ci)
  
  I_original <- apply(1 - x$Dt / x$V - x$D1 / x$V, 1, mean)              
  I_se <- apply(1 - x$Dt / x$V - x$D1 / x$V, 1, function(x) sqrt(var(x)/length(x)))
  I_min_ci <- apply(1 - x$Dt / x$V - x$D1 / x$V, 1, quantile, probs= c((1-x$conf)/2))
  I_max_ci <- apply(1 - x$Dt / x$V - x$D1 / x$V, 1, quantile, probs= c(1-(1-x$conf)/2))
  I <- data.frame(I_original, I_se, I_min_ci, I_max_ci)
                
  names(I) <- names(S) <- names(T) <- c("original", "std. error", "min. c.i.", "max. c.i.")
  row.names(I) <- row.names(S) <- row.names(T) <- dimnames(x$a)[[3]]
  
  x$V <- V
  x$D1 <- D1
  x$Dt <- Dt
  x$S <- S
  x$T <- T
  x$I <- I              
  
  assign(id, x, parent.frame())
}

print.rfast99 <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (! is.null(x$y)) {
    cat("\nModel runs:", dim(x$y)[1], "\n")
    cat("\nFirst order indices::\n")
    print(x$S)
    cat("\nInteraction indices::\n")
    print(x$I)
    cat("\nTotal order indices::\n")
    print(x$T)
  } else {
    cat("(empty)\n")
  }
}

plot.rfast99 <- function(x, ylim = c(0, 1), 
                             y_col = NULL, y_dim3 = NULL, ...) {
  if (!is.null(x$y)) {
    p <- dim(x$a)[3]
    pch = c(21, 24)
    
      nodeplot(x$S, xlim = c(1, p + 1), ylim = ylim, pch = pch[1])
      nodeplot(x$T, xlim = c(1, p + 1), ylim = ylim, labels = FALSE,
               pch = pch[2], at = (1:p)+.3, add = TRUE)
      
    legend(x = "topright", legend = c("main effect", "total effect"), pch = pch)
  }
}
                
nodeplot <- function(x, xlim = NULL, ylim = NULL, labels = TRUE,
                     col = par("col"), pch = 21, bg = "white",
                     add = FALSE, at = NULL, y_col = NULL, y_dim3 = NULL, ...) {
  n <- nrow(x)
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (add) {
    par(new = TRUE)
  }
  
  # axes
  
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", type = "n", ...)
  if (class(labels) == "logical") {
    if (labels) {
      axis(side = 1, at = at, labels = rownames(x))
    } else {
      axis(side = 1, at = at, labels = FALSE, tick = FALSE)
    }
  } else if (class(labels) == "character") {
    axis(side = 1, at = at, labels = labels)
  }
  axis(side = 2)
  box()
  
  # bias
  
  if ("bias" %in% dimnames(x)[[2]]) {
    if(is.null(y_col) && is.null(y_dim3)){
      xx <- x[["original"]] - x[["bias"]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      xx <- x[, "original", y_col] - x[, "bias", y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      xx <- x[, "original", y_col, y_dim3] - x[, "bias", y_col, y_dim3]
    }
  } else {
    if(is.null(y_col) && is.null(y_dim3)){
      xx <- x[["original"]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      xx <- x[, y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      xx <- x[, y_col, y_dim3]
    }
  }
  
  # confidence intervals
  
  if (("min. c.i." %in% dimnames(x)[[2]]) & "max. c.i." %in% dimnames(x)[[2]]) {
    if(is.null(y_col) && is.null(y_dim3)){
      min_ci <- x[["min. c.i."]]
      max_ci <- x[["max. c.i."]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      min_ci <- x[, "min. c.i.", y_col]
      max_ci <- x[, "max. c.i.", y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      min_ci <- x[, "min. c.i.", y_col, y_dim3]
      max_ci <- x[, "max. c.i.", y_col, y_dim3]
    }
    for (i in 1 : n) {
      lines(c(at[i], at[i]), c(min_ci[i], max_ci[i]),
            col = col)
    }
  }
  
  # points
  
  points(at, xx, col = col, pch = pch, bg = bg)
}
                
