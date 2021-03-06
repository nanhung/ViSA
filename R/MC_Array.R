# Transfer the MCSim output
MC_Array <- function(files, start_sampling = 0){
  
  if (class(files) == "character"){
    posterior <- lapply(files, read.delim)
  } else if (class(files) == "list"){
    posterior <- files
  }
  
  n_chains <- length(posterior)
  sample_number <- dim(posterior[[1]])[1] - start_sampling
  dim <- c(sample_number, n_chains, dim(posterior[[1]])[2])
  n_iter <- dim(posterior[[1]])[1]
  n_param <- dim(posterior[[1]])[2]
  x<-array(sample_number:( n_iter * n_chains * n_param ), dim = dim)
  
  for(i in 1 : n_chains){
    x[ , i, ]<-as.matrix(posterior[[i]][(start_sampling+1):n_iter, ])  
  }
  
  dimnames(x)[[3]] <- names(posterior[[1]])
  x
}

# Show parameter names 
MC_param_names <- function(x){
  dimnames(x)[[3]] 
}


# Select the input parameters
MC_param_select <- function(x, first, last){
  pars_name <- dimnames(x)[[3]]
  start <- which(pars_name == first)
  end <- which(pars_name == last)
  pars_name[start:end]  
}

# Density plot
MC_dens<-function(x, param, dprior = NULL, dprior.arg, n = 1000, lty = 2, 
                  ref_posterior = NULL){
  par(mar=c(2,2,3,1))
  df <- as.data.frame(x[,,param])
  dens <- apply(df, 2, density)
  if (is.character(param)==TRUE){
    plot(NA, 
         xlim=range(sapply(dens, "[", "x")), 
         ylim=range(sapply(dens, "[", "y")), 
         xlab="", ylab="", main = param)
  } else {
    plot(NA, 
         xlim=range(sapply(dens, "[", "x")), 
         ylim=range(sapply(dens, "[", "y")), 
         xlab="", ylab="", main = dimnames(x)[[3]][param])
  }
  mapply(lines, dens, col=1:length(dens))
  
  if (is.null(dprior) == FALSE){
    mapply(lines, dens, col=1:length(dens))
    plot.range <- seq(par("usr")[1], par("usr")[2], length.out = n)
    prior.dens <- do.call(dprior, c(list(x=plot.range), dprior.arg))
    lines(plot.range, prior.dens, lty=lty)
  } 
  if (is.null(ref_posterior) == FALSE){
    x1<-ref_posterior[1]
    x2<-ref_posterior[2]
    y1<-par("usr")[3]
    y2<-par("usr")[4]
    polygon(x = c(x1, x1, x2, x2), 
            y = c(y1, y2, y2, y1), 
            col = rgb(1, 0, 0, 0.1), lty = 2, border = "grey")
  }
}



MC_trace<-function(x, param){
  
  X <- x[, 1, "iter"]
  df <- as.data.frame(x[,,param])

  par(mar=c(2,2,3,1))
  if (is.character(param)==TRUE){
    plot(X, df[,1],type="l", 
         ylim=range(df), 
         xlab=" ", ylab="",
         main = param)
  } else{
    plot(X, df[,1],type="l", 
         ylim=range(df), 
         xlab=" ", ylab="",
         main = dimnames(x)[[3]][param])
  }
  for (i in 2:dim(x)[2]) {
    lines(X, df[,i], col=i)
  }
}

MC_running_mean <- function (x, param) {
  par(mar = c(2, 2, 3, 1))
  X <- x[ , 1, "iter"]
  Y <- apply(x[ , , param], 2, function(x) cumsum(x)/seq_along(x))
  if (is.character(param) == TRUE) {
    plot(X , Y[,1], type = "l", ylim = range(Y), xlab = " ", 
         ylab = "", main = param)
  } else {
    plot(X , Y[,1], type = "l", ylim = range(Y), xlab = " ", 
         ylab = "", main = dimnames(x)[[3]][param])
  }
  for (i in 2:dim(x)[2]) {
    lines(X, Y[, i], col = i)
  }
}
