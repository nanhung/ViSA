# Transfer the MCSim output
MC_Array <- function(files, start_sampling = 0){
  # files = list.files(pattern="*.mcmc.out")
  files <- files
  posterior <- lapply(files, read.delim)
  n.chains <- length(posterior)
  start_sampling <- start_sampling
  sample_number <- dim(posterior[[1]])[1] - start_sampling
  dim <- c(sample_number,n.chains,dim(posterior[[1]])[2])
  n_row <- dim(posterior[[1]])[1]
  n_iter <- dim(posterior[[1]])[2]
  x<-array(sample_number:(n_row*n.chains*n_iter), dim = dim)
  for(i in 1 : n.chains){
    x[,i,]<-as.matrix(posterior[[i]][(start_sampling+1):n_row,])  
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
MC_dens<-function(x, param, dprior = NULL, dprior.arg, n = 1000, lty = 2){
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
}

MC_trace<-function(x, param){
  df <- as.data.frame(x[,,param])
  
  par(mar=c(2,2,3,1))
  if (is.character(param)==TRUE){
    plot(df[,1],type="l", 
         ylim=range(df), 
         xlab=" ", ylab="",
         main = param)
  } else{
    plot(df[,1],type="l", 
         ylim=range(df), 
         xlab=" ", ylab="",
         main = dimnames(x)[[3]][param])
  }
  for (i in 2:dim(x)[2]) {
    lines(df[,i],col=i)
  }
}
