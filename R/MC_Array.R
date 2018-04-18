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

# Select the input parameters
MC_param_select <- function(x, first, last){
  pars_name <- dimnames(x)[[3]]
  start <- which(pars_name == first)
  end <- which(pars_name == last)
  pars_name[start:end]  
}
