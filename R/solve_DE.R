solve_DE <- function(fun = NULL, x, times, parameters, initState, 
                      func, jacfunc, initfunc, nout = 1, outnames){
  n <- length(x$s)  
  factors <- ifelse (class(x$factors) == "character", length(x$factors), factors) 
  replicate <- x$rep
  out <- length(times)
  y <- array(1:replicate*n*factors*out, 
             dim = c(n * factors, replicate, out))
  
  if (is.null(fun)){
    for (k in 1 : dim(y)[3]) { #outputs
      
      # Specific time or variable
      inputs = c(0, times[k])
      
      for (i in 1 : dim(y)[2]) { # replicate
        for (j in 1 : dim(y)[1]) { # Model evaluation
          for (p in 1 : factors) { # input individual factors
            parameters[p] <- x$a[j,i,p]
          }
          
          # Integrate
          tmp <- ode(initState, inputs, func = func, parms = parameters, 
                     jacfunc = jacfunc, dllname = mName, 
                     initfunc = initfunc, nout = nout, outnames = outnames)
          y[j,i,k] <- tmp[2, outnames]
        }
      }
    }
  } else(
    for (i in 1 : dim(y)[2]) { # Replicate
      for (j in 1 : dim(y)[1]) { # Model evaluation
        parameters <- x$a[j,i,]
        tmp <- fun(parameters, t)
        
        for (k in 1 : dim(y)[3]) { # Output time
          y[j,i,k] <- tmp[k]
        }
      }
    }
  )
  dimnames(y)[[3]]<-times
  return(y)
}
