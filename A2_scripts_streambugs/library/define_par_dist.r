define.par.dist <- function(par,rel.sd=0.1,dist="LognormalTrunc",parminmax=c(0,1e30))
{
  parout <- list()
  
  for (i in 1:length(par))
  {
    pari <- c(paste(dist),as.numeric(par[i]),as.numeric(par[i])*rel.sd,parminmax)
    parout[[i]] <- pari 
  }
 
  names(parout) <- names(par)

  return(parout)  
}

# parin <- c(a=1,b=2,d=3)
# parout <-  define.par.dist(parin)
# parout <-  define.par.dist(parin,parminmax=NA)